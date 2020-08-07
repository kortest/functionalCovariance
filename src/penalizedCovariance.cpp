#include <Rcpp.h>
#include <RcppEigen.h>
#include <RcppGSL.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_matrix.h>
#include <math.h>
#include "functionalCovariance_types.h"

using namespace Rcpp;

//typedef Eigen::MappedSparseMatrix<double> MSpMat;
//typedef Eigen::SparseMatrix<double, Eigen::RowMajor> SpMat;
//typedef Eigen::SparseMatrix<double, Eigen::ColMajor> SpMatC;
//typedef Eigen::Triplet<double, int> T;

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppGSL)]]


/*
 Given a numeric vector {y_i}, returns a vector with entries y_i * y_j
 with either i < j or i <= j
 Args: y: numeric vector
       mixedOnly: if TRUE, exclude products of the form y_i * y_i
 
 Returns: A numeric Vector with the product entries
 */

// [[Rcpp::export]]
NumericVector mixedMomentsC(NumericVector y, bool mixedOnly = false) {
  int n = y.size();
  // Booleans in Rcpp are translated to integers
  NumericVector out((n+1-mixedOnly) * (n - mixedOnly)  / 2);
  int counter = 0;
  for(int i = 0; i < n-mixedOnly; i++){
    for(int j = i+mixedOnly; j < n; j++){
      out[counter] = y[i]*y[j];
      counter++;
    }
  }
  return out;
}

/*
 Compute the b-spline design Matrix. 5 times faster than fda's 
 eval.basis (no error checking though) but most importantly, 
 allows to do everything in C without any conversions.
 
 Args: lower: lower interval limit for splines
       upper: lower interval limit for splines
       order: order of the bspline basis
       nbreaks: Number of internal knots
       x: Vector where to evaluate the basis
 
 Returns: The spline design matrix
 */

// [[Rcpp::export]]
NumericMatrix splineMatrixC(const double lower, const double upper, const int order, const int nbreaks, const NumericVector x){
  // The gsl library creates one object that stores 
  // all necessary b-spline basis attributes. This is this
  // workspace. 
  gsl_bspline_workspace *bw;
  gsl_vector *B;
  
  //Initialize...
  const int n = x.size();
  const int nsplines = nbreaks + order - 2; 
  bw = gsl_bspline_alloc(order, nbreaks);
  B = gsl_vector_alloc(nsplines);
  // Could potentially already be a sparse matrix
  // However, if other basis produce dense matrices,
  // then it would be problematic
  NumericMatrix X = NumericMatrix(n, nsplines);
  
  // This (and the function generally)
  // should probably be made more flexible in the future
  gsl_bspline_knots_uniform(lower, upper, bw);
  
  /* construct the fit matrix X */
  for (int i = 0; i < n; ++i){
    
    /* compute B_j(xi) for all j */
    gsl_bspline_eval(x[i], B, bw);
    
    /* fill in row i of X */
    for (int j = 0; j < nsplines; ++j){
      X(i,j) = gsl_vector_get(B, j);
    }
  }
  return X;
}

/*
 Compute the b-spline penalty design Matrix assuming that both
 b-spline basis used are the same. Takes as input the penalty
 matrix for one basis. Ideally, this would be computed in C++
 as well.
 
 Args: X: penalty matrix for the basis
       nzero: estimate for how many non zero values there are
 
 Returns: The combined spline penalty matrix
 */

// [[Rcpp::export]]
SpMat penaltyMatrixC(const NumericMatrix X, const int nzero) {
  
  const int n = X.ncol();
  
  const int new_rows = n * (n + 1) / 2; 
  
  std::vector<int> p(new_rows+1);
  std::vector<int> jind(new_rows * nzero);
  std::vector<double> x(new_rows * nzero);
  
  int rsc = -1;
  int csc = -1;
  int lc = 0;
  int limit = 0;
  int obsc = 0;
  int cco;
  
  for (int i = 0; i < new_rows; i++){
    if (i == limit){
      lc = 0;
      rsc++;
      limit = limit + n - rsc;
      csc = rsc; 
      x[obsc] = X(rsc, rsc);
    } else {
      x[obsc] = X(rsc, rsc) + X(rsc + lc, rsc + lc);
    }
    
    jind[obsc] = i;
    obsc++;
    
    for (int j = rsc + lc + 1; j < n; j++){
      if(X(rsc + lc, j)==0){
        continue;
      }
      x[obsc] = X(rsc + lc, j);
      jind[obsc] = i + j - rsc - lc;
      obsc++;
    }
    
    if (lc > 0){
      cco = limit;
      for (int j = rsc + 1; j < rsc + lc; j++){
        if (X(rsc, j)==0){
          cco += n - j;
        } else {
          x[obsc] = X(rsc, j);
          jind[obsc] = cco + lc + rsc - j;
          obsc++;
          cco += n - j;
        }
      }
      for (int j = 0; j < n - lc - rsc; j++){
        if (X(rsc, rsc + lc + j)==0){
        } else {
          x[obsc] = X(rsc, rsc + lc + j);
          jind[obsc] = cco + j;
          obsc++;
        }
      }
    }
    lc++;
    p[i+1] = obsc;
  }
  jind.resize(obsc);
  x.resize(obsc);
  Eigen::Map<Eigen::SparseMatrix<double, Eigen::RowMajor>> spMap(new_rows, new_rows, obsc, p.data(), jind.data(), x.data(), 0);
  Eigen::SparseMatrix<double, Eigen::RowMajor> res= spMap.eval();
  return res.selfadjointView<Eigen::Upper>();
}

/*
 Given data and a basis, the goal is to estimate coefficients
 of the covariance function in the product basis. To achieve this
 a design matrix needs to be built that corresponds to the 
 evaluation of the product basis at possible combinations of the 
 observation indexes 
 
 This function takes a sparse Matrix corresponding to basis
 evaluations and computes this design matrix. Since the covariance
 function is necessarily symmetric, the coefficients will also be 
 symmetric and thus, only a smaller design matrix is needed where
 entries corresponding to equal coefficients (due to symmetry) are
 summed.
 
 Data in a sparse (row) matrix in R is stored as three vectors, 
 a vector j indicating the columns of non-zero values, a
 vector x containing the corresponding non-zeros values and 
 a vector p indicating where new rows start in the vectors 
 j and x respectively. All the indexes are 0 based.
 For example: Let M be some matrix and a = M[10, 12] is the third
 non-zero entry in the tenth row and there were 30 non-zero 
 entries in the first 9 rows then 
 p[9] = 30 (Non-zero elements in the 10th row start at index 30)
 j[32] = 12 (The third non-zero element in row 10 is in column 12)
 x[32] = a (The values corresponding to the index [10,12])
 
 Args: dM: The design matrix in the original basis
       nnzero: An estimate of the nnzero entries per row in the new matrix,
               needs to be at least as large as the true value on average.
       mixedOnly: if TRUE, should exclude variances for estimation
 
 Returns: The new sparse design matrix. This will be in row major format
 
 Note, this assumes the basis evaluation matrix has been evaluated at an increasing
 sequence of indexes and that the basis behaves like a spline basis in the sense that
 if i < j, then the first non-zero column of the i-th row is smaller than or equal to
 the first non-zero column of the the j-th row
 */
// [[Rcpp::export]]
const SpMat singleDesignMatrixC(const SpMat dM, const int nnzero, const bool mixedOnly = false) {
  // The switch here is due to the fact that for this purpose
  // the column wise storage is not as good as row wise storage
  // so that I consider the transpose.
  const int n = dM.rows();
  const int ncol = dM.cols();
  
  const int *P = dM.outerIndexPtr();
  const int *J= dM.innerIndexPtr();
  const double *X = dM.valuePtr();
  
  // Overlap counts how many columns the two rows currently considered have in common
  int overlap;
  // Technical integer needed to correctly determine the start of a loop
  int start;
  int test;
  
  // Counters to keep track of where in the loop we are
  int linecounter = 0;
  int indcounter = 0;
  int valcounter = 0;
  
  // Integers counting how many non-zeros entries there are in the two rows currently considerd 
  int nzero1 = 0;
  int nzero2 = 0;
  
  // Compute the number of rows in the design matrix
  const int new_rows = (n+1-mixedOnly) * (n-mixedOnly) / 2;
  const int new_cols = (ncol+1) * ncol / 2;
  
  // Initialize the vectors to construct the new arma::sp_mat return values
  std::vector<int> p(new_rows+1);
  std::vector<int> jind(new_rows * nnzero);
  std::vector<double> x(new_rows * nnzero);
  // The first pointer must always 0
  p[0] = 0;
  
  // Start with the i-th row of the original matrix...
  for(int i = 0; i < n-mixedOnly; i++) {
    // Compute the number of non zero values in this row
    nzero1 = P[i+1] - P[i];
    // If this number is zero, then the n - i rows of potential 
    // products are all 0. But need to keep track of lines
    if (nzero1 == 0){
      for(int j = i+1; j<= n; j++){
        p[linecounter+1] = p[linecounter];
        linecounter++;
      }
      continue;
    }
    // For all rows after the i-th row...
    for(int j = i+mixedOnly; j < n; j++){
      // If row has only zeros continue, but add line
      nzero2 = P[j+1] - P[j];
      if (nzero2 == 0){
        p[linecounter+1] = p[linecounter];
        linecounter++;
        continue;
      }
      // Compute the overlap of columns non-zero indexes. If there is overlap, not all possible
      // products produce non-zeros entries, but some products contribute
      // to the same non-zero entry
      test = J[P[i+1] -1] - J[P[j]] + 1;
      overlap = std::max(0, test);
      p[linecounter+1] = p[linecounter] + (nzero1 * nzero2) - (overlap * (overlap - 1) / 2); 
      indcounter = 0;
      // The if here is to avoid unecessary computation in case there is no overlap
      // which should be the clear majority of cases
      if (overlap <= 1){
        // Iterate through the non-zeros elements of the two rows and compute
        // the new indeces and values. Since there is no overlap, no columns
        // need to be summed
        for(int k = 0; k < nzero1; k++){
          for(int l = 0; l < nzero2; l++){
            int ind = p[linecounter] + indcounter;
            jind[ind] = J[P[i] + k] * ncol - (J[P[i] + k] * (J[P[i] + k]-1) / 2) + J[P[j]+l] - J[P[i] + k];
            x[ind] = X[P[i] + k] * X[P[j] + l];
            indcounter++;
            valcounter++;
          }
        }
      } else {
        // Same as above, but now take care of adding columns since there is overlap
        for (int k = 0;  k < nzero1; k++){
          start = std::max(0, overlap - nzero1 + k);
          for(int l = start; l < nzero2; l++){
            int ind = p[linecounter] + indcounter;
            jind[ind] = J[P[i] + k] * ncol - (J[P[i] + k] * (J[P[i] + k]-1) / 2) + J[P[j]+l] - J[P[i] + k];
            if (l < overlap && l > start && k >= nzero1 - overlap){
              x[ind] = X[P[i] + k] * X[P[j] + l] + X[P[i] + l + nzero1 - overlap] * X[P[j] + k - nzero1 + overlap];
            } else {
              x[ind] = X[P[i] + k] * X[P[j] + l];
            }
            indcounter++;
            valcounter++;
          }
        }
      }
      // Next pair of rows...
      linecounter++;
    }
  }
  jind.resize(valcounter);
  x.resize(valcounter);
  Eigen::Map<Eigen::SparseMatrix<double, Eigen::RowMajor>> spMap(new_rows, new_cols, valcounter, p.data(), jind.data(), x.data(), 0);
  Eigen::SparseMatrix<double, Eigen::RowMajor> res= spMap.eval();
  return res;
}

/*
 Gives the whole design matrix
 
 Args: indexes: index points for the observation
       values: values observed at the index points
       lower: lower interval bound for the splines
       upper: upper interval bound for the splines
       nbasis: how man basis functions?
       order: order of splines
       mixedOnly: Whether to include the variances
 
 Returns: The t(Phi) %*% Phi design matrix as well as
          T(Phi) %*% Y where Y are the mixed moments 
          estimates.
 */
// [[Rcpp::export]]
List fullDesignC(const List indexes,
                 const List values,
                 const int lower,
                 const int upper,
                 const int nbasis,
                 const int order,
                 const Eigen::VectorXd weights,
                 const bool mixedOnly=false){
  // Initialize...
  // I think I haven't fully understood Eigen but sometimes
  // it seems to require intermediate assignments
  NumericMatrix spline_mat;
  Eigen::MatrixXd intm;
  SpMat Phi;
  SpMat res2;
  SpMat spline_sparse;
  Eigen::VectorXd res2Y;
  
  size_t n_funcs = indexes.size();
  
  // Initialize the return values and set them 0
  int res_size = (nbasis+1) * nbasis / 2;
  SpMat res = SpMat(res_size, res_size);
  Eigen::VectorXd resY = Eigen::VectorXd(res_size);
  res.setZero();
  resY.setZero();
  
  // Loop through all the profiles...
  // This somehow begs for paralleled implementation
  // but this really depends on the whole structure
  for (size_t i = 0; i < n_funcs; i++){
    spline_mat = splineMatrixC(lower, upper, order, nbasis - 2, indexes[i]);
    intm = Rcpp::as<Eigen::MatrixXd>(spline_mat);
    spline_sparse  = intm.sparseView();
    // The order squared will always be the least upper bound for Bsplines
    Phi = singleDesignMatrixC(spline_sparse, order*order, mixedOnly);
    // Calculate the matrix to solve for the coefficients
    // Time-wise by FAR the most expensive line
    // Eigen gives 7 seconds over Armadillo but there might sill be potential.
    // This would be a next step...
    res2 = res.selfadjointView<Eigen::Lower>().rankUpdate(Phi.adjoint(), weights[i]);
    res = res2;
    // Compute the transformation on the response vectors
    res2Y = Phi.adjoint() * Rcpp::as<Eigen::VectorXd>(mixedMomentsC(values[i], mixedOnly));
    resY = resY + weights[i] * res2Y;
  }
  return Rcpp::List::create(res, resY);
}


/*
 Computes the solution at once for a given 
 penalty
 
 Args: fullDesign: List containing t(Phi) %*% Phi
                   and t(Phi) %*% Y
       penalty: penalty Matrix for the dual basis
       lambda: penalty parameter
 
 Returns: The upper (lower) triangular part of the
          coefficient matrix as a vector
 */
// [[Rcpp::export]]
Eigen::VectorXd estimateCoefficientsC(const List fullDesign,
                                      const SpMat penalty,
                                      const double lambda){

  SpMat B;
  Eigen::VectorXd beta;
  
  // Get the design
  SpMat designX = fullDesign[0];
  Eigen::VectorXd designY = fullDesign[1];
  
  B = lambda * penalty +  designX;
  // Compute the cholesky decomposition to solve
  // There maybe some gains using the permutation matrices,
  // however, it would then later entail another matrix multiplication. 
  // Haven't checked which is faster. The inverse without permutation is
  // only slightly less sparse than the one with permutation.
  Eigen::SimplicialLLT<Eigen::SparseMatrix<double>, Eigen::Lower, Eigen::NaturalOrdering<int>> cholesky;
  cholesky.compute(B);
  beta = cholesky.solve(designY);
  return beta;
}

/*
 * Directly copied from R package sparseinv 
 * with some modifications for types of vectors
 */

NumericVector sparseinv2(int n, const int *Lp, const int *Li, const double *Lx, Eigen::VectorXd d, const int *Up, const int *Uj, const double *Ux, const int *Zp, const int *Zi) {
  
  double ljk, zkj ;
  int j, i, k, p, znz, pdiag, up, zp, flops = n ;
  double *z = (double *) calloc(n,sizeof(double));
  int *Zdiagp = (int *) malloc((n)*sizeof(int));
  int *Lmunch = (int *) malloc((n)*sizeof(int));
  
  /* ---------------------------------------------------------------------- */
  /* initializations */
  /* ---------------------------------------------------------------------- */
  
  int counter = 0;
  
  /* clear the numerical values of Z */
  int lengthZx = Zp[n] ;
  NumericVector Zx(lengthZx);
  znz = Zp[n] ;
  for (p = 0 ; p < znz ; p++)
  {
    Zx[p] = 0 ;
  }
  
  /* find the diagonal of Z and initialize it */
  for (j = 0 ; j < n ; j++)
  {
    pdiag = -1 ;
    for (p = Zp [j] ; p < Zp [j+1] && pdiag == -1 ; p++)
    {
      if (Zi [p] == j)
      {
        pdiag = p ;
        Zx [p] = 1 / d [j] ;
      }
    }
    Zdiagp [j] = pdiag ;
    if (pdiag == -1) return (-1) ;  /* Z must have a zero-free diagonal */
  }
  
  /* Lmunch [k] points to the last entry in column k of L */
  for (k = 0 ; k < n ; k++)
  {
    Lmunch [k] = Lp [k+1] - 1 ;
  }
  /* ---------------------------------------------------------------------- */
  /* compute the sparse inverse subset */
  /* ---------------------------------------------------------------------- */
  
  for (j = (n)-1 ; j >= 0 ; j--)
  {
    /* ------------------------------------------------------------------ */
    /* scatter Z (:,j) into z workspace */
    /* ------------------------------------------------------------------ */
    
    /* only the lower triangular part is needed, since the upper triangular
     part is all zero */
    for (p = Zdiagp [j] ; p < Zp [j+1] ; p++)
    {
      z [Zi [p]] = Zx [p] ;
    }
    
    /* ------------------------------------------------------------------ */
    /* compute the strictly upper triangular part of Z (:,j) */
    /* ------------------------------------------------------------------ */
    
    /* for k = (j-1):-1:1 but only for the entries Z(k,j) */
    for (p = Zdiagp [j]-1 ; p >= Zp [j] ; p--)
    {
      /* Z (k,j) = - U (k,k+1:n) * Z (k+1:n,j) */
      k = Zi [p] ;
      zkj = 0 ;
      flops += (Up [k+1] - Up [k]) ;
      for (up = Up [k] ; up < Up [k+1] ; up++)
      {
        /* skip the diagonal of U, if present */
        i = Uj [up] ;
        if (i > k)
        {
          zkj -= Ux [up] * z [i] ;
        }
      }
      z [k] = zkj ;
    }
    
    /* ------------------------------------------------------------------ */
    /* left-looking update to lower triangular part of Z */
    /* ------------------------------------------------------------------ */
    /* for k = (j-1):-1:1 but only for the entries Z(k,j) */
    for (p = Zdiagp [j]-1 ; p >= Zp [j] ; p--)
    {
      
      k = Zi [p] ;
      
      /* ljk = L (j,k) */
      if (Lmunch [k] < Lp [k] || Li [Lmunch [k]] != j)
      {
        /* L (j,k) is zero, so there is no work to do */
        continue ;
      }
      ljk = Lx [Lmunch [k]--] ;
      
      /* Z (k+1:n,k) = Z (k+1:n,k) - Z (k+1:n,j) * L (j,k) */
      flops += (Zp [k+1] - Zdiagp [k]) ;
      for (zp = Zdiagp [k] ; zp < Zp [k+1] ; zp++)
      {
        Zx [zp] -= z [Zi [zp]] * ljk ;
      }
    }
    
    /* ------------------------------------------------------------------ */
    /* gather Z (:,j) back from z workspace */
    /* ------------------------------------------------------------------ */
    
    for (p = Zp [j] ; p < Zp [j+1] ; p++)
    {
      i = Zi [p] ;
      Zx [p] = z [i] ;
      z [i] = 0 ;
    }
  }
  
  free(z);
  free(Zdiagp);
  free(Lmunch);
  
  return Zx;
}

/*
 * Computes the inverse using the Takahashi Davis formula
 * Essentially, this a translation (together with the
 * the function above) from their package. It is a bit messy
 * because Eigen is hard to deal with at times. It should 
 * probably be cleaned up at some point.
 * 
 * Args: L: The L matrix from the Cholesky decomposition
 * 
 * Returns: The Takahashi-Davis sparse inverse
 */

// [[Rcpp::export]]
SpMatC takahashiDavisC(const SpMat L){
  
  const int nrows = L.rows();
  SpMat Lnew;
  SpMatC LnewC;
  SpMat D; 
  Eigen::VectorXd d;
  Eigen::VectorXd di;
  Eigen::VectorXi dp(nrows);
  Eigen::VectorXi jj;
  Eigen::VectorXi ii;
  const int *p;
  int *Li;
  int dpsum;
  int nzero;
  int LnewNzero;
  std::vector<T> tripletList;
  int i;
  
  // The following lines basically follow the sparseinv package and
  // prepare the call to the above function
  // Extract the L Matrix
  d = L.diagonal();
  di = d.array().inverse();
  D = SpMat( di.asDiagonal() );
  Lnew = L * D;
  d = d.array().square();
  Lnew.diagonal().setZero();
  Lnew = Lnew.pruned();
  LnewC = SpMatC(Lnew);
  p = LnewC.outerIndexPtr();
  for (i = 0; i < nrows; i++){
    dp[i] = p[i+1] - p[i];
  }
  dpsum = dp.sum();
  
  // These vectors ii and jj will be needed in the next step
  LnewNzero = LnewC.nonZeros();
  jj = Eigen::VectorXi(dpsum);
  ii = Eigen::VectorXi(LnewNzero);
  Li = LnewC.innerIndexPtr();
  Eigen::VectorXi vec_joined(LnewNzero + dpsum);
  Eigen::VectorXi vec_joined2(LnewNzero + dpsum);
  
  int counter = 0;
  for (int j = 0; j < nrows; j++){
    for (i = 0; i < dp[j]; i++){
      jj[counter + i] = j;
    }
    counter += dp[j];
  }
  
  for(i = 0; i < LnewNzero; i++){
    ii[i] = Li[i];
  }
  
  // Join ii and jj as in the R package
  vec_joined << ii, jj;
  vec_joined2 << jj, ii;
  
  // Basically, here the index pattern for the inverse is set up,
  // thus, nzero counts the number of non-zeros in the inverse.
  nzero = dpsum + LnewNzero + nrows;
  
  // Should memory allocation always happen at the beginning?
  tripletList.reserve(nzero);
  
  // Create the triplets to create a sparse matrix
  // The values do not matter, really, this is just
  // a convenient way to get the correct column and
  // row pointers from the un-ordered arrays above
  for(i = 0; i < LnewNzero + dpsum; i++){
    tripletList.push_back(T(vec_joined[i], vec_joined2[i], 1.0));
  }
  for(i = 0; i < nrows; i++){
    tripletList.push_back(T(i, i, 1.0));
  }
  SpMatC Z(nrows, nrows);
  Z.setFromTriplets(tripletList.begin(), tripletList.end());
  
  // Now compute the values of the inverse using the above formula
  Eigen::VectorXd Zx = Rcpp::as<Eigen::VectorXd>(sparseinv2(nrows, LnewC.outerIndexPtr(), LnewC.innerIndexPtr(), LnewC.valuePtr(), d, LnewC.outerIndexPtr(), LnewC.innerIndexPtr(), LnewC.valuePtr(), Z.outerIndexPtr(), Z.innerIndexPtr()));
  
  // Now construct the sparseMatrix to be returned, unfortunately col major....
  Eigen::Map<Eigen::SparseMatrix<double, Eigen::ColMajor>> spMap(nrows, nrows, Z.nonZeros(), Z.outerIndexPtr(), Z.innerIndexPtr(), Zx.data(), 0);
  Eigen::SparseMatrix<double, Eigen::ColMajor> res= spMap.eval();
  return res;
}


/*
 * Computes the prediction in a fast way. May not be necessary
 * but is faster than constructing the single design and then
 * multiplying it by the coefficient vector
 * 
 * Args: dM: Spline design matrix
 *       beta: coefficient vector
 *       mixedOnly: should variances be included?
 * 
 * Returns: The vector of predictions
 */
// [[Rcpp::export]]
Eigen::VectorXd computePredC(const SpMat dM, const Eigen::VectorXd beta, const bool mixedOnly = false) {
  // The switch here is due to the fact that for this purpose
  // the column wise storage is not as good as row wise storage
  // so that I consider the transpose.
  const int n = dM.rows();
  const int ncol = dM.cols();
  
  double sum = 0;
  Eigen::VectorXd result(n*(n-1)/2);
  
  const int *P = dM.outerIndexPtr();
  const int *J= dM.innerIndexPtr();
  const double *X = dM.valuePtr();
  
  // Overlap counts how many columns the two rows currently considered have in common
  int overlap;
  // Technical integer needed to correctly determine the start of a loop
  int start;
  int test;
  
  // Counters to keep track of where in the loop we are
  int linecounter = 0;
  
  int testcounter = 0;
  // Integers counting how many non-zeros entries there are in the two rows currently considerd 
  int nzero1 = 0;
  int nzero2 = 0;
  
  // Start with the i-th row of the original matrix...
  for(int i = 0; i < n-mixedOnly; i++) {
    // Compute the number of non zero values in this row
    nzero1 = P[i+1] - P[i];
    // If this number is zero, then the n - i rows of potential 
    // products are all 0. But need to keep track of lines
    if (nzero1 == 0){
      for(int j = i+1; j<= n; j++){
        linecounter++;
      }
      continue;
    }
    // For all rows after the i-th row...
    for(int j = i+mixedOnly; j < n; j++){
      // If row has only zeros continue, but add line
      nzero2 = P[j+1] - P[j];
      if (nzero2 == 0){
        linecounter++;
        continue;
      }
      // Compute the overlap of columns non-zero indexes. If there is overlap, not all possible
      // products produce non-zeros entries, but some products contribute
      // to the same non-zero entry
      test = J[P[i+1] -1] - J[P[j]] + 1;
      overlap = std::max(0, test); 
      // The if here is to avoid unecessary computation in case there is no overlap
      // which should be the clear majority of cases
      if (overlap <= 1){
        // Iterate through the non-zeros elements of the two rows and compute
        // the new indeces and values. Since there is no overlap, no columns
        // need to be summed
        for(int k = 0; k < nzero1; k++){
          for(int l = 0; l < nzero2; l++){
            sum += beta[J[P[i] + k] * ncol - (J[P[i] + k] * (J[P[i] + k]-1) / 2) + J[P[j]+l] - J[P[i] + k]] * X[P[i] + k] * X[P[j] + l];
          }
        }
      } else {
        // Same as above, but now take care of adding columns since there is overlap
        for (int k = 0;  k < nzero1; k++){
          start = std::max(0, overlap - nzero1 + k);
          for(int l = start; l < nzero2; l++){
            if (l < overlap && l > start && k >= nzero1 - overlap){
              sum += beta[J[P[i] + k] * ncol - (J[P[i] + k] * (J[P[i] + k]-1) / 2) + J[P[j]+l] - J[P[i] + k]] * (X[P[i] + k] * X[P[j] + l] + X[P[i] + l + nzero1 - overlap] * X[P[j] + k - nzero1 + overlap]);
            } else {
              sum += beta[J[P[i] + k] * ncol - (J[P[i] + k] * (J[P[i] + k]-1) / 2) + J[P[j]+l] - J[P[i] + k]] * X[P[i] + k] * X[P[j] + l];
            }
          }
        }
      }
      result[linecounter] = sum;
      sum = 0;
      // Next pair of rows...
      linecounter++;
    }
  }
  return result;
}

/*
 * Computes the Sum of Squared errors. The disadvantage of never
 * having the whole design matrix as one is that this will 
 * always have to reconstruct it. However, this should be 
 * relatively quick since the predictions are computed at the
 * same time
 * 
 * Args: indexes: index points for the observation
 *       values: values observed at the index points
 *       lower: lower interval bound for the splines
 *       upper: upper interval bound for the splines
 *       nbasis: how many basis functions should be used
 *       order: order of the splines
 *       beta: vector of coefficients
 *       mixedOnly: should variances be included?
 * 
 * Returns: The sum of squared errors
 */
// [[Rcpp::export]]
double computeSSEC(const List indexes,
                   const List values,
                   const int lower,
                   const int upper,
                   const int nbasis,
                   const int order,
                   const Eigen::VectorXd beta,
                   const Eigen::VectorXd weights,
                   const bool mixedOnly=false){
  Eigen::VectorXd preds;
  Eigen::MatrixXd intm;
  NumericMatrix spline_mat;
  SpMat spline_sparse;
  SpMat Phi;
  size_t n_funcs = indexes.size();
  double res = 0; 
  
  for (size_t i = 0; i < n_funcs; i++){
    spline_mat = splineMatrixC(lower, upper, order, nbasis - 2, indexes[i]);
    intm = Rcpp::as<Eigen::MatrixXd>(spline_mat);
    spline_sparse  = intm.sparseView();
    preds = computePredC(spline_sparse, beta, mixedOnly);
    intm = preds - Rcpp::as<Eigen::VectorXd>(mixedMomentsC(values[i], mixedOnly));
    res += weights[i] * intm.array().square().sum();
  }
  
  return res;
}


/*
 * Function to be used in optimize to compute the
 * optimal penalty coefficient using the generalized
 * cross validation method
 * 
 * Args: lambda: penalty parameter
 *       indexes: index points for the observation
 *       values: values observed at the index points
 *       lower_int: lower interval bound for the splines
 *       upper_int: upper interval bound for the splines
 *       nbasis: how many basis functions should be used
 *       order: order of the splines
 *       designX: t(Phi) %*% Phi
 *       designY: t(Phi) %*% Y
 *       penalty: joint penalty matrix
 *       multiplier: change scale of penalty parameter
 *       mixedOnly: should variances be included?
 * 
 * Returns: The gcv score
 */
// [[Rcpp::export]]
double generalizedCrossValidationC(double lambda,
                                   const List indexes,
                                   const List values,
                                   const int lower_int,
                                   const int upper_int,
                                   const int nbasis,
                                   const int order,
                                   const SpMat designX,
                                   const Eigen::VectorXd designY,
                                   const SpMat penalty,
                                   const double multiplier,
                                   const Eigen::VectorXd weights,
                                   const bool mixedOnly=false){
  double sse;
  double trace;
  double res; 
  int vec_size;
  Eigen::VectorXd beta;
  Eigen::VectorXd vec;
  SpMat L;
  SpMatC B_inverse;
  SpMat intm;
  SpMat test;
  Eigen::MatrixXd dense1;
  Eigen::MatrixXd dense2;
  
  size_t n_funcs = indexes.size();
  double n_obs = 0;
  double sumw = weights.sum();
  
  // Compute the total amount of observations used
  for (int i = 0; i < n_funcs; i++){
    vec = values[i];
    vec_size = vec.size();
    n_obs += (vec_size+1-mixedOnly) * (vec_size-mixedOnly) / 2;
  }
  
  lambda = 1/multiplier * pow(16, (lambda * 6 - 2));
  SpMat B = lambda * penalty + designX; 
  
  // Compute the cholesky decomposition to solve
  // There maybe some gains using the permuatation matrices,
  // however, it would then later entail another matrix multiplication. 
  // Haven't checked which is faster. The inverse without permutation is
  // only slightly less sparse than the one with permutation.
  Eigen::SimplicialLLT<Eigen::SparseMatrix<double>, Eigen::Lower, Eigen::NaturalOrdering<int>> cholesky;
  cholesky.compute(B);
  
  beta = cholesky.solve(designY);
  
  sse = computeSSEC(indexes, values, lower_int, upper_int, nbasis, order, beta, weights, mixedOnly);
  
  L = SpMat( cholesky.matrixL() );
  
  B_inverse = takahashiDavisC(L);
  
  test = designX * B_inverse;
  trace = test.diagonal().sum();
  
  Rprintf("SSE %f ", sse);
  Rprintf("Trace %f ", trace);
  Rprintf("Lambda %f ", lambda);
  res = (1/sumw) * sse / ((1/n_obs)  * (n_obs - trace) * (1/n_obs)  * (n_obs - trace));
  Rprintf("GCV %f ", res);
  return res;
}



