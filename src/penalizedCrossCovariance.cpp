#include <Rcpp.h>
#include <RcppEigen.h>
#include <gsl/gsl_bspline.h>
#include <math.h>

using namespace Rcpp;

typedef Eigen::MappedSparseMatrix<double> MSpMat;
typedef Eigen::SparseMatrix<double, Eigen::RowMajor> SpMat;
typedef Eigen::SparseMatrix<double, Eigen::ColMajor> SpMatC;
typedef Eigen::Triplet<double, int> T;

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
NumericVector mixedMomentsCrC(NumericVector x, NumericVector y) {
  int n1 = x.size();
  int n2 = y.size();
  // Booleans in Rcpp are translated to integers
  NumericVector out(n1 * n2);
  int counter = 0;
  for(int i = 0; i < n1; i++){
    for(int j = 0; j < n2; j++){
      out[counter] = x[i]*y[j];
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
NumericMatrix splineMatrixCrC(const double lower, const double upper, const int order, const int nbreaks, const NumericVector x){
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

// [[Rcpp::export]]
const SpMat singleDesignMatrixCrC(const SpMat dM1, const SpMat dM2) {
  return Eigen::kroneckerProduct(dM1, dM2).eval();
}

// [[Rcpp::export]]
SpMat penaltyMatrixCrC(const NumericMatrix X, const NumericMatrix Y){
  int n1 = X.cols();
  int n2 = Y.cols();
  
  SpMat res1 = SpMat(n1, n1);
  SpMat res2 = SpMat(n2, n2);
  res1.setIdentity();
  res2.setIdentity();
  
  SpMat int1 = Rcpp::as<Eigen::MatrixXd>(X).sparseView();
  SpMat int2 = Rcpp::as<Eigen::MatrixXd>(Y).sparseView();;
  
  res2 = Eigen::kroneckerProduct(int1, res2).eval();
  res1 = Eigen::kroneckerProduct(res1, int2).eval();
  
  return res1 + res2;
}

// [[Rcpp::export]]
List fullDesignCrC(const List indexes1,
                 const List values1,
                 const List indexes2,
                 const List values2,
                 const NumericVector intv1,
                 const NumericVector intv2,
                 const IntegerVector nbasis,
                 const IntegerVector order,
                 const Eigen::VectorXd weights){
  // Initialize...
  // I think I haven't fully understood Eigen but sometimes
  // it seems to require intermediate assignments
  NumericMatrix spline_mat1;
  NumericMatrix spline_mat2;
  Eigen::MatrixXd intm1;
  Eigen::MatrixXd intm2;
  SpMat Phi;
  SpMat res2;
  SpMat spline_sparse1;
  SpMat spline_sparse2;
  Eigen::VectorXd res2Y;
  
  size_t n_funcs = indexes1.size();
  
  // Initialize the return values and set them 0
  int res_size = nbasis[0] * nbasis[1];
  SpMat res = SpMat(res_size, res_size);
  Eigen::VectorXd resY = Eigen::VectorXd(res_size);
  res.setZero();
  resY.setZero();
  
  // Loop through all the profiles...
  // This somehow begs for paralleled implementation
  // but this really depends on the whole structure
  for (size_t i = 0; i < n_funcs; i++){
    spline_mat1 = splineMatrixCrC(intv1[0], intv1[1], order[0], nbasis[0] - 2, indexes1[i]);
    spline_mat2 = splineMatrixCrC(intv2[0], intv2[1], order[0], nbasis[1] - 2, indexes2[i]);
    intm1 = Rcpp::as<Eigen::MatrixXd>(spline_mat1);
    intm2 = Rcpp::as<Eigen::MatrixXd>(spline_mat2);
    spline_sparse1  = intm1.sparseView();
    spline_sparse2  = intm2.sparseView();
    // The order squared will always be the least upper bound for Bsplines
    Phi = singleDesignMatrixCrC(spline_sparse1, spline_sparse2);
    // Calculate the matrix to solve for the coefficients
    // Time-wise by FAR the most expensive line
    // Eigen gives 7 seconds over Armadillo but there might sill be potential.
    // This would be a next step...
    res2 = res.selfadjointView<Eigen::Lower>().rankUpdate(Phi.adjoint(), weights[i]);
    res = res2;
    // Compute the transformation on the response vectors
    res2Y = Phi.adjoint() * Rcpp::as<Eigen::VectorXd>(mixedMomentsCrC(values1[i], values2[i]));
    resY = resY + weights[i] * res2Y;
  }
  return Rcpp::List::create(res, resY);
}

// [[Rcpp::export]]
Eigen::VectorXd estimateCoefficientsCrC(const List fullDesign,
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

NumericVector sparseinv2Cr(int n, const int *Lp, const int *Li, const double *Lx, Eigen::VectorXd d, const int *Up, const int *Uj, const double *Ux, const int *Zp, const int *Zi) {
  
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
SpMatC takahashiDavisCrC(const SpMat L){
  
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
  Eigen::VectorXd Zx = Rcpp::as<Eigen::VectorXd>(sparseinv2Cr(nrows, LnewC.outerIndexPtr(), LnewC.innerIndexPtr(), LnewC.valuePtr(), d, LnewC.outerIndexPtr(), LnewC.innerIndexPtr(), LnewC.valuePtr(), Z.outerIndexPtr(), Z.innerIndexPtr()));
  
  // Now construct the sparseMatrix to be returned, unfortunately col major....
  Eigen::Map<Eigen::SparseMatrix<double, Eigen::ColMajor>> spMap(nrows, nrows, Z.nonZeros(), Z.outerIndexPtr(), Z.innerIndexPtr(), Zx.data(), 0);
  Eigen::SparseMatrix<double, Eigen::ColMajor> res= spMap.eval();
  return res;
}



// [[Rcpp::export]]
double computeSSECrC(const List indexes1,
                   const List values1,
                   const List indexes2,
                   const List values2,
                   const NumericVector intv1,
                   const NumericVector intv2,
                   const IntegerVector nbasis,
                   const IntegerVector order,
                   const Eigen::VectorXd beta,
                   const Eigen::VectorXd weights){
  Eigen::VectorXd preds;
  Eigen::VectorXd mm;
  Eigen::MatrixXd intm1;
  Eigen::MatrixXd intm2;
  NumericMatrix spline_mat1;
  NumericMatrix spline_mat2;
  SpMat spline_sparse1;
  SpMat spline_sparse2;
  SpMat Phi;
  size_t n_funcs = indexes1.size();
  double res = 0; 
  
  for (size_t i = 0; i < n_funcs; i++){
    spline_mat1 = splineMatrixCrC(intv1[0], intv1[1], order[0], nbasis[0] - 2, indexes1[i]);
    spline_mat2 = splineMatrixCrC(intv2[0], intv2[1], order[1], nbasis[1] - 2, indexes2[i]);
    intm1 = Rcpp::as<Eigen::MatrixXd>(spline_mat1);
    intm2 = Rcpp::as<Eigen::MatrixXd>(spline_mat2);
    spline_sparse1  = intm1.sparseView();
    spline_sparse2  = intm2.sparseView();
    Phi = singleDesignMatrixCrC(spline_sparse1, spline_sparse2);
    mm = Rcpp::as<Eigen::VectorXd>(mixedMomentsCrC(values1[i], values2[i]));
    preds = Phi * beta;
    res += weights[i] * (mm - preds).array().square().sum();
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
double generalizedCrossValidationCrC(double lambda,
                                   const List indexes1,
                                   const List values1,
                                   const List indexes2,
                                   const List values2,
                                   const NumericVector intv1,
                                   const NumericVector intv2,
                                   const IntegerVector nbasis,
                                   const IntegerVector order,
                                   const SpMat designX,
                                   const Eigen::VectorXd designY,
                                   const SpMat penalty,
                                   const double multiplier,
                                   const Eigen::VectorXd weights){
  double sse;
  double trace;
  double res; 
  int vec_size1;
  int vec_size2;
  Eigen::VectorXd beta;
  Eigen::VectorXd vec1;
  Eigen::VectorXd vec2;
  SpMat L;
  SpMatC B_inverse;
  SpMat intm;
  SpMat test;
  Eigen::MatrixXd dense1;
  Eigen::MatrixXd dense2;
  
  size_t n_funcs = indexes1.size();
  double n_obs = 0;
  
  // Compute the total amount of observations used
  for (int i = 0; i < n_funcs; i++){
    vec1 = values1[i];
    vec_size1 = vec1.size();
    vec2 = values2[i];
    vec_size2 = vec2.size();
    n_obs += vec_size1 * vec_size2;
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
  
  sse = computeSSECrC(indexes1, values1, indexes2, values2, intv1, intv2, nbasis, order, beta, weights);
  
  L = SpMat( cholesky.matrixL() );
  
  B_inverse = takahashiDavisCrC(L);
  
  test = designX * B_inverse;
  trace = test.diagonal().sum();
  
  Rprintf("Result1 %f ", sse);
  Rprintf("Result2 %f ", trace);
  Rprintf("Result3 %f ", lambda);
  res = sse / ((1/n_obs)  * (n_obs - trace) * (1/n_obs)  * (n_obs - trace));
  Rprintf("Result4 %f ", res);
  return res;
}

