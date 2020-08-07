#include <Rcpp.h>
#include <RcppEigen.h>

typedef Eigen::MappedSparseMatrix<double> MSpMat;
typedef Eigen::SparseMatrix<double, Eigen::RowMajor> SpMat;
typedef Eigen::SparseMatrix<double, Eigen::ColMajor> SpMatC;
typedef Eigen::Triplet<double, int> T;
