# A suite of unit tests for the penalized Cross Covariance estimation
library(testthat)
library(fda)
library(Matrix)
library(mvtnorm)
library(sparseinv)
library(Rcpp)

# - Mixed Moments -----------------------------------------------
# Test for mixedMomentsC, simple test for a simple Function

test_that("mixedMoments Computations", {
  expect_equal(mixedMomentsCrC(c(1,2,3), c(1,2,3)), c(1,2,3,2,4,6,3,6,9))
})


# - Spline Basis -----------------------------------------------
# Once there is error handling, this should probably be included
# in the testing 

test_that("Spline Design Matrix Computations", {
  knots <- seq(0, 2000, length.out = 100)
  basis <- create.bspline.basis(knots)
  p = seq(1,1950, length.out = 1500)
  true <- matrix(eval.basis(p, basis), ncol = 100+2)
  test = splineMatrixCrC(0, 2000, 4, 100, p)
  expect_equal(true, test)
})

# - Kronecker Product -----------------------------------------------

test_that("Kronecker Product", {
  knots <- seq(0, 2000, length.out = 100)
  basis <- create.bspline.basis(knots)
  p1 = seq(1,1950, length.out = 1500)
  p2 = sort(runif(1000, min=0,max=2000))
  dM1 = as(matrix(eval.basis(p1, basis), ncol = 100+2), 'dgRMatrix')
  dM2 = as(matrix(eval.basis(p2, basis), ncol = 100+2), 'dgRMatrix')
  true = Matrix::kronecker(dM1, dM2)
  test = singleDesignMatrixCrC(dM1, dM2)
  expect_equal(as(true, 'CsparseMatrix'), as(test, 'CsparseMatrix'))
})

# - Penalty Matrix -----------------------------------------------
# Tests for Penalty Matrix  there should probably be more tests

test_that("Penalty Matrix Computations", {
  knots1 <- seq(0, 2000, length.out = 100)
  knots2 <- seq(50, 100, length.out = 200)
  basis1 <- create.bspline.basis(knots1)
  basis2 <- create.bspline.basis(knots2)
  pen1 = as(bsplinepen(basis1), 'CsparseMatrix')
  pen2 = as(bsplinepen(basis2), 'CsparseMatrix')
  true = Matrix::kronecker(pen1, as(diag(1, ncol(pen2)), 'CsparseMatrix')) + Matrix::kronecker(as(diag(1, ncol(pen1)), 'CsparseMatrix'), pen2)
  test = penaltyMatrixCrC(bsplinepen(basis1), bsplinepen(basis2))
  expect_equal(as(true, 'CsparseMatrix'), as(test, 'CsparseMatrix'))
})

# - Complete Design -----------------------------------------------

test_that("Complete Design Computations", {
  n_funcs = 10
  n1 = 10
  n2 = 12
  indexes1 = lapply(seq(1,n_funcs,1), function(x){return(sort(runif(10, min = 1, max = 2)))})
  indexes2 = lapply(seq(1,n_funcs,1), function(x){return(sort(runif(10, min = 4, max = 5)))})
  values1 = lapply(seq(1,n_funcs,1), function(x){return(sort(runif(10, min = -2, max = 2)))})
  values2 = lapply(seq(1,n_funcs,1), function(x){return(sort(runif(10, min = 4, max = 8)))})
  weights = as.numeric(seq(1,1, length.out = 10))
  knots1 <- seq(1, 2, length.out = n1-2)
  knots2 <- seq(4, 5, length.out = n2-2)
  basis1 <- create.bspline.basis(knots1)
  basis2 <- create.bspline.basis(knots2)
  mats1 = lapply(seq(1,n_funcs,1), function(x){return(as(eval.basis(indexes1[[x]], basis1), 'CsparseMatrix'))})
  mats2 = lapply(seq(1,n_funcs,1), function(x){return(as(eval.basis(indexes2[[x]], basis2), 'CsparseMatrix'))})
  prods = lapply(seq(1,n_funcs,1), function(x){return(Matrix::kronecker(mats1[[x]], mats2[[x]]))})
  mm_true = do.call(c, lapply(seq(1,n_funcs,1), function(x){return(Matrix::kronecker(values1[[x]], values2[[x]]))}))
  Phi_true = do.call(rbind, prods)
  test  = Matrix::t(Phi_true) %*% Phi_true
  result = fullDesignCrC(indexes1, values1, indexes2, values2, as.numeric(c(1,2)),
                       as.numeric(c(4,5)), as.integer(c(10, 12)), as.integer(c(4, 4)),
                       weights)
  expect_equal(as(result[[1]], 'CsparseMatrix'), as(test, 'CsparseMatrix'))
})

test_that("SSE Computations", {
  n_funcs = 100
  nobs1 = 20
  nobs2 = 30
  n1 = 10
  n2 = 12
  indexes1 = lapply(seq(1,n_funcs,1), function(x){return(sort(runif(nobs1, min = 1, max = 2)))})
  indexes2 = lapply(seq(1,n_funcs,1), function(x){return(sort(runif(nobs2, min = 4, max = 5)))})
  values1 = lapply(seq(1,n_funcs,1), function(x){return(sort(runif(nobs1, min = -2, max = 2)))})
  values2 = lapply(seq(1,n_funcs,1), function(x){return(sort(runif(nobs2, min = 4, max = 8)))})
  knots1 <- seq(1, 2, length.out = n1-2)
  knots2 <- seq(4, 5, length.out = n2-2)
  basis1 <- create.bspline.basis(knots1)
  basis2 <- create.bspline.basis(knots2)
  mats1 = lapply(seq(1,n_funcs,1), function(x){return(as(eval.basis(indexes1[[x]], basis1), 'CsparseMatrix'))})
  mats2 = lapply(seq(1,n_funcs,1), function(x){return(as(eval.basis(indexes2[[x]], basis2), 'CsparseMatrix'))})
  prods = lapply(seq(1,n_funcs,1), function(x){return(Matrix::kronecker(mats1[[x]], mats2[[x]]))})
  weights = as.numeric(seq(1,1, length.out = n_funcs))
  b = runif(n1 * n2)
  mm_true = do.call(c, lapply(seq(1,n_funcs,1), function(x){return(Matrix::kronecker(values1[[x]], values2[[x]]))}))
  Phi_true = do.call(rbind, prods)
  st  = proc.time()
  preds_true = Phi_true %*% b
  sse_true = sum((preds_true - mm_true)^2)
  e1 = proc.time()
  sse_test = computeSSECrC(indexes1, values1, indexes2, values2, as.numeric(c(1,2)),
                           as.numeric(c(4,5)), as.integer(c(10,12)), as.integer(c(4,4)),
                           b, weights)
  e2 = proc.time()
  expect_equal(sse_true, sse_test)
})

