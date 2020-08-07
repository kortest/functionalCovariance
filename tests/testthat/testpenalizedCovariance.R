# A suite of unit tests for the penalized Covariance estimation
library(testthat)
library(fda)
library(Matrix)
library(mvtnorm)
library(sparseinv)
library(Rcpp)


# - Mixed Moments -----------------------------------------------
# Test for mixedMomentsC, simple test for a simple Function

test_that("mixedMoments Computations", {
  expect_equal(mixedMomentsC(c(1,2,3), TRUE), c(2,3,6))
  expect_equal(mixedMomentsC(c(1,2,3), FALSE), c(1,2,3,4,6,9))
  expect_equal(mixedMomentsC(c(1,2,3)), c(1,2,3,4,6,9))
})


# - Spline Basis -----------------------------------------------
# Once there is error handling, this should probably be included
# in the testing 

test_that("Spline Design Matrix Computations", {
  knots <- seq(0, 2000, length.out = 100)
  basis <- create.bspline.basis(knots)
  p = seq(1,1950, length.out = 1500)
  true <- matrix(eval.basis(p, basis), ncol = 100+2)
  test = splineMatrixC(0, 2000, 4, 100, p)
  expect_equal(true, test)
})


# - Design Matrix -----------------------------------------------
# Tests for singleDesignMatrixC, there should probably be more tests

test_that("Single Design Matrix Computations", {
  n = 50
  p_length = 300
  # Create a symmetric coefficients matrix
  b = matrix(runif(n^2, min = -1, max = 1), ncol = n)
  b = t(b) %*% b
  # Only take the lower triangle
  b_t = as.vector(tril(b)[tril(b) != 0])
  # As coefficient vectors
  b = as.vector(t(b))
  # Create the spline matrix
  knots <- seq(0, 2000, length.out = n-2)
  basis <- create.bspline.basis(knots)
  p = seq(1,1950, length.out = p_length)
  splineDesign = as(eval.basis(p, basis), 'dgRMatrix')
  # Create the design matrix via Kronecker product 
  # and remove the unneeded lines
  kr = Matrix::kronecker(splineDesign, splineDesign)
  indexes = integer()
  counter = 1
  indcounter = 1
  for (i in 1:p_length){
    for (j in 1:p_length){
      if (j > i){
        indexes[counter] = indcounter
        counter = counter + 1
      }
      indcounter = indcounter + 1
    }
  }
  kr = kr[indexes,]  # 11.5 MB
  # Create the design matrix using the symmetry of the coefficients
  test = singleDesignMatrixC(splineDesign, 16, TRUE)   # 8.7 MB
  expect_equal(kr %*% b, test %*% b_t)
})


# - Penalty Matrix -----------------------------------------------
# Tests for Penalty Matrix  there should probably be more tests

test_that("Penalty Matrix Computations", {
  n = 50
  p_length = 300
  # Create a symmetric coefficients matrix
  b = matrix(runif(n^2, min = -1, max = 1), ncol = n)
  b = t(b) %*% b
  # Only take the lower triangle
  b_t = as.vector(tril(b)[tril(b) != 0])
  # As coefficient vectors
  b = as.vector(t(b))
  # Create the spline penalty matrix
  knots <- seq(0, 2000, length.out = n-2)
  basis <- create.bspline.basis(knots)
  p = seq(1,1950, length.out = p_length)
  splinePenalty = bsplinepen(basis)
  # Create the design matrix via Kronecker product 
  penalty_true = as(Matrix::kronecker(splinePenalty, diag(1, ncol = n, nrow = n)) + Matrix::kronecker(splinePenalty, diag(1, ncol = n, nrow = n)), 'CsparseMatrix') 
  # Create the smaller design matrix
  test = penaltyMatrixC(splinePenalty, 20)
  # The resulting penalty matrix has a factor 2 in every element. 
  # To speed up, I ignored this
  expect_equal(t(b) %*% penalty_true %*% b, 2*(t(b_t) %*% test %*% b_t))
})


# - Complete Matrix -----------------------------------------------
# Test for t(Phi) %*% Phi, where Phi is the complete design matrix for 
# all functional observations. The tests become increasingly complex
# to set up.

test_that("Complete Design Computations", {
  n = 50
  p_length = 200
  n_funcs = 10
  order = 4
  # Should at some point be amended to more meaningful weights
  weights = as.numeric(seq(1,1, length.out = n_funcs))
  # Create randomized observation indexes
  func_indexes = lapply(seq(1,n_funcs,1), function(x){ return(sort(runif(p_length, min = 0, max = 1)))}) 
  # Simulate from a Brownian motion
  funcs = lapply(func_indexes, function(x){
    Sigma = matrix(nrow = p_length, ncol = p_length)
    for (i in 1:p_length){
      for (j in i:p_length){
        Sigma[i,j] = x[i]
        Sigma[j,i] = x[i]
      }
    }
    return(rmvnorm(1, sigma = Sigma))
  })
  # Compute the mixed moment estimates. Only strictly increasing as with ARGO 
  # but could of course also be done differently.
  mixed_moments = lapply(funcs, function(x){
    y = numeric()
    counter=1
    for (i in 1:p_length){
      for (j in 1:p_length){
        if (i != j){
          y[counter] = x[i] * x[j]
          counter = counter + 1  
        }
      }
    }
    return(y)
  })
  # Create the spline matrix
  knots <- seq(0, 1, length.out = n-2)
  basis <- create.bspline.basis(knots)
  p = seq(0.00001,0.9999, length.out = p_length)
  splineDesigns = lapply(func_indexes, function(x){return(as(eval.basis(x, basis), 'CsparseMatrix'))})
  # Get the relevant indexes
  indexes = integer()
  counter = 1
  indcounter = 1
  for (i in 1:p_length){
    for (j in 1:p_length){
      if (j != i){
        indexes[counter] = indcounter
        counter = counter + 1
      }
      indcounter = indcounter + 1
    }
  }
  Phi_true = do.call(rbind,lapply(splineDesigns, function(x){
    kr = Matrix::kronecker(x, x)
    return(kr[indexes,])
  }))
  Y_true = as.matrix(do.call(c, mixed_moments))
  Y_true = Matrix::t(Phi_true) %*% Y_true
  Phi_true = Matrix::crossprod(Phi_true) # t(Phi_true) %*% Phi_true
  beta_true = matrix(Matrix::solve(Phi_true, Y_true), ncol=50)
  # The estimate without using the symmetry should still be a symmetric matrix
  # Oddly, R sometimes doesn't say it's not symmetric but doesn't complain later 
  # when compared to a forced symmetric matrix. Commented out for now
  # expect_true(isSymmetric(beta_true))
  # Create the design matrix using the symmetry of the coefficients
  result = fullDesignC(func_indexes, funcs, 0.0, 1.0, n, 4, weights, TRUE) 
  # Solve the new system of equations and transform to matrix
  beta_test = matrix(0, n, n)
  beta_test[lower.tri(beta_test, diag=TRUE)] = solve(as(result[[1]], 'CsparseMatrix'), result[[2]])
  beta_test = as.matrix(forceSymmetric(beta_test, 'L'))
  # The resulting estimates should be equal.
  # as.numeric because R has so much flexibility that 
  # it is hard to compare
  expect_equal(as.numeric(beta_test), as.numeric(beta_true))
})


# - Check output of the solver -----------------------------------------------
# Comparing against my own code here because it is not obvious to me 
# what to compare it against. But the used code is checked above.

test_that("Solver comparison", {
  n = 50
  p_length = 200
  n_funcs = 10
  order = 4
  weights = as.numeric(seq(1,1,length.out = n_funcs))
  # Create randomized observation indexes
  func_indexes = lapply(seq(1,n_funcs,1), function(x){ return(sort(runif(p_length, min = 0, max = 1)))}) 
  # Simulate from a Brownian motion
  funcs = lapply(func_indexes, function(x){
    Sigma = matrix(nrow = p_length, ncol = p_length)
    for (i in 1:p_length){
      for (j in i:p_length){
        Sigma[i,j] = x[i]
        Sigma[j,i] = x[i]
      }
    }
    return(rmvnorm(1, sigma = Sigma))
  })
  # Create the design matrix using the symmetry of the coefficients
  fullDesign = fullDesignC(func_indexes, funcs, 0.0, 1.0, n, 4, weights, TRUE) 
  # Solve the system of equations and use R's solver
  knots <- seq(0, 1, length.out = n-2)
  basis <- create.bspline.basis(knots)
  penalty = bsplinepen(basis)
  jointPenalty = penaltyMatrixC(penalty, 2*order*order)
  beta_true = solve(as(fullDesign[[1]], 'CsparseMatrix') + jointPenalty, fullDesign[[2]])
  # Check output of Eigens solver
  beta_test = estimateCoefficientsC(fullDesign, jointPenalty, 1.0)
  expect_equal(as.numeric(beta_test), as.numeric(beta_true))
})

# - Check the inverse computation -----------------------------------------------
# Comparing my output against that of the R package.

test_that("Takahashi Davis Inverse", {
  n = 50
  p_length = 200
  n_funcs = 10
  order = 4
  weights = as.numeric(seq(1,1, length.out = n_funcs))
  # Create randomized observation indexes
  func_indexes = lapply(seq(1,n_funcs,1), function(x){ return(sort(runif(p_length, min = 0, max = 1)))}) 
  # Simulate from a Brownian motion
  funcs = lapply(func_indexes, function(x){
    Sigma = matrix(nrow = p_length, ncol = p_length)
    for (i in 1:p_length){
      for (j in i:p_length){
        Sigma[i,j] = x[i]
        Sigma[j,i] = x[i]
      }
    }
    return(rmvnorm(1, sigma = Sigma))
  })
  # Create the design matrix
  result = fullDesignC(func_indexes, funcs, 0.0, 1.0, n, 4, weights, TRUE) 
  # Create the penalty
  knots <- seq(0, 1, length.out = n-2)
  basis <- create.bspline.basis(knots)
  pen = penaltyMatrixC(bsplinepen(basis), 20)
  B = result[[1]] + pen
  # Apparently there is transformation from dgC to dgR directly
  B = as(as.matrix(B), 'dgRMatrix')
  chol = Cholesky(B, perm=FALSE)
  L = expand(chol)$L
  Lr = as(as.matrix(L), 'dgRMatrix')
  # Because the takahashi forumla only really needs L,
  # I wrote it this way, I commented on this in the c++ file
  inv_test = takahashiDavisC(Lr)
  # I force no permuatation here because I think this faster 
  # overall in our use case and thus implemented it that way
  # on the C side of things.
  inv_true = Takahashi_Davis(B, L, P = diag(1, nrow=nrow(B)))
  expect_equal(as.matrix(inv_test), as.matrix(inv_true))
})

# - Compute Sum Squared Error -----------------------------------------------
# Check the SSE computations

test_that("SSE Computations", {
  n = 50
  p_length = 200
  n_funcs = 10
  order = 4
  weights = as.numeric(seq(1,1, length.out = n_funcs))
  # Get random coefficients
  b = matrix(runif(n^2, min = -1, max = 1), ncol = n)
  b = t(b) %*% b
  # Only take the lower triangle
  b_t = as.vector(tril(b)[tril(b) != 0])
  # As coefficient vectors
  b = as.vector(t(b))
  # Create randomized observation indexes
  func_indexes = lapply(seq(1,n_funcs,1), function(x){ return(sort(runif(p_length, min = 0, max = 1)))}) 
  # Simulate from a Brownian motion
  funcs = lapply(func_indexes, function(x){
    Sigma = matrix(nrow = p_length, ncol = p_length)
    for (i in 1:p_length){
      for (j in i:p_length){
        Sigma[i,j] = x[i]
        Sigma[j,i] = x[i]
      }
    }
    return(rmvnorm(1, sigma = Sigma))
  })
  # Compute the mixed moment estimates. Only strictly increasing as with ARGO 
  # but could of course also be done differently.
  mixed_moments = lapply(funcs, function(x){
    y = numeric()
    counter=1
    for (i in 1:p_length){
      for (j in 1:p_length){
        if (i != j){
          y[counter] = x[i] * x[j]
          counter = counter + 1  
        }
      }
    }
    return(y)
  })
  # Create the spline matrix
  knots <- seq(0, 1, length.out = n-2)
  basis <- create.bspline.basis(knots)
  p = seq(0.00001,0.9999, length.out = p_length)
  splineDesigns = lapply(func_indexes, function(x){return(as(eval.basis(x, basis), 'CsparseMatrix'))})
  # Get the relevant indexes
  indexes = integer()
  counter = 1
  indcounter = 1
  for (i in 1:p_length){
    for (j in 1:p_length){
      if (j != i){
        indexes[counter] = indcounter
        counter = counter + 1
      }
      indcounter = indcounter + 1
    }
  }
  Phi_true = do.call(rbind,lapply(splineDesigns, function(x){
    kr = Matrix::kronecker(x, x)
    return(kr[indexes,])
  }))
  preds_true = Phi_true %*% b 
  Y_true = do.call(c, mixed_moments)
  sse_true = sum((preds_true-Y_true)^2)
  # Multiply by two because the R version used every obs twice.
  sse_test = 2*computeSSEC(func_indexes, funcs, 0.0, 1.0, n, 4, b_t, weights, TRUE)
  expect_equal(sse_true, sse_test)
})

