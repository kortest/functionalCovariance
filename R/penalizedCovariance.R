#' Compute the coefficients for the covariance kernel in a product bspline basis
#' 
#' This function takes lists of indexes and values of observations of random functions
#' and estimates the marginal covariance kernel in a product bspline basis
#' 
#' @param indexes List of indexes of the quantity observed as a random function
#' @param values List of values of the quantity observed as a random function at indexes
#' @param lower Lower interval boundary for the splines
#' @param upper Upper interval boundary for the splines
#' @param nbasis Number of basis functions 
#' @param splinePenalty Penalty Matrix for one bspline basis
#' @param order The order of the bsplines (default 4)
#' @param weights Vector of weights for the function observations (default equal weights)
#' @param lambda Penalty parameter (default 10)
#' @param use_cv Whether to use generalized cross validation to determine the penalty parameter (default FALSE)
#' @param mixedOnly Whether to exclude estimates for the second moments(i.e. only values[indexes[i]]^2 or only values[indexes[i]] * values[indexes[i]], i != j) for the estimation
#' @param lower_optim Parameter for optimization during gcv, (default = 10^0.25) 
#' @param upper_optim Parameter for optimization during gcv, (upper_optim = 10^6) 
#' @return A matrix of coefficients for the product bspline basis
#' @export
penalizedCovariance <- function(indexes, values, lower, upper, nbasis, order=4, splinePenalty, weights=NULL, lambda = 10, use_cv=FALSE, mixedOnly=FALSE,
                                lower_optim = 10^0.25, upper_optim = 10^6){
  
  # Check the arguments first 
  .check_arguments(indexes, values, lower, upper, nbasis, order, weights=NULL, lambda = 10)
  
  # Use weights = 1 if not given
  if (is.null(weights)){
    weights = as.numeric(seq(1,1, length.out = length(indexes)))
  }
  
  # Compute the design matrices
  fullDesign = fullDesignC(indexes, values, lower, upper, nbasis, order, weights, mixedOnly)
  # This is ugly, but can't be changed until there is c++ implementation
  # of the spline penalties
  # knots <- seq(lower, upper, length.out = nbasis-2)
  # basis <- create.bspline.basis(knots, norder = order)
  # splinePenalty = bsplinepen(basis)
  penalty = penaltyMatrixC(splinePenalty, 2*order^2)
  
  # Check if cv is wanted
  if (use_cv){
    # Perform the optimization
    mult = sum(Matrix::diag(penalty)) / sum(Matrix::diag(fullDesign[[1]]))
    #lower_optim = 10^0.25; upper_optim = 10^6
    tol = 1e-3
    opt_result = stats::optimize(generalizedCrossValidationC,
                          lower = (log(base = 16, mult*lower_optim) + 2)/6,
                          upper = (log(base = 16, mult*upper_optim) + 2)/6,
                          tol = tol, indexes = indexes, values = values, lower_int=lower,
                          upper_int=upper, nbasis = nbasis, order = order, designX = fullDesign[[1]],
                          designY = fullDesign[[2]], penalty=penalty, multiplier = mult,
                          weights = weights, mixedOnly=TRUE)
    lambda = opt_result[[1]]
  }
  
  return(estimateCoefficientsC(fullDesign, penalty, lambda))
}


.check_arguments <- function(indexes, values, lower, upper, nbasis, order, weights=NULL, lambda = 10, use_cv=FALSE){
  if (length(indexes) != length(values)){
    stop("The list of indexes and the list of functions must be of the same length.")
  }
  for(i in 1:length(indexes)){
    if (length(indexes[i]) != length(values[i])){
      stop("Length of indexes and values at indexes should be of same length.")
    }
    # Could probably be sorted inside the function if not already
    if (is.unsorted(indexes[i])){
      stop("Indexes need to be in increasing order.")
    }
    if (!is.numeric(indexes[[i]])){
      stop("Indexes need to be numeric.")
    }
    if (!is.numeric(values[[i]])){
      stop("Values need to be numeric.")
    }
    if (lower > min(indexes[[i]])){
      stop("lower should be smaller than the smallest index")
    }
    if (upper < max(indexes[[i]])){
      stop("upper should be larger than the biggest index")
    }
  }
  if (nbasis - 4 + order <= 0){
    stop("nbasis needs to be larger for this order for splines.")
  }
  if (lower >= upper){
    stop("Lower interval boundary should be strictly smaller than upper boundary.")
  }
}
