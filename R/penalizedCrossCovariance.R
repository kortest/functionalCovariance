#' Compute the coefficients for the cross-covariance kernel in a product bspline basis
#' 
#' This function takes lists of indexes and values of observations of random functions
#' and estimates the marginal cross-covariance kernel in a product bspline basis
#' 
#' @param indexes1 List of indexes of the first quantity observed as a random function
#' @param values1 List of values of the first quantity observed as a random function
#' @param indexes2 List of indexes of the second quantity observed as a random function
#' @param values2 List of values of the second quantity observed as a random function
#' @param intv1 A vector of length two containing the interval boundaries for the splines for first quantity
#' @param intv2 A vector of length two containing the interval boundaries for the splines for second quantity
#' @param nbasis Number of basis functions (currently the same for both quantities)
#' @param order Integer vector of length two determining the order of the bsplines (default c(4,4))
#' @param weights Vector of weights for the function observations (default equal weights)
#' @param lambda Penalty parameter (default 10)
#' @param use_cv Whether to use generalized cross validation to determine the penalty parameter (default FALSE)
#' @param lower_optim Parameter for optimization during gcv, (default = 10^0.25) 
#' @param upper_optim Parameter for optimization during gcv, (upper_optim = 10^6) 
#' @return A matrix of coefficients for the product bspline basis
#' @export
penalizedCrossCovariance <- function(indexes1, values1, indexes2, values2,
                                     intv1, intv2, nbasis, order=as.integer(c(4,4)),
                                     weights=NULL, lambda = 10, use_cv=FALSE,
                                     lower_optim = 10^0.25, upper_optim = 10^6){
  
  # Check the arguments first 
  .check_arguments(indexes1, values1, indexes2, values2, intv1, intv2, nbasis, order, weights=NULL, lambda = 10)
  
  # Use weights = 1 if not given
  if (is.null(weights)){
    weights = as.numeric(seq(1,1, length.out = length(indexes1)))
  }
  
  # Compute the design matrices
  fullDesign = fullDesignCrC(indexes1, values2, indexes1, values2, intv1, intv2, nbasis, order, weights)
  # This is ugly, but can't be changed until there is c++ implementation
  # of the spline penalties
  knots1 <- seq(intv1[1], intv1[2], length.out = nbasis[1]-2)
  basis1 <- fda::create.bspline.basis(knots1, norder = order[1])
  splinePenalty1 = fda::bsplinepen(basis1)
  knots2 <- seq(intv2[1], intv2[2], length.out = nbasis[2]-2)
  basis2 <- fda::create.bspline.basis(knots2, norder = order[2])
  splinePenalty2 = fda::bsplinepen(basis2)
  penalty = penaltyMatrixCrC(splinePenalty1, splinePenalty2)
  
  # Check if cv is wanted
  if (use_cv){
    # Perform the optimization
    mult = sum(diag(penalty)) / sum(diag(fullDesign[[1]]))
 #   lower_optim = 10^0.25; upper_optim = 10^6
    tol = 1e-3
    opt_result = stats::optimize(generalizedCrossValidationCrC,
                          lower = (log(base = 16, mult*lower_optim) + 2)/6,
                          upper = (log(base = 16, mult*upper_optim) + 2)/6,
                          tol = tol, indexes1 = indexes1, values1 = values1, indexes2 = indexes2,
                          values2 = values2, intv1=intv1, intv2=intv2, nbasis = nbasis, order = order,
                          designX = fullDesign[[1]], designY = fullDesign[[2]], penalty=penalty,
                          multiplier = mult, weights = weights)
    lambda = opt_result[[1]]
  }
  
  return(estimateCoefficientsCrC(fullDesign, penalty, lambda))
}


.check_arguments <- function(indexes1, values1, indexes2, values2, intv1,
                             intv2, nbasis, order, weights=NULL, lambda = 10, use_cv=FALSE){
  if (length(indexes1) != length(values1)){
    stop("The list of indexes1 and the list of functions1 must be of the same length.")
  }
  if (length(indexes2) != length(values2)){
    stop("The list of indexes2 and the list of functions2 must be of the same length.")
  }
  if (length(indexes1) != length(indexes1)){
    stop("The list of indexes1 and the list of indexes2 must be of the same length.")
  }
  if(length(intv1) != 2){
    stop("Intv1 needs to be be a numeric vector of length 2")
  }
  if(length(nbasis) != 2){
    stop("The number of basis functions needs to be provided for each of the observation lists")
  }
  if(length(order) != 2){
    stop("The order needs to be provided for each of the observation lists")
  }
  if(intv1[1] > intv1[2]){
    stop("The first entry of intv1 needs to be smaller or equal to the second entry")
  }
  if(intv2[1] > intv2[2]){
    stop("The first entry of intv2 needs to be smaller or equal to the second entry")
  }
  for(i in 1:length(indexes1)){
    if (length(indexes1[i]) != length(values1[i])){
      stop("Length of indexes and values at indexes should be of same length.")
    }
    if (length(indexes2[i]) != length(values2[i])){
      stop("Length of indexes and values at indexes should be of same length.")
    }
    # Could probably be sorted inside the function if not already
    if (is.unsorted(indexes1[i])){
      stop("Indexes need to be in increasing order.")
    }
    if (!is.numeric(indexes1[[i]])){
      stop("Indexes need to be numeric.")
    }
    if (is.unsorted(indexes2[i])){
      stop("Indexes need to be in increasing order.")
    }
    if (!is.numeric(indexes2[[i]])){
      stop("Indexes need to be numeric.")
    }
    if (!is.numeric(values1[[i]])){
      stop("Values need to be numeric.")
    }
    if (!is.numeric(values2[[i]])){
      stop("Values need to be numeric.")
    }
    if (intv1[1] > min(indexes1[[i]])){
      stop("lower should be smaller than the smallest index")
    }
    if (intv2[1] > min(indexes2[[i]])){
      stop("lower should be smaller than the smallest index")
    }
    if (intv1[2] < max(indexes1[[i]])){
      stop("upper should be larger than the biggest index")
    }
    if (intv2[2] < max(indexes2[[i]])){
      stop("upper should be larger than the biggest index")
    }
  }
  if (nbasis[1] - 4 + order[1] <= 0){
    stop("nbasis needs to be larger for this order for splines.")
  }
  if (nbasis[2] - 4 + order[2] <= 0){
    stop("nbasis needs to be larger for this order for splines.")
  }
}
