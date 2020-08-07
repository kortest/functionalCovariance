#' @export
cov_est_example <-  function(df, mi_max = 50, knots, lambda, basis, penalty_mat) {
  # reduce size of data by randomly selecting 50 observations
  a <- proc.time()
  df_split <- split(df, df$profile)
  df_lens <- sapply(df_split, nrow)
  new_rows <- list()
  # for each profile
  for (q in 1:length(df_split)) {
    if (df_lens[q] < mi_max) {
      new_rows[[q]] <- df_split[[q]]
      next
    }
    new_rows[[q]] <- df_split[[q]][sample(1:(df_lens[q]), size = mi_max),]
  }
  rm(df_split)
  df <- dplyr::bind_rows(new_rows)
  
  n <- nrow(df) 
  nbasis <- length(knots) + 2
  design_mat_prelim <- as(eval.basis(df$pressure, 
                                     basis), 'sparseMatrix')
  # weight matrix
  W <- lapply(1:max(df$profile), function(x) {
    return(Diagonal(n = sum(df$profile == x)^2,
                    sqrt(kronecker(df$wt[df$profile == x],
                                   df$wt[df$profile == x]))))})
  W <- do.call(bdiag, W)
  
  # response vector
  Y <-  lapply(1:max(df$profile), function(x) {
    return(kronecker(df$temperature[df$profile == x],
                     df$temperature[df$profile == x]))}) 
  Y <- do.call(c, Y)
  
  # Create design matrix
  Phi_prelim <- lapply(1:max(df$profile), function(x) {
    return(Matrix::kronecker(design_mat_prelim[df$profile == x,],
                             design_mat_prelim[df$profile == x,]))})  
  prof_lengths <- aggregate(df$profile, list(df$profile), length)
  Phi_init <- lapply(Phi_prelim, Matrix::summary) 
  prod_lengths <- sapply(Phi_prelim, nrow)
  rm(Phi_prelim)
  prod_lengths_sum <- cumsum(prod_lengths)
  Phi_init <- lapply(1:length(Phi_init),function(x) {
    if (x ==1) {
      return(Phi_init[[1]])
    } else {
      Phi_init[[x]][,1] <- Phi_init[[x]][,1] + prod_lengths_sum[x-1]
    }
    return(Phi_init[[x]])
  })
  Phi_init <- bind_rows(Phi_init) # 177 mb
  Phi_init <- sparseMatrix(i = Phi_init[,1], j= Phi_init[,2], #-88.6, #2.04 overall
                           x = Phi_init[,3], dims = c(sum(prod_lengths), nbasis^2))
  
  # Penalty Matrix
  Omega <- Matrix::kronecker(penalty_mat, Diagonal(n = nrow(penalty_mat), x = 1))+
    Matrix::kronecker(Diagonal(n = nrow(penalty_mat), x = 1),penalty_mat)
  
  # Remove where j = k
  check_variance <- lapply(1:max(df$profile), function(x) {
    return(expand.grid('p1' = c(paste0(x, '.', 1:sum(df$profile == x))),
                       'p2' = c(paste0(x, '.', 1:sum(df$profile == x)))))}) 
  check_variance <- bind_rows(check_variance)
  rm('df') # - 1Mb
  
  Y <- Y[check_variance$p1 != check_variance$p2]
  Phi_init <- Phi_init[check_variance$p1 != check_variance$p2,]
  W<- W[check_variance$p1 != check_variance$p2,check_variance$p1 != check_variance$p2]
  
  # Compute preliminary matrices for cross validation
  mod <- model_set_up(covariates = matrix(nrow = length(Y), ncol = 1, 1), 
                      lambdas = 10^c(rep(lambda, nbasis))/35^2,
                      Phi = Phi_init,Psi = Omega, W = W,Y = Y) #268 MB Phi_all (same as Phi_init)
  rm('Phi_init')
  rm('Omega')
  # cross validate
  spline_results <- complete_optim(U = mod[['U']], V = mod[['V']], W = W, 
                                   Omega = mod[['Omega']], Y = Y, 
                                   Phi_all = mod[['Phi_all']],
                                   lower_optim = 10^0.25, upper_optim = 10^6, tol = 10^-1, 
                                   selection_vector = 1:length(mod[['reorder']]),
                                   n_basis = ncol(design_mat_prelim), knots = knots) #nothing
  b <- proc.time()
  print(c(n, (b-a)[3],length(Y)))
  return(list(spline_results,
              c(n, (b-a)[3], length(Y))))
}
