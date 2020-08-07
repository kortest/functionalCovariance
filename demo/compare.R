# implement grid
# library(argofda) # requires package from paper
library(Matrix)
library(fda)
library(sparseinv)
library(dplyr)
library(Rcpp)
library(RcppGSL)
library(mvtnorm)

# Set seed, change if required. 
set.seed(100)

n = 102
p_length = 50
n_funcs = 500
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

# choice of knots
knots <- seq(0, 1, length.out = n-2)
basis <- create.bspline.basis(knots)
penalty_mat <- bsplinepen(basis)

profiles = rep(seq(1,n_funcs,1),seq(p_length,p_length, length.out = n_funcs))
weights = seq(1,1, length.out = n_funcs)

df <- data.frame(pressure=do.call(c, func_indexes),
                 temperature=do.call(c, funcs), 
                 profile=profiles, 
                 wt=weights) 

profile_lengths <- table(df$profile)

# number of naive cross products
sum(profile_lengths * (profile_lengths-1))

# number of cross products after reducing to 50 from each profile
profile_lengths_50 <- ifelse(profile_lengths > 50, 50, profile_lengths)
sum(profile_lengths_50 * (profile_lengths_50-1))

st = proc.time()
test1 <- cov_est_example(df = df, mi_max = 50, knots = knots, lambda = 3,
                         basis = basis, penalty_mat = penalty_mat)
e1 = proc.time()
test2 = penalizedCovariance(func_indexes, funcs, 0.0, 1.0, n, order=4, weights=NULL, lambda = 10, use_cv=TRUE, mixedOnly=FALSE)
e2 = proc.time()

# Look at matrix of coefficients
beta = test1[[1]][[1]]
beta = do.call(c, beta)
beta_mat = matrix(beta, ncol=102, nrow = 102)
image(beta_mat)

beta2_mat = matrix(0, n, n)
beta2_mat[lower.tri(beta2_mat, diag=TRUE)] = test2
beta2_mat = as.matrix(forceSymmetric(beta2_mat, 'L'))
image(beta2_mat)

# Create matrix of actual covariances
test_length = 100

# Design matrix
p = seq(0, 1, length.out = test_length)
spline = as(eval.basis(p, basis), 'dgRMatrix')
kronecker = Matrix::kronecker(spline, spline)

# Old method
res1 = kronecker %*% beta
res1mat = matrix(res1, test_length, test_length)
image(res1mat)

# New method
beta2_vec = as.vector(beta2_mat)
res2 = kronecker %*% beta2_vec
res2mat = matrix(res2, test_length, test_length)
image(res2mat)

# Compare the actual numbers:
true = matrix(0, test_length, test_length)
for (i in 1:test_length){
  for (j in 1:test_length){
    true[i,j] = min(p[i], p[j])
  }
}

frob1 = norm(true - res1mat)
frob2 = norm(true - res2mat)
print(frob1)
print(frob2)

compare_results <- data.frame(
  expand.grid(p1 = p, p2 = p),
  res1 = as.vector(res1mat),
  res2 = as.vector(res2mat),
  true = as.vector(true)
)

library(ggplot2)

a <- ggplot(data = compare_results, aes(x = p1, y = p2, fill = res1-true))+
  geom_raster()+
  scale_fill_gradient2(limits = c(-.1, .1))
b <- ggplot(data = compare_results, aes(x = p1, y = p2, fill = res2-true))+
  geom_raster()+
  scale_fill_gradient2(limits = c(-.1, .1))

library(patchwork)
a+b
ggplot(data = compare_results, aes(x = p1, y = p2, 
                                   fill = abs(res1-true)-abs(res2-true)))+
  geom_raster()+
  scale_fill_gradient2(limits = c(-.08, .08))+
  labs(fill = 'Comparative Error',
       subtitle = 'Positive values means using all data does better')

