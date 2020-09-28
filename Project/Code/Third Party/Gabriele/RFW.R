library(ramify) # For randn
library(tictoc) # For time tracking
library(tibble)
rm(list = ls())

#### Dataset creation according to Kedreux et al. paper
generate_datset = function(n, d, overlap_size, group_size)
{
  set.seed(42)
  
  beta = 14 # Regularizing parameter of the paper.
  gt_nonzero = n * 0.01
  
  # Caligraphic G which in the literature of LGL
  # is the set of all groups
  cal_G = c()
  
  # The number of groups.
  number_of_groups = as.integer((d - group_size) / (group_size - overlap_size))
  
  # The partitioning consists of creating groups
  # where each element is labeled with its own index.
  # We also manage the overlapping of the groups.
  
  # Group index
  i = 1
  cal_G[[i]] = c(1:group_size)
  
  # Current covariate index, which after the creation of
  #  the first group is equal to the group size itself.
  k = group_size
  
  for(i in 2:number_of_groups)
  {
    k = (k - overlap_size + 1)
    k_end = (k + group_size-1)
    cal_G[[i]] = c(k : k_end)
    
    i = i + 1
    k = k_end
  }
  
  "From the paper: We chose the ground truth parameter vector w0 ∈ conv(A) with a fraction of 0.01 of 
  nonzero coefficients, where on eachactive group, the coefficients are generated 
  from a Gaussian distribution
  "
  
  w0 = zeros(d)
  g0 = c() 
  # For the indices related to the groups with non zero coef
  # of the ground truth parameter vector
  
  # Selecting groups with nonzero coef
  w0_indices = sample(1 : length(cal_G), gt_nonzero) # 0.01 of 1000
  
  # w0_indices are what the paper call the "active groups", whose values
  # come frome a gaussian distribution
  
  # We now update the nonzero coefficients related to
  # the groups extracted
  for(i in seq(1:length(w0_indices)))
  {
    # Generate gaussian values
    gt_weight_values = randn(group_size, mean = 0, sd = 1)
    
    # We get the weights of each group randomly chosen and 
    # we update them 
    chosen_group = cal_G[[w0_indices[i]]]
    w0[chosen_group] = w0[chosen_group] + gt_weight_values
    g0[[i]] = gt_weight_values
  }
  
  # Matrix A of GL in the Obiozinski literature is referred as X,
  # but the core remains the same
  
  A = randn(n, d, mean = 0, sd = 1)
  
  # We suppose that the "some additive noise" the paper refers to
  # is the one previously mentioned for GL
  additive_noise = randn(n, 1, mean = 0, sd = 1)
  
  b = A%*%w0 + additive_noise
  
  return(list("G" = cal_G, "A" = A, "b" = b, "w0" = w0))
}

#### LMO as described in the Sub-sampling section in
#### Kedreux et al.
#### n_groups in this case is the cardinality of Gp
LMO = function(A, G, G_p, x_t, grad, n_groups, group_size = 10, overlap_size = 3, gp_size = 14, gp_first_covariate)
{
  n = dim(A)[1]
  d = dim(A)[2]
  s = zeros(d) # LMO vertex
  
  # Instead of iterating over ALL groups as we have done for FW, we have to 
  # iterate over the Gp groups and compute the euclidean norm on them
  
  Gp_gradients = zeros(gp_size)
  for(i in seq(1:gp_size))
  {
    # Given the group index at position in in Gp, we want
    # to retrieve all its covariates associated from the 
    # original group..
    group_covariate_indices = G[[G_p[i]]]
    
    # But, of course, the index is absolute right now. We want to map the number
    # of the indices with respect to the covariates related with gp
    # Computing the euclidean norm for each g in Gp. For example:
    # Gp's covariates starts from 300 to 450. This means that our Gp's gradient
    # has a length equal to 50, but its fist element refers to 300, the 2nd to
    # 301 and so on. Thus, we need to obtain the offset, by taking the overlap size
    # in account.
    
    relative_indices = (group_covariate_indices - gp_first_covariate) + 1 #+1 because R starts indexing from 1.
    Gp_gradients[i] = sqrt(sum(grad[relative_indices]^2))
  }
  
  lmo_index = which.max(Gp_gradients)
  
  # For the experiment section, we have to store the number of 
  # gradients computed at each iteration. This can be obtained
  # by the following math:
  # n_coef_involved = group_size * |Gp| - (|Gp| * overlap_size)
  n_coef_involved = (gp_size *  group_size) - ((gp_size-1) * overlap_size) # -1 because we have gp_size-1 intersections
  
  # Again, we retrieve the relative indices, but at this time for the
  # group with the maximum gradient value
  max_group_covariate_indices = G[[G_p[lmo_index]]]
  max_group_covariates = (max_group_covariate_indices - gp_first_covariate) + 1
  s[max_group_covariate_indices] = 14 * grad[max_group_covariates] / Gp_gradients[lmo_index]
  # to add: dual gap
  return(list("s" = s, "n_coef" = n_coef_involved))
}

#### Gradient, starting from the equation 25 in section 7
complete_gradient = function(A, b, x_t, residual)
{
  # As in the Lasso case, the derivative is similar: -X^t(Y-X*x_t)
  # for each group, so it will be multiplied by
  
  n = dim(A)[1]
  d = dim(A)[2]
  
  g = t(A)%*%(residual)
  
  return(g)
}

#### Gradient, starting from the equation 25 in section 7
gradient = function(A, b, x_t, residual, first_covariate, last_covariate)
{
  # As in the Lasso case, the derivative is similar: -X^t(Y-X*x_t)
  # for each group, so it will be multiplied by
  
  n = dim(A)[1]
  d = dim(A)[2]
  
  # Subgradient.
  g = t(A[,first_covariate:last_covariate])%*%(residual)
  
  return(g)
}

# It this the same as the one in FW, for which we perform only a 
# subselection of the G_p groups after having computed all the
# L2 norms
complete_LMO = function(A, G, G_p, p, x_t, grad, gp_first_group)
{
  n = dim(A)[1]
  d = dim(A)[2]
  n_groups = length(G)
  s = zeros(d) # LMO vertex
  
  groups_gradient_regularized = zeros(n_groups)
  
  for(i in seq(1:n_groups))
  {
    groups_gradient_regularized[i] = sqrt(sum(grad[G[[i]]]^2))
  }
  # We are only interested in the groups obtained from the subsampling.
  regularized_grad_Gp = groups_gradient_regularized[G_p]
  
  lmo_index = which.max(groups_gradient_regularized)
  lmo_gp_index = which.max(regularized_grad_Gp)
  
  # Checking for boundaries:
  lmo_gp_index = min(gp_first_group + lmo_gp_index, length(G))
  lmo_group = G[[lmo_gp_index]]

  s[lmo_group] = 14 * grad[lmo_group] / groups_gradient_regularized[lmo_gp_index]
  dual_gap = 14 * groups_gradient_regularized[lmo_index] + -t(grad) %*% x_t
  
  return(list("s" = s, "gap" = dual_gap))
}

#### Line search
line_search = function(A, b, s, residual)
{
  n = dim(A)[1]
  d = dim(A)[2]
  
  s_v = which(s != 0) # Getting the active coefficients indices from the LMO vertex
  first_g_index = s_v[1] # First and last coefficients
  last_g_index = s_v[length(s_v)]
  
  As = A[, first_g_index:last_g_index] %*% s[first_g_index:last_g_index]
  
  # This is equal to As - b + b - Ax = As-Ax
  d = As - b + residual
  
  gamma_max = (t(residual) %*% d) / (t(d) %*% d)
  gamma_t = max(0, min(1, gamma_max)) # It's the clip described in the paper
  
  return(gamma_t)
}

#### RFW Algorithm
RFW = function(G, A, b, max_iterations, tollerance, x_star, group_size = 10, overlap_size = 3, eta = 0.1, verbose = FALSE)
{
  n = dim(A)[1]
  d = dim(A)[2]
  
  # Return values we are interested for plots.
  of_values = zeros(max_iterations)
  gaps = zeros(max_iterations)
  size_support = zeros(max_iterations)
  recovered_coeff = zeros(max_iterations)
  iteration_times = zeros(max_iterations)
  
  x_t = zeros(d)
   
  residual = b - A%*%x_t

  
  
  for(t in seq(1:max_iterations))
  {
    tic = Sys.time()
    # p refers to the index of the subsampling of the papers, as descrieb in Kedreux et al.
    p = floor(eta*length(G))
    
    "
    Again in the paper, talking about the LMO in RFW, it is also explained that
    this means that we only need to compute the gradient on the gp index.
    In practice, we have to compute a subgradient by picking at random a number of
    groups equal to the eta percentage and then compute the gradient over it.
  "
    random_group_index = sample(1:length(G), 1)
    
    # After that, we select the p consecutive groups, creating our G_p
    # described in the paper
    last_random_group_index = random_group_index + p
    G_p = random_group_index : min(length(G), last_random_group_index) # Because we can go out of bounds
    
    # Once we have the subsample of groups, we have to retrieve its covariates index to compute
    # the subgradient. Since groups are consecutive, we pick the first one in the first group and
    # the last one in the last, and then we create a range
    
    first_group_first_covariate = G[[G_p[1]]][1]
    last_group_last_covariate = G[[G_p[length(G_p)]]][group_size] # the 10th element is the last.
    
    subsample_covariates = seq(first_group_first_covariate:last_group_last_covariate)
    
    # Gradient
    g = gradient(A, b, x_t, residual, first_group_first_covariate, last_group_last_covariate)
    
    # LMO
    lmo_result = LMO(A, G, G_p, x_t, g, group_size = group_size, overlap_size, gp_size = p, gp_first_covariate = first_group_first_covariate)
    s = lmo_result$s
    n_coef = lmo_result$n_coef
    
    gamma_t = line_search(A, b, s, residual)
    
    x_t = (1-gamma_t)*x_t + gamma_t*s
    
    # Update residual
    residual = gamma_t*b + (1-gamma_t)*residual - gamma_t*A%*%s
    
    toc = Sys.time()
    
    of_values[t] = sum(residual^2)
    iteration_times[t] = as.numeric(toc - tic)
    recovered_coeff[t] = lmo_result$n_coef
    
    # Computing full gradient and gap in order to retrieve a measure
    # for the comparison plots. It will not be measured
    g = complete_gradient(A, b, x_t, residual)
    s = lmo_result$s
    
    lmo_result = complete_LMO(A, G, G_p, p, x_t, g, random_group_index)
    gap = lmo_result$gap
    gaps[t] = lmo_result$gap
    
    if(verbose == TRUE)
    {
      cat("t = ", t, "loss = ", of_values[t], "gap = ", gaps[t], "step = ", gamma_t, "time" = iteration_times[t], "\n\n")
    }
  }
  
  return(list("time" = iteration_times, "loss" = of_values, "gaps" = gaps, "G" = G, "n_coef" = recovered_coeff))
}

#### Main
setwd("C:\\Users\\ettag\\Documents\\GitHub\\Optimization-Project")
d = 10000
dataset = generate_datset(n = 1000, d = d, overlap_size = 3, group_size = 10)

G = dataset$G
A = dataset$A
b = dataset$b
w0 = dataset$w0


max_iterations = 500 
eta = 0.1 # As described in the experiments section

result = RFW(G = G, A = A, b = b, x_star = w0, tollerance = 1e-6,  max_iterations = 500, group_size = 10, overlap_size = 3, eta = 0.1, verbose = T)

df_rfw_gaps_by_time = tibble(time = result$time, gap = result$gaps)
df_rfw_gaps_by_time$time = cumsum(df_rfw_gaps_by_time$time)

df_rfw_gaps_by_coef = tibble(nbr_coef_per_grad = result$n_coef, gap = result$gaps)
df_rfw_gaps_by_coef$nbr_coef_per_grad = cumsum(df_rfw_gaps_by_coef$nbr_coef_per_grad)

save(df_rfw_gaps_by_coef, file = ".\\data\\LGL\\RFW_gaps_by_coef.RData")
save(df_rfw_gaps_by_time, file = ".\\data\\LGL\\RFW_gaps_by_time.RData")

  