library(ramify) # For randn
library(tictoc) # For time tracking

rm(list = ls())

#### Dataset creation according to Kedreux et al. paper
generate_dataset = function(n, d, overlap_size, group_size)
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
  
  # From the paper: We chose the ground truth parameter vector w0 ∈ conv(A) with a fraction of 0.01 of 
  # nonzero coefficients, where on eachactive group, the coefficients are generated 
  # from a Gaussian distribution
  w0 = zeros(d)
  g0 = c() 
  
  # Selecting groups with nonzero coef for the solution
  w0_indices = sample(1 : length(cal_G), gt_nonzero) # 0.01 of d
  
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
  
    # Matrix A of GL in the Obozinski literature is referred as X,
    # but the core remains the same. The combination of A with the
    # management of indices that will result in the opertor \obigtimes
    # which is the covariate duplicator operator.
    A = randn(n, d, mean = 0, sd = 1)
    
    # We suppose that the "some additive noise" the paper refers to
    # is the one previously mentioned for GL
    additive_noise = randn(n, 1, mean = 0, sd = 1)
    
    b = A%*%w0 + additive_noise
    
    return(list("G" = cal_G, "A" = A, "b" = b, "w0" = w0))
}

#### LMO as described in the Sub-sampling section in
#### Kedreux et al.
LMO = function(A, G, x_t, grad)
{
  n = dim(A)[1]
  d = dim(A)[2]
  n_groups = length(G)
  s = zeros(d) # LMO vertex
  
  # According to the solution for group lasso proposed in Tibshirani,
  # we compute the euclidean norm of the groups gradient in advance
  # in order to establish the LMO.
  groups_gradient_regularized = zeros(n_groups)
  
  for(i in seq(1:n_groups))
  {
    groups_gradient_regularized[i] = sqrt(sum(grad[G[[i]]]^2))
  }
  
  # Again, according to the Tibshirani finds the LMO,
  # to find the minimum we only have to maximize the denominator, i.e.,
  # the euclidean norm of the groups gradient
  lmo_index = which.max(groups_gradient_regularized)
  lmo_group = G[[lmo_index]]
  
  # LMO vertex
  s[lmo_group] = 14 * grad[lmo_group] / groups_gradient_regularized[lmo_index]
  
  # If we perform the canonical dual gap for FW, i.e., -grad*(s-x_t), we did not
  # obtain a good gap. Thus, we thought about adding the gradient of the group
  # corresponding the LMO index, because we define it as the contribution provided
  # by the Group Lasso regularization term for the LMO vertex x and its corresponding
  # active elements.
  dual_gap = 14 * + groups_gradient_regularized[lmo_index] -t(grad) %*% x_t
  
  # We added the regularization term on both because otherwise it did not converge
  return(list("s" = s, "gap" = dual_gap))
}

# Gradient of the empirical risk function, which is the LSE in our case.
gradient = function(A, b, x_t, residual)
{
  n = dim(A)[1]
  d = dim(A)[2]
  
  g = t(A)%*%(residual)
  
  return(g)
}

#### Line search
line_search = function(A, b, s, residual)
{
  n = dim(A)[1]
  d = dim(A)[2]
  
  s_v = which(s != 0) # We obtain the indices corresponding to the active elements of s.
  first_g_index = s_v[1] # First and last column indices for A.
  last_g_index = s_v[length(s_v)]
  
  As = A[, first_g_index:last_g_index] %*% s[first_g_index:last_g_index]
  
  # Since we defined the residual as b - Ax, this is an overengineering
  # to obtain As - Ax = s - x, the direction.
  d = As - b + residual
  
  gamma_max = (t(residual) %*% d) / (t(d) %*% d)
  gamma_t = max(0, min(1, gamma_max)) # It's the clip described in the paper
  
  return(gamma_t)
}

#### FW Algorithm
FW = function(G, A, b, max_iterations, tollerance, x_star, verbose = FALSE)
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
    g = gradient(A, b, x_t, residual)
    
    lmo_result = LMO(A, G, x_t, g)
    gap = lmo_result$gap
    s = lmo_result$s
    
    gamma_t = line_search(A, b, s, residual)
    
    x_t = (1-gamma_t)*x_t + gamma_t*s
    
    # Update residual
    residual = gamma_t*b + (1-gamma_t)*residual - gamma_t*A%*%s
    
    toc = Sys.time()
    
    gaps[t] = lmo_result$gap
    of_values[t] = sum(residual^2)
    iteration_times[t] = as.numeric(toc - tic)
    
    if(verbose == TRUE)
    {
      cat("t = ", t, "loss = ", of_values[t], "gap = ", gaps[t], "step = ", gamma_t, "time" = iteration_times[t], "\n\n")
    }
  }
  
  return(list("time" = iteration_times, "loss" = of_values, "gaps" = gaps, "G" = G))
}

#### Main
setwd("C:\\Users\\ettag\\Documents\\GitHub\\Optimization-Project")
d = 10000
dataset = generate_dataset(n = 1000, d = d, overlap_size = 3, group_size = 10)

G = dataset$G
A = dataset$A
b = dataset$b
w0 = dataset$w0

max_iterations = 500

result = FW(G = G, A = A, b = b, x_star = w0, tollerance = 1e-6, max_iterations = 500, verbose = T)

df_fw_gaps_by_time = tibble(time = result$time, gap = result$gaps)
df_fw_gaps_by_time$time = cumsum(df_fw_gaps_by_time$time)

df_fw_gaps_by_coef = tibble(nbr_coef_per_grad = d*2 * ones(max_iterations), gap = result$gap)
df_fw_gaps_by_coef$nbr_coef_per_grad = cumsum(df_fw_gaps_by_coef$nbr_coef_per_grad)

save(df_fw_gaps_by_coef, file = ".\\data\\LGL\\FW_gaps_by_coef.RData")
save(df_fw_gaps_by_time, file = ".\\data\\LGL\\FW_gaps_by_time.RData")

