
rm(list = ls())

library(ramify) # For randn
library(ggplot2) # For plotting
library(tibble) # For dataframe creation
library(data.table) # For data table manipulation s
library(scales) # For trans_breaks() in plots
library(Matrix) # For sparse initial point
library(dplyr) # For in place data.table manipulations


#### Synthetic dataset creation
# For the extended description about the creation of the 
# synthetic dataset, please refer to:
# http://papers.nips.cc/paper/5925-on-the-global-linear-convergence-of-frank-wolfe-optimization-variants.pdf
create_synthetic_dataset = function(m, n)
{
  set.seed(42)
  A = randn(m, n, mean = 0, sd = 1)
  
  # x_star is a vector with 10% of nonzero coefficients and 
  # values in {−1, +1}
  # The 10% of 500 is 50, which is divided by 2 since our atoms
  # in the dataset are specular.
  non_zero_coef_indices = sample(seq(1:500), 51)
  
  x_star = zeros(n)
  x_star[non_zero_coef_indices[1:25]] = 1
  x_star[non_zero_coef_indices[26:50]] = -1
  
  # 10% of additive noise:
  additive_noise = randn(m, 1, mean = 0, sd = 1)
  
  b = A%*%x_star + additive_noise
  
  # In the L1 norm, we have a symmetric simplex,
  # thus we have to add the counterpoints to the 
  # dataset
  # 40 is the regularization parameter
  A = 40 * cbind(A, -A)
  
  return(list("A" = A, "b" = b, "x_star" = x_star))
}

#### Objective Function cost
function_cost = function(Ax,b)
{
  cost = sum((Ax-b)^2)/2
  return(cost)
}

gradient = function(Ax, b)
{
  # Because the only time that we need the matrix A itself
  # is when we calculate the LMO on the entire convex set.
  # For the computation of the direction, we will compute s-x,
  # where s and x are retrieved from A multiplied with their 
  # corresponding canonical basis.
  # For the gap computation, instead, we will have something like:
  # d*grad, which in our case will be:
  # (s-x)*(Ax-b) = (A^t*ej - A^T*ek)(Ax-b) = (ej-ek)A^t(Ax-b),
  # where A^t(Ax-b) is the canonical gradient.
  gradient = (Ax - b)
  
  return(gradient)
}

LMO = function(A, gradient)
{
  complete_gradient = t(A) %*% gradient
  lmo_index = which.max(abs(complete_gradient))
  sign_index = sign(-complete_gradient[lmo_index])
  return(list("lmo_index" = lmo_index, "sign" = sign_index))
}

away_step_calculation = function(A, gradient, S_t)
{
  # Calculating the away step. It is called "absolute"
  # since we have no information in the result about the
  # atom index they refer to.
  complete_gradient = t(A[, S_t]) %*% gradient
  vt_index_absolute = which.max(complete_gradient)
  sign_index = sign(complete_gradient[vt_index_absolute])
  
  # The previous calculation returns as much rows as the 
  # atoms in the active set, but we have no clue about who
  # they refer to. Therefore, we proceed to perform a "mapping
  # step" on the active set S_t in order to obtain the original 
  # index.
  vt_index = S_t[vt_index_absolute]
  
  return(list("vt_index" = vt_index, "sign" = sign_index))
}

PFW = function(A, b, max_iterations, tollerance, x_star, verbose = FALSE, synthetic = FALSE)
{
  m = dim(A)[1]
  n = dim(A)[2]
  
  # Initialization, by following the paper (https://arxiv.org/pdf/1511.05932.pdf)
  # and its explanation of PFW since FW with subsampling oracle cites it
  
  # coefficients for our solution
  alpha_t = zeros(n)
  
  if(synthetic == FALSE)
  {
    alpha_t[5192] = 1
  } else {
    alpha_t[1] = 1
  }
  
  # Initial atom x_0
  x_t = NULL
  
  if(synthetic == FALSE)
  {
    x_t = A[, 5192]
  } else {
    x_t = A[, 1]  # It is equal as multiplying A * e1.
  }
  
  # List of active atoms S_t. It contains the atoms
  # with an associated weight > 0 (see line 13 of paper
  # http://papers.nips.cc/paper/5925-on-the-global-linear-convergence-of-frank-wolfe-optimization-variants.pdf). 
  # Therefore, at the begin it is only composed of the first atom.
  S_t = which(alpha_t > 0)
  
  # Return values we are interested for plots.
  of_values = zeros(max_iterations)
  gaps = zeros(max_iterations)
  size_support = zeros(max_iterations)
  recovered_coeff = zeros(max_iterations)
  
  x_star_indices = which(x_star != 0)
  for(t in seq(1:max_iterations))
  {
    # Calculation the cost function
    of_values[t] = function_cost(x_t, b)
    
    # Calculating the gradient:
    gradient = gradient(x_t, b)
    
    # Calculating the LMO for our direction.
    lmo_result = LMO(A, gradient)
    index_s = lmo_result$lmo_index
    sign = lmo_result$sign
    
    s_t = (sign) * A[, index_s]
    
    # Calculation of the duality gap.
    # It follows the same thinking explained in the gradient
    # function
    d = s_t - x_t
    
    gap = - d %*% gradient
    gaps[t] = gap[1]
    
    # Away step calculation:
    away_step = away_step_calculation(A, gradient, S_t)
    index_vt = away_step$vt_index
    sign = away_step$sign
    
    v_t = (sign) * A[, index_vt]
    
    # Calculating the direction. 
    # We remember that As = s and Ax = x_t,
    # it's just to remember that they come from a
    # matrix A multiplied by a canonical basis.
    d_t = s_t - v_t 
    
    # Replacing as explained in (https://arxiv.org/pdf/1511.05932.pdf)
    alpha_vt = alpha_t[index_vt]
    gamma_max = alpha_vt
    
    if(gap < tollerance)
    {
      print("Gap is lower than the tollerance. Exiting..")
      break
    }
    
    # Line search
    # We have 3 options: fixed stepsize, find it with exact
    # line search or use armijo. The paper (https://arxiv.org/pdf/1511.05932.pdf)
    # proposes the exact line search, so we go for it.
    step = - (t(gradient) %*% d_t) / (t(d_t) %*% d_t)
    gamma_t = max(0, min(gamma_max, step)) # It's the clip described in the paper
    
    if(is.nan(gamma_t))
    {
      print("Stepsize is too low. Quitting..")
      return(list("gaps" = gaps, "f_costs" = of_values, "size_support" = size_support, "recovered_coeff" = recovered_coeff))
    }
    # Managing the so called "bad drop step" explained at page
    # 4 of paper (http://papers.nips.cc/paper/5925-on-the-global-linear-convergence-of-frank-wolfe-optimization-variants.pdf)
    if(gamma_t == gamma_max)
    {
      print("Drop step")
      alpha_t[index_vt] = 0
      S_t = S_t[S_t != index_vt]
    } 
    else 
    {
      alpha_t[index_vt] = alpha_t[index_vt] - step
    }
    
    # Updating weights as explained at page 4 of (https://arxiv.org/pdf/1511.05932.pdf)
    alpha_t[index_s] = alpha_t[index_s] + step
    
    # Looking for recovered coefficients
    alpha_t_nonzeros = which(alpha_t != 0)
    recovered_coeff[t] = length(base::intersect(alpha_t_nonzeros, x_star_indices))
    
    # Cut off threshold for size support
    alpha_t[alpha_t < 1e-3] = 0
    size_support[t] = length(alpha_t[alpha_t > 0])
    
    S_t = which(alpha_t > 0)
    
    # Moving onto the feasible direction to find the new point
    x_t = x_t + gamma_t * d_t
    
    
    if(verbose == TRUE)
    {
      cat("t = ", t, "loss = ", of_values[t], "gap = ", gap, "step = ", step, "\n\n")
    }
  }
  return(list("gaps" = gaps, "f_costs" = of_values, "size_support" = size_support, "recovered_coeff" = recovered_coeff))
}

#### Main
setwd("C:\\Users\\ettag\\Documents\\GitHub\\Optimization-Project")
tollerance = 1e-12
max_iterations = 20000
m = 200
n = 500

synthetic_dataset = T
datasets = NULL
A = NULL
b = NULL
x_star = c()

if(synthetic_dataset == TRUE)
{
  datasets = create_synthetic_dataset(m, n)
  
  A = datasets$A
  b = datasets$b
  x_star = datasets$x_star
} else {
  max_iterations = 150
  A = read.csv(".\\data\\Lasso\\datasets\\E2006_merged.csv", sep = ",")
  A = A[1:8000, 2:8000]
  b = read.csv(".\\data\\Lasso\\datasets\\E2006_target.csv", sep = ",")
  b = b[1:8000, 2]
}

print("Executing Pairwise Frank Wolfe..")
result = PFW(A, b, max_iterations, tollerance, x_star, verbose = TRUE, synthetic = synthetic_dataset)
print("Done")

gaps = result$gaps
f_costs = result$f_costs
size_support = result$size_support
recovered_coeff = result$recovered_coeff

df_pfw_gaps = tibble(nbr_coef_per_grad = n*2 * ones(max_iterations), gap = gaps[1:max_iterations])
df_pfw_gaps$nbr_coef_per_grad[1] = 1
df_pfw_gaps$nbr_coef_per_grad = cumsum(df_pfw_gaps$nbr_coef_per_grad)

df_pfw_gaps_by_iterations = tibble(iteration = seq(1:max_iterations), gap = gaps[1:max_iterations])
df_pfw_costs = tibble(iteration = seq(1:max_iterations), function_cost = f_costs[1:max_iterations])
df_pfw_size_support = tibble(iteration = seq(1:max_iterations), size_support = size_support)
df_pfw_recovered_coeff = tibble(iteration = seq(1:max_iterations), recovered_coeff = recovered_coeff)

if(synthetic_dataset == TRUE)
{
  save(df_pfw_gaps, file = ".\\data\\Lasso\\PFW_gaps.RData")
  save(df_pfw_gaps_by_iterations, file = ".\\data\\Lasso\\PFW_gaps_by_iter.RData")
  save(df_pfw_costs, file = ".\\data\\Lasso\\PFW_function_costs.RData")
  save(df_pfw_size_support, file = ".\\data\\Lasso\\PFW_support_size.RData")
  save(df_pfw_recovered_coeff, file = ".\\data\\Lasso\\PFW_recovered_coeff.RData")
} else {
  save(df_pfw_gaps, file = ".\\data\\Lasso\\PFW_gaps_real.RData")
  save(df_pfw_gaps_by_iterations, file = ".\\data\\Lasso\\PFW_gaps_by_iter_real.RData")
  save(df_pfw_costs, file = ".\\data\\Lasso\\PFW_function_costs_real.RData")
  save(df_pfw_size_support, file = ".\\data\\Lasso\\PFW_support_size_real.RData")
  save(df_pfw_recovered_coeff, file = ".\\data\\Lasso\\PFW_recovered_coeff_real.RData")

}