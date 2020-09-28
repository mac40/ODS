
rm(list = ls())

library(ramify) # For randn
library(ggplot2) # For plotting
library(tibble) # For dataframe creation
library(data.table) # For data table manipulation s
library(scales) # For trans_breaks() in plots


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
  # 
  # 
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
  # We do not calculate the canonical gradient A^T(Ax - b) 
  # Why?
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

RFW_LMO = function(A, gradient, subsample_indices)
{
  local_gradient = t(A[, subsample_indices]) %*% gradient
  lmo_cal_At_Absolute = which.max(abs(local_gradient)) # We could have used max(abs(A^t*g))
  lmo_cal_At = subsample_indices[lmo_cal_At_Absolute]
  sign_index = sign(-local_gradient[lmo_cal_At_Absolute])
  return(list("lmo_cal_At" = lmo_cal_At, "sign" = sign_index))
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

RAFW = function(A, b, max_iterations, tollerance, eta, x_star, t_0, verbose = FALSE, synthetic = FALSE)
{
  m = dim(A)[1]
  n = dim(A)[2]
  
  # Subsample size of our Calligraphic A_t
  subsample_size = n * eta
  
  # Initialization, by following the paper (https://arxiv.org/pdf/1511.05932.pdf)
  # and its explanation of AFW since FW with subsampling oracle cites it
  
  # Weights of the covariates
  alpha_t = zeros(n)
  

  cat("T_0: ", t_0, "\n\n")
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
    x_t = A[1] # It is equal as multiplying A * e1.
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
  drop_steps = 0
  away_steps = 0
  fw_steps = 0
  
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
    
    # Away step calculation:
    away_step = away_step_calculation(A, gradient, S_t)
    index_vt = away_step$vt_index
    sign = away_step$sign
    
    v_t = (sign) * A[, index_vt]
  
    # Calculating the LMO for our subsample Calligraphic A_t.
    # We first need to extract indices with the same probability
    cal_At_indices = sample(seq(1:n), subsample_size)
    rfw_lmo_result = RFW_LMO(A, gradient, cal_At_indices)
    
    # In the paper FW with subsampling oracle, for RAFW we have
    # a set of atoms corresponding of the union between the active
    # ones and the others obtained from subsampling, i.e., A_t U S_t
    St_U_At = union(S_t, cal_At_indices)
    index_cal_At = rfw_lmo_result$lmo_cal_At
    sign = rfw_lmo_result$sign
    s_t_rfw = sign * A[, index_cal_At]
    
    # Calculating the direction. 
    # We remember that As = s and Ax = x_t,
    # it's just to remember that they come from a
    # matrix A multiplied by a canonical basis.
    d_t = s_t - x_t 
    
    # RFW direction
    d_rfw_t = s_t_rfw - x_t
    
    # AS direction
    d_as_t = x_t - v_t
    
    # Calculation of the duality gap.
    # It follows the same thinking explained in the gradient
    # function
    gap_rfw = - d_rfw_t %*% gradient
    gap_as = -d_as_t %*% gradient
    gap = -d_t %*% gradient
    gaps[t] = gap
    
    d = NULL
    gamma_max = NULL
    direction_chosen = NULL
    
    if(gap_rfw >= gap_as)
    {
      # RFW case
      d = d_rfw_t
      direction_chosen = "RFW"
      gamma_max = 1
    }
    else
    {
      # AS case
      # The alpha_vt in algorithm the paper (http://papers.nips.cc/paper/5925-on-the-global-linear-convergence-of-frank-wolfe-optimization-variants.pdf)
      # is the weight corresponding to the atom chosen in the AS
      alpha_vt = alpha_t[index_vt]
      d = d_as_t
      direction_chosen = "AS"
      gamma_max = alpha_vt / (1-alpha_vt)
    }

    # Line search
    # We have 3 options: fixed stepsize, find it with exact
    # line search or use armijo. The paper (https://arxiv.org/pdf/1511.05932.pdf)
    # proposes the exact line search, so we go for it.
    step = - (t(gradient) %*% d) / (t(d) %*% d)
    gamma_t = max(0, min(gamma_max, step)) # It's the clip described in the paper
    
    
    # Updating weights as explained in pag 3 of the paper (http://papers.nips.cc/paper/5925-on-the-global-linear-convergence-of-frank-wolfe-optimization-variants.pdf)
    
    if(direction_chosen == "RFW")
    {
      alpha_t = (1 - gamma_t) * alpha_t
      alpha_t[index_cal_At] = alpha_t[index_cal_At] + gamma_t
      fw_steps = fw_steps + 1
    }
    else
    {
      alpha_t = (1 + gamma_t) * alpha_t
      
      # Managing the so called "bad drop step" explained at page
      # 4 of paper (http://papers.nips.cc/paper/5925-on-the-global-linear-convergence-of-frank-wolfe-optimization-variants.pdf)
      if(gamma_t == gamma_max)
      {
        print("Drop step!")
        alpha_t[index_vt] = 0
        S_t = S_t[S_t != index_vt]
        drop_steps = drop_steps + 1 
      }
      else
      {
        away_steps = away_steps + 1
        alpha_t[index_vt] = alpha_t[index_vt] - gamma_t
      }
    }
    
    # Looking for recovered coefficients
    alpha_t_nonzeros = which(alpha_t != 0)
    recovered_coeff[t] = length(base::intersect(alpha_t_nonzeros, x_star_indices))
    
    # # Cut off threshold for size support
    alpha_t[alpha_t < 1e-3] = 0
    
    # 
    # List of active atoms.
    S_t = which(alpha_t > 0)
    size_support[t] = length(S_t)
    # Moving onto the feasible direction to find the new point
    x_t = x_t + gamma_t * d
    
    if(verbose == TRUE)
    {
      cat("t = ", t, "loss = ", of_values[t], "gap = ", gap, "step = ", gamma_t, "\n\n")
    }
    
    if(gap < tollerance)
    {
      print("Gap is lower than the tollerance. Exiting..")
      return(list("gaps" = NULL, "f_costs" = NULL, "size_support" = NULL, "recovered_coeff" = NULL))
    }
  }
  
  return(list("gaps" = gaps, "f_costs" = of_values, "size_support" = size_support, "recovered_coeff" = recovered_coeff, "drop_steps" = drop_steps,
              "away_steps" = away_steps, "fw_steps" = fw_steps))
}

#### Main
setwd("C:\\Users\\ettag\\Documents\\GitHub\\Optimization-Project")
tollerance = 1e-12
max_iterations = 20000
m = 200
n = 500
eta = 0.25

synthetic_dataset = T
datasets = NULL
A = NULL
b = NULL
x_star = c()

if(synthetic_dataset == TRUE)
{
  max_iterations = 20000
  eta = 0.03
  datasets = create_synthetic_dataset(m, n)
  
  A = datasets$A
  b = datasets$b
  x_star = datasets$x_star
} else {
  max_iterations = 150
  eta = 0.25
  A = read.csv(".\\data\\Lasso\\datasets\\E2006_merged.csv", sep = ",")
  A = 40 * cbind(A, -A)
  b = read.csv(".\\data\\Lasso\\datasets\\E2006_target.csv", sep = ",")
  b = b[1:8000, 2]
}

result = NULL

print("Executing Frank Wolfe..")
if(synthetic_dataset == FALSE)
{
  result = RAFW(A, b, max_iterations, tollerance, eta, x_star, t_0 = 5192, verbose = TRUE, synthetic_dataset)
} else {
  result = RAFW(A, b, max_iterations, tollerance, eta, x_star, t_0 = 1, verbose = TRUE, synthetic_dataset)
}

print("Done")

gaps = result$gaps
f_costs = result$f_costs
size_support = result$size_support
recovered_coeff = result$recovered_coeff

df_rafw_gaps = tibble(nbr_coef_per_grad = eta*(n*2) * ones(max_iterations), gap = gaps[1:max_iterations])
df_rafw_gaps$nbr_coef_per_grad[1] = 1
df_rafw_gaps$nbr_coef_per_grad = cumsum(df_rafw_gaps$nbr_coef_per_grad)

df_rafw_gaps_by_iterations = tibble(iteration = seq(1:max_iterations), gap = gaps[1:max_iterations])
df_rafw_costs = tibble(iteration = seq(1:max_iterations), function_cost = f_costs[1:max_iterations])
df_rafw_size_support = tibble(iteration = seq(1:max_iterations), size_support = size_support)
df_rafw_recovered_coeff = tibble(iteration = seq(1:max_iterations), recovered_coeff = recovered_coeff)
df_rafw_fw_steps = tibble(step_value = result$fw_steps)
df_rafw_drop_steps = tibble(step_value = result$drop_steps)
df_rafw_away_steps = tibble(step_value = result$away_steps)

if(synthetic_dataset == TRUE)
{
  save(df_rafw_gaps, file = ".\\data\\Lasso\\RAFW_gaps.RData")
  save(df_rafw_gaps_by_iterations, file = ".\\data\\Lasso\\RAFW_gaps_by_iter.RData")
  save(df_rafw_costs, file = ".\\data\\Lasso\\RAFW_function_costs.RData")
  save(df_rafw_size_support, file = ".\\data\\Lasso\\RAFW_support_size.RData")
  save(df_rafw_recovered_coeff, file = ".\\data\\Lasso\\RAFW_recovered_coeff.RData")
  save(df_rafw_fw_steps, file = ".\\data\\Lasso\\RAFW_FW_steps.RData")
  save(df_rafw_drop_steps, file = ".\\data\\Lasso\\RAFW_drop_steps.RData")
  save(df_rafw_away_steps, file = ".\\data\\Lasso\\RAFW_away_steps.RData")
} else {
  save(df_rafw_gaps, file = ".\\data\\Lasso\\RAFW_gaps_real.RData")
  save(df_rafw_gaps_by_iterations, file = ".\\data\\Lasso\\RAFW_gaps_by_iter_real.RData")
  save(df_rafw_costs, file = ".\\data\\Lasso\\RAFW_function_costs_real.RData")
  save(df_rafw_size_support, file = ".\\data\\Lasso\\RAFW_support_size_real.RData")
  save(df_rafw_recovered_coeff, file = ".\\data\\Lasso\\RAFW_recovered_coeff_real.RData")
  save(df_rafw_fw_steps, file = ".\\data\\Lasso\\RAFW_FW_steps_real.RData")
  save(df_rafw_drop_steps, file = ".\\data\\Lasso\\RAFW_drop_steps_real.RData")
  save(df_rafw_away_steps, file = ".\\data\\Lasso\\RAFW_away_steps_real.RData")
  }