library(ggplot2) # For plotting
library(tibble) # For dataframe creation
library(data.table) # For data table manipulation
library(scales) # For trans_breaks() in plots


rm(list = ls())

#### Loading data
setwd("C:\\Users\\ettag\\Documents\\GitHub\\Optimization-Project")

load(".\\data\\Lasso\\FW_gaps_real.RData")
load(".\\data\\Lasso\\FW_gaps_by_iter_real.RData")
load(".\\data\\Lasso\\FW_function_costs_real.RData")
load(".\\data\\Lasso\\FW_support_size_real.RData")
load(".\\data\\Lasso\\FW_recovered_coeff_real.RData")

load(".\\data\\Lasso\\RFW_gaps_real.RData")
load(".\\data\\Lasso\\RFW_gaps_by_iter_real.RData")
load(".\\data\\Lasso\\RFW_function_costs_real.RData")
load(".\\data\\Lasso\\RFW_support_size_real.RData")
load(".\\data\\Lasso\\RFW_recovered_coeff_real.RData")



# TODO: Move then into dataframe creation
df_fw_gaps$Method = "FW"
df_fw_costs$Method = "FW"
df_fw_size_support$Method = "FW"
df_fw_recovered_coeff$Method = "FW"
df_fw_gaps_by_iterations$Method = "FW"

df_rfw_gaps$Method = "RFW"
df_rfw_costs$Method = "RFW"
df_rfw_size_support$Method = "RFW"
df_rfw_recovered_coeff$Method = "RFW"
df_rfw_gaps_by_iterations$Method = "RFW"


df_gaps = rbind(df_fw_gaps, df_rfw_gaps)
df_objective_functions = rbind(df_fw_costs, df_rfw_costs)
df_size_support = rbind(df_fw_size_support, df_rfw_size_support)
df_recovered_coeff = rbind(df_fw_recovered_coeff, df_rfw_recovered_coeff)
df_gaps_by_iter = rbind(df_fw_gaps_by_iterations, df_rfw_gaps_by_iterations)

#### Plotting gaps over iterations
plot_gaps_by_coef = ggplot(df_gaps, aes(x = nbr_coef_per_grad, y = gap, color = Method)) +
  geom_line() +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))


#### Plotting gap over iterations
plot_gaps_by_iter = ggplot(df_gaps_by_iter, aes(x = iteration, y = gap, color = Method)) +
  geom_line() +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))
#### Plotting function costs over iterations
plot_of_by_niterations = ggplot(df_objective_functions, aes(x = iteration, y = function_cost, color = Method)) +
  geom_line() +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))


#### Plotting the size of the support over iterations
plot_size_support_by_iterations = ggplot(df_size_support, aes(x = iteration, y = size_support, color = Method)) +
  geom_line() +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))

# Plotting recovered coeff
recovered_coeff_by_iterations = ggplot(df_recovered_coeff, aes(x = iteration, y = recovered_coeff, color = Method)) +
  geom_line() +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))

library(patchwork)

(plot_gaps_by_iter / plot_gaps_by_coef)
