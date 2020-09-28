library(ggplot2) # For plotting
library(tibble) # For dataframe creation
library(data.table) # For data table manipulation
library(scales) # For trans_breaks() in plots


rm(list = ls())

#### Loading data
setwd("C:\\Users\\ettag\\Documents\\GitHub\\Optimization-Project")
load(".\\data\\Lasso\\RAFW_gaps_real.RData")
load(".\\data\\Lasso\\AFW_gaps_real.RData")
load(".\\data\\Lasso\\AFW_gaps_by_iter_real.RData")
load(".\\data\\Lasso\\RAFW_gaps_by_iter_real.RData")
load(".\\data\\Lasso\\AFW_function_costs_real.RData")
load(".\\data\\Lasso\\RAFW_function_costs_real.RData")
load(".\\data\\Lasso\\AFW_support_size_real.RData")
load(".\\data\\Lasso\\RAFW_support_size_real.RData")
load(".\\data\\Lasso\\AFW_recovered_coeff_real.RData")
load(".\\data\\Lasso\\RAFW_recovered_coeff_real.RData")
load(".\\data\\Lasso\\PFW_gaps_real.RData")
load(".\\data\\Lasso\\PFW_gaps_by_iter_real.RData")
load(".\\data\\Lasso\\PFW_function_costs_real.RData")
load(".\\data\\Lasso\\PFW_support_size_real.RData")
load(".\\data\\Lasso\\PFW_recovered_coeff_real.RData")

# TODO: Move then into dataframe creation
df_afw_gaps$Method = "AFW"
df_afw_costs$Method = "AFW"
df_rafw_gaps$Method = "RAFW"
df_rafw_costs$Method = "RAFW"
df_afw_size_support$Method = "AFW"
df_rafw_size_support$Method = "RAFW"
df_afw_recovered_coeff$Method = "AFW"
df_rafw_recovered_coeff$Method = "RAFW"
df_afw_gaps_by_iterations$Method = "AFW"
df_rafw_gaps_by_iterations$Method = "RAFW"

# df_rafw_fw_steps$Method = "RAFW"
# df_rafw_drop_steps$Method = "RAFW"
# df_rafw_away_steps$Method = "RAFW"
# 
# df_afw_fw_steps$Method = "AFW"
# df_afw_drop_steps$Method = "AFW"
# df_afw_away_steps$Method = "AFW"

# df_rafw_fw_steps$Category = "FW Step"
# df_rafw_drop_steps$Method = "Drop Step"
# df_rafw_away_steps$Method = "Away Step"
# 
# df_afw_fw_steps$Method = "FW Step"
# df_afw_drop_steps$Method = "Drop Step"
# df_afw_away_steps$Method = "Away Step"

df_pfw_gaps$Method = "PFW"
df_pfw_costs$Method = "PFW"
df_pfw_size_support$Method = "PFW"
df_pfw_recovered_coeff$Method = "PFW"
df_pfw_gaps_by_iterations$Method = "PFW"

df_gaps = rbind(df_afw_gaps, df_rafw_gaps, df_pfw_gaps)
df_objective_functions = rbind(df_afw_costs, df_rafw_costs, df_pfw_costs)
df_size_support = rbind(df_afw_size_support, df_rafw_size_support, df_pfw_size_support)
df_recovered_coeff = rbind(df_afw_recovered_coeff, df_rafw_recovered_coeff, df_pfw_recovered_coeff)
df_gaps_by_iter = rbind(df_afw_gaps_by_iterations, df_rafw_gaps_by_iterations, df_pfw_gaps_by_iterations)
# df_fw_steps = rbind(df_afw_fw_steps, df_rafw_fw_steps)
# df_drop_steps = rbind(df_afw_drop_steps, df_rafw_drop_steps)
# df_away_steps = rbind(df_afw_away_steps, df_rafw_away_steps)

# df_steps = rbind(df_fw_steps, df_drop_steps, df_away_steps)
#### Plotting gaps by n_coef
setwd("C:\\Users\\ettag\\Documents\\GitHub\\Optimization-Project\\figures\\Lasso")
plot_gaps_by_coef = ggplot(df_gaps, aes(x = nbr_coef_per_grad, y = gap, color = Method)) +
  geom_line(size = 2) +
  scale_color_manual(values = c("#3C5688", "darkgreen", "darkorange")) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_minimal() + 
  labs(x = "No. of grad. coefficients", y = "Dual gap") + 
  theme(legend.position = "right",
        aspect.ratio = 1, 
        title = element_text(size = 28),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 28),
        axis.text.y = element_text(size = 28),
        axis.title.y = element_text(size = 28),
        axis.title.x = element_text(size = 28),
        legend.key = element_blank(), legend.key.size = unit(3,"line"),
        legend.title=element_text(size=28), legend.text=element_text(size=28))
ggsave(".\\plot_gaps_by_coef_real.pdf", plot_gaps_by_coef, device="pdf", width = 28, height = 28, units = "cm")

#### Plotting gap over iterations
plot_gaps_by_iter = ggplot(df_gaps_by_iter, aes(x = iteration, y = gap, color = Method)) +
  geom_line(size = 2) +
  scale_color_manual(values = c("#3C5688", "darkgreen", "darkorange")) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  labs(x = "Iteration", y = "Dual Gap") +
  theme_minimal() + 
  theme(legend.position = "right",
        aspect.ratio = 1, 
        title = element_text(size = 28),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 28),
        axis.text.y = element_text(size = 28),
        axis.title.y = element_text(size = 28),
        axis.title.x = element_text(size = 28),
        legend.key = element_blank(), legend.key.size = unit(3,"line"),
        legend.title=element_text(size=28), legend.text=element_text(size=28))
ggsave(".\\plot_gaps_by_iter_real.pdf", plot_gaps_by_iter, device="pdf", width = 28, height = 28, units = "cm")

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