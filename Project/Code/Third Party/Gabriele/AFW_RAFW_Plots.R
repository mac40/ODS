library(ggplot2) # For plotting
library(tibble) # For dataframe creation
library(data.table) # For data table manipulation
library(scales) # For trans_breaks() in plots


rm(list = ls())

#### Loading data
setwd("C:\\Users\\ettag\\Documents\\GitHub\\Optimization-Project")
load(".\\data\\Lasso\\RAFW_gaps.RData")
load(".\\data\\Lasso\\AFW_gaps.RData")
load(".\\data\\Lasso\\AFW_gaps_by_iter.RData")
load(".\\data\\Lasso\\RAFW_gaps_by_iter.RData")
load(".\\data\\Lasso\\AFW_function_costs.RData")
load(".\\data\\Lasso\\RAFW_function_costs.RData")
load(".\\data\\Lasso\\AFW_support_size.RData")
load(".\\data\\Lasso\\RAFW_support_size.RData")
load(".\\data\\Lasso\\AFW_recovered_coeff.RData")
load(".\\data\\Lasso\\RAFW_recovered_coeff.RData")
load(".\\data\\Lasso\\PFW_gaps.RData")
load(".\\data\\Lasso\\PFW_gaps_by_iter.RData")
load(".\\data\\Lasso\\PFW_function_costs.RData")
load(".\\data\\Lasso\\PFW_support_size.RData")
load(".\\data\\Lasso\\PFW_recovered_coeff.RData")
load(".\\data\\Lasso\\RAFW_FW_steps.RData")
load(".\\data\\Lasso\\RAFW_drop_steps.RData")
load(".\\data\\Lasso\\RAFW_away_steps.RData")
load(".\\data\\Lasso\\AFW_FW_steps.RData")
load(".\\data\\Lasso\\AFW_drop_steps.RData")
load(".\\data\\Lasso\\AFW_away_steps.RData")

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
df_pfw_gaps$Method = "PFW"
df_pfw_costs$Method = "PFW"
df_pfw_size_support$Method = "PFW"
df_pfw_recovered_coeff$Method = "PFW"
df_pfw_gaps_by_iterations$Method = "PFW"

df_rafw_fw_steps$Method = "RAFW"
df_rafw_drop_steps$Method = "RAFW"
df_rafw_away_steps$Method = "RAFW"

df_afw_fw_steps$Method = "AFW"
df_afw_drop_steps$Method = "AFW"
df_afw_away_steps$Method = "AFW"

df_rafw_fw_steps$Category = "FW Step"
df_rafw_drop_steps$Category = "Drop Step"
df_rafw_away_steps$Category = "Away Step"

df_afw_fw_steps$Category = "FW Step"
df_afw_drop_steps$Category = "Drop Step"
df_afw_away_steps$Category = "Away Step"

df_gaps = rbind(df_afw_gaps, df_rafw_gaps, df_pfw_gaps)
df_objective_functions = rbind(df_afw_costs, df_rafw_costs, df_pfw_costs)
df_size_support = rbind(df_afw_size_support, df_rafw_size_support, df_pfw_size_support)
df_recovered_coeff = rbind(df_afw_recovered_coeff, df_rafw_recovered_coeff, df_pfw_recovered_coeff)
df_gaps_by_iter = rbind(df_afw_gaps_by_iterations, df_rafw_gaps_by_iterations, df_pfw_gaps_by_iterations)

df_fw_steps = rbind(df_afw_fw_steps, df_rafw_fw_steps)
df_drop_steps = rbind(df_afw_drop_steps, df_rafw_drop_steps)
df_away_steps = rbind(df_afw_away_steps, df_rafw_away_steps)

df_steps = rbind(df_fw_steps, df_drop_steps, df_away_steps)

setwd("C:\\Users\\ettag\\Documents\\GitHub\\Optimization-Project\\figures\\Lasso")
#### Plotting gaps by n_coef
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
ggsave(".\\plot_gaps_by_coef.pdf", plot_gaps_by_coef, device="pdf", width = 28, height = 28, units = "cm")

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
ggsave(".\\plot_gaps_by_iter.pdf", plot_gaps_by_iter, device="pdf", width = 28, height = 28, units = "cm")

#### Plotting function costs over iterations
plot_of_by_niterations = ggplot(df_objective_functions, aes(x = iteration, y = function_cost, color = Method)) +
  geom_line(size = 2) +
  scale_color_manual(values = c("#3C5688", "darkgreen", "darkorange")) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
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
  


#### Plotting the size of the support over iterations
plot_size_support_by_iterations = ggplot(df_size_support, aes(x = iteration, y = size_support, color = Method)) +
  geom_line(size = 2) +
  scale_color_manual(values = c("#3C5688", "darkgreen", "darkorange")) + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+ 
  labs(x = "Iteration", y = "Size of the support") +
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
ggsave(".\\plot_size_support_by_iterations.pdf", plot_size_support_by_iterations, device="pdf", width = 28, height = 28, units = "cm")

# Plotting recovered coeff
recovered_coeff_by_iterations = ggplot(df_recovered_coeff, aes(x = iteration, y = recovered_coeff, color = Method)) +
  geom_line(size = 2) +
  scale_color_manual(values = c("#3C5688", "darkgreen", "darkorange")) + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + 
  labs(x = "Iteration", y = "Recovered coefficients") +
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
ggsave(".\\recovered_coeff_by_iterations.pdf", recovered_coeff_by_iterations, device="pdf", width = 28, height = 28, units = "cm")

steps_categories = ggplot(df_steps) +
  geom_bar(stat="identity", position = "dodge", aes(x = Category, y = step_value, fill = Method)) +
  scale_fill_manual(values = c("#3C5688", "darkorange")) +
  labs(x = "Step Type", y = "Number of Steps") +
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
steps_categories
ggsave(".\\steps_categories.pdf", steps_categories, device="pdf", width = 28, height = 28, units = "cm")

library(patchwork)

(plot_gaps_by_iter + plot_size_support_by_iterations) / (plot_gaps_by_coef + recovered_coeff_by_iterations)
