library(ggplot2) # For plotting
library(tibble) # For dataframe creation
library(data.table) # For data table manipulation
library(scales) # For trans_breaks() in plots
library(ggsci)
rm(list = ls())

#### Loading data
setwd("C:\\Users\\ettag\\Documents\\GitHub\\Optimization-Project")
fw_times_python = read.csv(".\\data\\Lasso\\datasets\\fw_times_python.csv")
rfw_cs_times_python = read.csv(".\\data\\Lasso\\datasets\\rfw_times_cs_python.csv")
rfw_ncs_times_python = read.csv(".\\data\\Lasso\\datasets\\rfw_times_ncs_python.csv")

load(".\\data\\Lasso\\FW_iteration_times.RData")
load(".\\data\\Lasso\\RFW_iteration_times.RData")
#### Computing cumulative times

df_fw_times_python = tibble(iteration = fw_times_python$X, time = fw_times_python$X0, language = "Python", method = "FW")
df_fw_times_python$time = cumsum(df_fw_times_python$time)

#df_rfw_cs_times_python = tibble(iteration = rfw_cs_times_python$X, time = rfw_cs_times_python$X0, language = "Python", method = "RFW_cs")
#df_rfw_cs_times_python$time = cumsum(df_rfw_cs_times_python$time)

df_rfw_ncs_times_python = tibble(iteration = rfw_ncs_times_python$X, time = rfw_ncs_times_python$X0, language = "Python", method = "RFW")
df_rfw_ncs_times_python$time = cumsum(df_rfw_ncs_times_python$time)

df_fw_times_r = tibble(iteration = df_fw_iteration_times$iteration, time = df_fw_iteration_times$time, language = "R", method = "FW")
df_fw_times_r$time = cumsum(df_fw_times_r$time)

df_rfw_times_r = tibble(iteration = df_rfw_iteration_times$iteration, time = df_rfw_iteration_times$time, language = "R", method = "RFW")
df_rfw_times_r$time = cumsum(df_rfw_times_r$time)


df_times = rbind(df_fw_times_python, df_rfw_ncs_times_python, df_fw_times_r, df_rfw_times_r)

#### Plotting
setwd("C:\\Users\\ettag\\Documents\\GitHub\\Optimization-Project\\figures\\Lasso")
#### Plotting gaps over iterations
plot_iteration_time_python = ggplot(df_times, aes(x = iteration, y = time, color = language)) +
  geom_line(size = 2) +
  facet_grid(cols = vars(method)) + 
  scale_color_aaas(name = "Language") + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)), 
                limits = c(1, 10000)) +
  theme_minimal() + 
  labs(x = "Iteration", y = "Cumulative Time") + 
  theme(legend.position = "right",
        aspect.ratio = 1, 
        title = element_text(size = 24),
        strip.text.x = element_text(size = 28),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 24),
        axis.text.y = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        axis.title.x = element_text(size = 24),
        panel.spacing = unit(2, "lines"),
        legend.key = element_blank(), legend.key.size = unit(3,"line"),
        legend.title=element_text(size=28), legend.text=element_text(size=28))
plot_iteration_time_python
ggsave(".\\plot_iteration_time_languages.pdf", plot_iteration_time_python, device="pdf", width = 28, height = 15, units = "cm")
