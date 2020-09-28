library(ggplot2) # For plotting
library(tibble) # For dataframe creation
library(data.table) # For data table manipulation
library(scales) # For trans_breaks() in plots

rm(list = ls())

#### Loading data
setwd("C:\\Users\\ettag\\Documents\\GitHub\\Optimization-Project")
load(".\\data\\LGL\\FW_gaps_by_coef.RData")
load(".\\data\\LGL\\FW_gaps_by_time.RData")
load(".\\data\\LGL\\RFW_gaps_by_time.RData")
load(".\\data\\LGL\\RFW_gaps_by_coef.RData")


df_fw_gaps_by_time$Method = "FW"
df_fw_gaps_by_coef$Method = "FW"
df_rfw_gaps_by_coef$Method = "RFW"
df_rfw_gaps_by_time$Method = "RFW"

df_gaps_by_time = rbind(df_fw_gaps_by_time, df_rfw_gaps_by_time)
df_gaps_by_coef = rbind(df_fw_gaps_by_coef, df_rfw_gaps_by_coef)

setwd("C:\\Users\\ettag\\Documents\\GitHub\\Optimization-Project\\figures\\LGL")

#### Plotting gaps over iterations
plot_gaps_by_coef = ggplot(df_gaps_by_coef, aes(x = nbr_coef_per_grad, y = gap, color = Method)) +
  geom_line(size = 2) +
  scale_color_manual(values = c("#3C5688", "darkorange")) + 
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
plot_gaps_by_coef
ggsave(".\\plot_gaps_by_coef.pdf", plot_gaps_by_coef, device="pdf", width = 28, height = 28, units = "cm")

#### Plotting gap over iterations
plot_gaps_by_time = ggplot(df_gaps_by_time, aes(x = time, y = gap, color = Method)) +
  geom_line(size = 2) +
  scale_color_manual(values = c("#3C5688", "darkorange")) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_minimal() + 
  labs(x = "Time (s)", y = "Dual gap") + 
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

plot_gaps_by_time
ggsave(".\\plot_gaps_by_time.pdf", plot_gaps_by_time, device="pdf", width = 28, height = 28, units = "cm")

library(patchwork)

(plot_gaps_by_coef + plot_gaps_by_time)
