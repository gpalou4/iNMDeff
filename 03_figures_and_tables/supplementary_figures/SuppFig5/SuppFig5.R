library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(viridis)
library(grid)
library(gridExtra)
library(corrplot)
library(ggcorrplot)

source("/home/gpalou/projects/NMD/scripts/NMD_efficiency/Figures/ggplot_themes.R")

################################
############ PLOT A ############
################################

input_figureA <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig5/panel_A.RData")
input_figureA$dataset <- ifelse(input_figureA$dataset == "GTEx","GTex","TCGA")
input_figureA$type <- factor(input_figureA$type, levels = c("threshold","subset"))
# Plot
plot_A <- ggplot(input_figureA, aes(x = factor(threshold), y = corr, fill = factor(type), 
        color = factor(highlight), linetype = factor(highlight)) ) +
  geom_bar(stat = "identity", position = position_dodge(), size = 0.5) +
  geom_text(aes(label = sample_size), position = position_dodge(width = 0.9), vjust = 0.5, hjust = 1.5, size = 2, angle = 90) +
  # geom_text(aes(label = sample_size_perc), position = position_dodge(width = 0.9), vjust = -0.5, size = 2) +
  geom_text(aes(label = significance), position = position_dodge(width = 0.9), vjust = -2, size = 2) +
  labs(x = "Number of PTCs per sample (ASE iNMDeff method)", 
        y = "ASE vs ETG iNMDeff correlation", 
        fill = "Samples used") +
#   facet_wrap(dataset ~ VAF, scale = "free") +
  facet_grid(VAF ~ dataset, scales = "free_y") +
  scale_fill_manual(
    name = "Samples used",
    labels = c("thresholded subset", "random, equal-sized subset"),
    values = c(
      "threshold" = "#3cdaa0",
      "subset" = "#91afe2"
    )) +
  scale_color_manual(
    values = c("normal" = "black", "highlight" = "red")
  ) +
  scale_linetype_manual(
    values = c("normal" = "solid", "highlight" = "dashed")
  ) +
  guides(color = "none", linetype = "none") +  # Remove the fill legend
  ggplot_theme() +
  theme(
      axis.text.x = element_text(),
      legend.position = "top",
      strip.background = element_blank(),  # Remove panel background
      strip.text = element_text(size = 12) # Keep titles
        )

################################
############ PLOT B ############
################################

input_figureB <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig5/panel_B.RData")

input_figureB$VAF <- ifelse(input_figureB$VAF == "0.01","AF <= 0.01","AF <= 0.2")
input_figureB$VAF <- factor(input_figureB$VAF, levels = c("AF <= 0.2","AF <= 0.01"))
 
plot_B <- ggplot(input_figureB, aes(x = threshold)) +
  # geom_line(aes(y = mean_correlation, color = "Mean correlation"), size = 1) +
  geom_errorbar(
    aes(
      ymin = Q1, # Lower bound of whiskers (Q1)
      ymax = Q3, # Upper bound of whiskers (Q3),
      color = "Q1-Q3 range"
    ),
    width = 0.2, # Width of whiskers
    size = 0.8 # Thickness of the whisker lines
  ) +
  geom_line(aes(y = median_correlation, color = "Median correlation"), size = 1, linetype = "dashed") +
  # geom_line(aes(y = positive_correlation_percentage / 100, color = "% of tissues with positive correlation"), size = 1, linetype = "dotted") +
  geom_bar(aes(y = positive_correlation_percentage / 100, fill = "% of tissues with positive correlation"), 
           stat = "identity", position = "dodge", alpha = 0.6) +
  facet_grid(VAF ~ dataset, scales = "free_y") +
#   facet_wrap(dataset ~ VAF, scale = "free") +
  scale_x_continuous(
    breaks = 1:10  # Explicitly set X-axis breaks from 1 to 10
  ) +
  scale_y_continuous(
    breaks = c(0,0.25,0.5,0.75,1)  # Explicitly set X-axis breaks from 1 to 10
  ) +  
  scale_color_manual(
    values = c("Mean correlation" = "blue", 
               "Median correlation" = "green", 
               "% of tissues with positive correlation" = "red")
  ) +
  coord_cartesian(ylim = c(0, 1)) + # Zoom into the y-axis range between 0 and 1
  labs(
    title = "",
    fill = "",
    x = "Number of PTCs per sample (ASE iNMDeff method)",
    y = "ASE vs ETG iNMDeff correlation",
    color = ""
  ) +
  ggplot_theme() +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5),
    strip.text = element_text(size = 12)
  )

##################################################
############ PANEL UP --> PLOTS A ################
##################################################

fig_up <- plot_grid(plot_A, nrow = 1, ncol = 1, labels = c("A"), 
                    label_size = 20, vjust = 1, rel_widths = c(0.5,0.5))

######################################################
############ PANEL BOTTOM --> PLOTS B ################
######################################################

fig_bottom <- plot_grid(plot_B, nrow = 1, ncol = 1, labels = c("B"), label_size = 20, vjust = 1)

##################################################
############ FINAL FIG --> SUPP FIG 2 ############
##################################################

SuppFig_final <- plot_grid(fig_up, fig_bottom, nrow = 2, rel_heights = c(0.55,0.45))

final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/SuppFig/SuppFig5_complete.png"
ggsave(final_figure_path, SuppFig_final, width = 160, height = 225, units = "mm")
