library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(viridis)
library(grid)
library(gridExtra)
library(corrplot)
library(ggcorrplot)
library(dplyr)
library(ggpubr)
library(ggrepel)
library(scales)
library(ggh4x)

source("/home/gpalou/projects/NMD/scripts/NMD_efficiency/Figures/ggplot_themes.R")

################################
############ PLOT A ############
################################

input_figureA <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig24/panel_A.RData")
colnames(input_figureA)[11] <- "NMD_method_FDR"
# Plot

# Combine significant and moderate points into one dataset with a category column
labeled_points <- input_figureA %>%
  mutate(significance = case_when(
    SKATO_p_value_FDR_adjusted < 0.05 ~ "Significant",
    SKATO_p_value_FDR_adjusted >= 0.05 & SKATO_p_value_FDR_adjusted < 0.2 ~ "Moderate",
    TRUE ~ NA_character_)) %>%
  filter(!is.na(significance))

# Plot
plot_A <- input_figureA %>%
  ggplot(aes(y = gene, x = coefficient)) +
  geom_boxplot() +
  geom_point(size = 2, alpha = 0.3) +
  facet_nested(~ NMD_method, scales = "free_x") +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", cex = 0.8, alpha = 0.5) +
  labs(x = "Association effect size", y = "", fill = "", title = "TCGA") + 
  scale_color_brewer(palette = "Paired") +
  ggplot_theme_bw() + 
  xlim(-3.5, 3.5) +
  theme(legend.position = "none") +
  # Add combined labels with arrows pointing to the points
  geom_text_repel(data = labeled_points, aes(x = coefficient, label = tissue, color = significance), 
                   nudge_y = 0.1, size = 2, max.overlaps = Inf, 
                   arrow = arrow(length = unit(0.03, "npc"), type = "closed"), # Add arrows
                   label.size = 0.2) +  # Optional: Adjust label box size
  
  # Set colors for significant and moderate points
  scale_color_manual(values = c("Significant" = "#6b6bb6", "Moderate" = "#eb3131"))

####################################################
############ PANEL UP --> PLOTS A ##################
####################################################

fig_up <- plot_grid(plot_A, nrow = 1, ncol = 1, labels = c("A"), rel_widths = c(0.5))

################################
############ PLOT B ############
################################

input_figureB <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig24/panel_B.RData")
colnames(input_figureB)[11] <- "NMD_method_FDR"

# Combine significant and moderate points into one dataset with a category column
labeled_points <- input_figureB %>%
  mutate(significance = case_when(
    SKATO_p_value_FDR_adjusted < 0.05 ~ "Significant",
    SKATO_p_value_FDR_adjusted >= 0.05 & SKATO_p_value_FDR_adjusted < 0.2 ~ "Moderate",
    TRUE ~ NA_character_)) %>%
  filter(!is.na(significance))

# Plot
plot_B <- input_figureB %>%
  ggplot(aes(y = gene, x = coefficient)) +
  geom_boxplot() +
  geom_point(size = 2, alpha = 0.3) +
  facet_nested(~ NMD_method, scales = "free_x") +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", cex = 0.8, alpha = 0.5) +
  labs(x = "Association effect size", y = "", fill = "", title = "GTex") + 
  scale_color_brewer(palette = "Paired") +
  ggplot_theme_bw() + 
  xlim(-2.5, 2.5) +
  theme(legend.position = "none") +
  # Add combined labels with arrows pointing to the points
  geom_text_repel(data = labeled_points, aes(x = coefficient, label = tissue, color = significance), 
                   nudge_y = 0.1, size = 2, max.overlaps = Inf, 
                   arrow = arrow(length = unit(0.03, "npc"), type = "closed"), # Add arrows
                   label.size = 0.2) +  # Optional: Adjust label box size
  
  # Set colors for significant and moderate points
  scale_color_manual(values = c("Significant" = "#6b6bb6", "Moderate" = "#eb3131"))

########################################################
############ PANEL BOTTOM --> PLOTS B ##################
########################################################

fig_bottom <- plot_grid(plot_B, nrow = 1, ncol = 1, labels = c("B"), rel_widths = c(0.5))

########################################
############ FINAL SUPP FIG ############
########################################

SuppFig_final <- plot_grid(fig_up, fig_bottom, nrow = 2, ncol = 1, rel_heights = c(0.5,0.5))

final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/SuppFig/SuppFig24_complete.png"
ggsave(final_figure_path, SuppFig_final, width = 175, height = 175, units = "mm")
final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/SuppFig/SuppFig24_complete.pdf"
ggsave(final_figure_path, SuppFig_final, width = 175, height = 175, units = "mm")
