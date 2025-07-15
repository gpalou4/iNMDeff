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

source("/home/gpalou/projects/NMD/scripts/NMD_efficiency/Figures/ggplot_themes.R")

################################
############ PLOT A ############
################################

input_figureA <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig18/panel_A.RData")

# Plot
plot_A <- ggplot(input_figureA, aes(y = round(samples_percentage,2), x = cancer_type, fill = CNA_PC_bins)) + 
    geom_bar(position="stack", stat="identity") +
    labs(title = paste0("CNA-PC3"), x = "", y = "% of individuals", fill = paste0("Group")) + 
    scale_fill_brewer(labels = c("Low", "Mid", "High"), palette = "Set2", direction = -1) +
    ggplot_theme() +
    theme(axis.text.x = element_text(angle = 90, hjust=0.5, vjust=0.5, size = 9),
            legend.position = "top")

################################
############ PLOT B ############
################################

input_figureB <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig18/panel_B.RData")

# Plot
plot_B <- ggplot(input_figureB, aes(y = round(samples_percentage,2), x = cancer_type, fill = CNA_PC_bins)) + 
    geom_bar(position="stack", stat="identity") +
    labs(title = paste0("CNA-PC86"), x = "", y = "% of individuals", fill = paste0("Group")) + 
    scale_fill_brewer(labels = c("Low", "Mid", "High"), palette = "Set2", direction = -1) +
    ggplot_theme() +
    theme(axis.text.x = element_text(angle = 90, hjust=0.5, vjust=0.5, size = 9),
            legend.position = "top")

######################################################
############ PANEL UP --> PLOTS A-B ##################
######################################################

fig_up <- plot_grid(plot_A, plot_B, nrow = 1, ncol = 2, labels = c("A","B"))

########################################################
############ PANEL BOTTOM --> PLOTS C ##################
########################################################

fig_bottom <- plot_grid(plot_C, ggplot() + theme_classic(), nrow = 1, ncol = 2, labels = c("C",""))

########################################
############ FINAL SUPP FIG ############
########################################

SuppFig_final <- plot_grid(fig_up, nrow = 1, ncol = 1, rel_heights = c(0.5))

final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/SuppFig/SuppFig18_complete.png"
ggsave(final_figure_path, SuppFig_final, width = 200, height = 100, units = "mm")
final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/SuppFig/SuppFig18_complete.pdf"
ggsave(final_figure_path, SuppFig_final, width = 200, height = 100, units = "mm")
