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

input_figureA <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig15/panel_A.RData")
# Plot
sig_CNV_PCs <- c(3,52,86)
plot_A <- input_figureA %>%
    filter(CNV_PCs %in% sig_CNV_PCs) %>%
    ggplot(aes(y = cancer, x = coefficient, color = factor(CNV_PCs))) +
        geom_point(size = 3, shape = 16) +
        facet_wrap( ~ CNV_PCs + NMD_method, scales = "free_x", nrow = 1) +
        geom_errorbarh(aes(xmin = CI_2.5_values, xmax = CI_97.5_values), size = 1, height = 0.1) +
        geom_vline(xintercept = 0, color = "red", linetype = "dashed", cex = 0.8, alpha = 0.5) +
        labs(x = "Association effect size", y = "", color = "PCs",
                title = paste0("Significant pan-cancer CNA-PCs at FDR < 10%")) + 
        scale_color_brewer(palette = "Paired") +
        ggplot_theme() + 
        theme(legend.position = "top",
                panel.spacing = grid::unit(1, "lines"),
                panel.background = element_rect(colour = "black"),
                strip.background = element_rect(colour = "black"),
                strip.text = element_text(colour = "black", size = 12)) +
        scale_x_continuous(labels = label_number(scale = 1, accuracy = 0.1)) #+

####################################################
############ PANEL UP --> PLOTS A ##################
####################################################

fig_up <- plot_grid(plot_A, nrow = 1, ncol = 1, labels = c(""))

########################################
############ FINAL SUPP FIG ############
########################################

SuppFig_final <- plot_grid(fig_up, nrow = 1, ncol = 1, rel_heights = c(0.5))

final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/SuppFig/SuppFig15_complete.png"
ggsave(final_figure_path, SuppFig_final, width = 275, height = 150, units = "mm")

