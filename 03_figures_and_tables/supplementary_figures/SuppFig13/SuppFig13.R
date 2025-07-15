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

source("/home/gpalou/projects/NMD/scripts/NMD_efficiency/Figures/ggplot_themes.R")

################################
############ PLOT A ############
################################

input_figureA <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig13/panel_A.RData")

# Plot    
lambda <- unique(input_figureA$df_non_random$lambda)
plot_A_1 <- ggplot(input_figureA$df_non_random, aes(x = -log10(p_value_expected), y = -log10(p_value))) +
    geom_point(size = 2) + ylim(c(0,3.5)) +
    geom_abline(intercept = 0, slope = 1) +
    labs(y = expression("-log"[10] * "(P-values)"), x = expression("-log"[10] * "(expected P-values)"), 
        title = "Cancer and NMD genes") +
    annotate("text", x = 0.25, y = 3, label = paste0("Lambda \n= ",lambda), hjust = 0, vjust = 0) +
    ggplot_theme()
lambda <- unique(input_figureA$df_random$lambda)
plot_A_2 <- ggplot(input_figureA$df_random, aes(x = -log10(eval(parse(text="p_value_expected"))), y = -log10(eval(parse(text="p_value"))) )) +
    geom_point(size = 2) + ylim(c(0,3.5)) +
    geom_abline(intercept = 0, slope = 1) +
    labs(y = expression("-log"[10] * "(P-values)"), x = expression("-log"[10] * "(expected P-values)"), 
        title = "Random genes") +
    annotate("text", x = 0.25, y = 3, label = paste0("Lambda \n= ",lambda), hjust = 0, vjust = 0) +
    ggplot_theme()

####################################################
############ PANEL UP --> PLOTS A ##################
####################################################

plot_A <- ggarrange(plot_A_1,plot_A_2, labels = c("A",""),
                font.label = list(size = 20),
                ncol=2, nrow=1, common.legend = TRUE, legend = "bottom")
fig_up <- annotate_figure(plot_A, top = text_grob("ETG iNMDeff", 
                            size = 20, face = "plain")) 

################################
############ PLOT B ############
################################

input_figureB <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig13/panel_B.RData")

# Plot    
lambda <- unique(input_figureB$df_non_random$lambda)
plot_B_1 <- ggplot(input_figureB$df_non_random, aes(x = -log10(p_value_expected), y = -log10(p_value))) +
    geom_point(size = 2) + ylim(c(0,3.5)) +
    geom_abline(intercept = 0, slope = 1) +
    labs(y = expression("-log"[10] * "(P-values)"), x = expression("-log"[10] * "(expected P-values)"), 
        title = "Cancer and NMD genes") +
    annotate("text", x = 0.5, y = 3, label = paste0("Lambda = ",lambda), hjust = 0, vjust = 0) +
    ggplot_theme()
lambda <- unique(input_figureB$df_random$lambda)
plot_B_2 <- ggplot(input_figureB$df_random, aes(x = -log10(eval(parse(text="p_value_expected"))), y = -log10(eval(parse(text="p_value"))) )) +
    geom_point(size = 2) + ylim(c(0,3.5)) +
    geom_abline(intercept = 0, slope = 1) +
    labs(y = expression("-log"[10] * "(P-values)"), x = expression("-log"[10] * "(expected P-values)"), 
        title = "Random genes") +
    annotate("text", x = 0.5, y = 3, label = paste0("Lambda = ",lambda), hjust = 0, vjust = 0) +
    ggplot_theme()

########################################################
############ PANEL BOTTOM --> PLOTS B ##################
########################################################

plot_B <- ggarrange(plot_B_1,plot_B_2, labels = c("B",""),
                font.label = list(size = 20),
                ncol=2, nrow=1, common.legend = TRUE, legend = "bottom")
fig_bottom <- annotate_figure(plot_B, top = text_grob("ASE iNMDeff", 
                            size = 20, face = "plain")) 

# final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/SuppFig/test.png"
# ggsave(final_figure_path, plot_B, width = 400, height = 350, units = "mm")

########################################
############ FINAL SUPP FIG ############
########################################

SuppFig_final <- plot_grid(fig_up, fig_bottom, nrow = 2, ncol = 1, rel_heights = c(0.5,0.5))

final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/SuppFig/SuppFig13_complete.png"
ggsave(final_figure_path, SuppFig_final, width = 150, height = 175, units = "mm")
