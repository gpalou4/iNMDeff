library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(viridis)
library(grid)
library(gridExtra)
library(stringr)
library(dplyr)
library(png)
library(ggpubr)
library(ggbreak)
library(scales)
library(tibble)
library(corrplot)
library(ggcorrplot)
library(MASS)
library(ggh4x)
library(survminer)

source("/home/gpalou/projects/NMD/scripts/NMD_efficiency/Figures/ggplot_themes.R")

# Function to save ggsurv plots with "ggsave"
ggsave_workaround <- function(g){survminer:::.build_ggsurvplot(x = g,
                                                               surv.plot.height = NULL,
                                                               risk.table.height = NULL,
                                                               ncensor.plot.height = NULL)} 

################################
############ PLOT A ############
################################

input_figureA <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig27/panel_A_H.RData")

surv_obj <- input_figureA[["SKCM_all_OS_ETG"]]$surv_obj
km_fit <- input_figureA[["SKCM_all_OS_ETG"]]$km_fit
df <- input_figureA[["SKCM_all_OS_ETG"]]$df

plot_A <- ggsurvplot(km_fit, 
        data = df, risk.table = FALSE,
        pval = TRUE, conf.int = FALSE,
        ylim = c(0, 1), #xlim = c(0, 22), 
        xlab = "Time (years)", ylab = "Survival probability",
        title = "TCGA - SKCM\nP20",
        ggtheme = ggplot_theme_bw(), surv.median.line = "v", 
        legend.labs = c("High", "Low"), 
        legend.title = paste0("ETG iNMDeff"),
        legend = "top") 

plot_A_to_save <- ggsave_workaround(plot_A)

################################
############ PLOT B ############
################################

input_figureB <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig27/panel_A_H.RData")

surv_obj <- input_figureB[["PRAD_all_OS_ETG"]]$surv_obj
km_fit <- input_figureB[["PRAD_all_OS_ETG"]]$km_fit
df <- input_figureB[["PRAD_all_OS_ETG"]]$df

plot_B <- ggsurvplot(km_fit, 
        data = df, risk.table = FALSE,
        pval = TRUE, conf.int = FALSE,
        ylim = c(0, 1), #xlim = c(0, 22), 
        xlab = "Time (years)", ylab = "Survival probability",
        title = "TCGA - PRAD\nP35",
        ggtheme = ggplot_theme_bw(), surv.median.line = "v", 
        legend.labs = c("High", "Low"), 
        legend.title = paste0("ETG iNMDeff"),
        legend = "top") 

plot_B_to_save <- ggsave_workaround(plot_B)

################################
############ PLOT C ############
################################

input_figureC <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig27/panel_A_H.RData")

surv_obj <- input_figureC[["ESCA_ac_all_OS_ETG"]]$surv_obj
km_fit <- input_figureC[["ESCA_ac_all_OS_ETG"]]$km_fit
df <- input_figureC[["ESCA_ac_all_OS_ETG"]]$df

plot_C <- ggsurvplot(km_fit, 
        data = df, risk.table = FALSE,
        pval = TRUE, conf.int = FALSE,
        ylim = c(0, 1), #xlim = c(0, 22), 
        xlab = "Time (years)", ylab = "Survival probability",
        title = "TCGA - ESCA_ac\nP30",
        ggtheme = ggplot_theme_bw(), surv.median.line = "v", 
        legend.labs = c("High", "Low"), 
        legend.title = paste0("ETG iNMDeff"),
        legend = "top") 

plot_C_to_save <- ggsave_workaround(plot_C)

################################
############ PLOT D ############
################################

input_figureD <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig27/panel_A_H.RData")

surv_obj <- input_figureC[["LUAD_all_OS_ETG"]]$surv_obj
km_fit <- input_figureC[["LUAD_all_OS_ETG"]]$km_fit
df <- input_figureC[["LUAD_all_OS_ETG"]]$df

plot_D <- ggsurvplot(km_fit, 
        data = df, risk.table = FALSE,
        pval = TRUE, conf.int = FALSE,
        ylim = c(0, 1), #xlim = c(0, 22), 
        xlab = "Time (years)", ylab = "Survival probability",
        title = "TCGA - LUAD\nP10",
        ggtheme = ggplot_theme_bw(), surv.median.line = "v", 
        legend.labs = c("High", "Low"), 
        legend.title = paste0("ETG iNMDeff"),
        legend = "top") 

plot_D_to_save <- ggsave_workaround(plot_D)

################################
############ PLOT E ############
################################

input_figureE <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig27/panel_A_H.RData")

surv_obj <- input_figureA[["SKCM_all_OS_ASE"]]$surv_obj
km_fit <- input_figureA[["SKCM_all_OS_ASE"]]$km_fit
df <- input_figureA[["SKCM_all_OS_ASE"]]$df

plot_E <- ggsurvplot(km_fit, 
        data = df, risk.table = FALSE,
        pval = TRUE, conf.int = FALSE,
        ylim = c(0, 1), #xlim = c(0, 22), 
        xlab = "Time (years)", ylab = "Survival probability",
        title = "TCGA - SKCM\nP20",
        ggtheme = ggplot_theme_bw(), surv.median.line = "v", 
        legend.labs = c("High", "Low"), 
        legend.title = paste0("ASE iNMDeff"),
        legend = "top") 

plot_E_to_save <- ggsave_workaround(plot_E)

################################
############ PLOT F ############
################################

input_figureF <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig27/panel_A_H.RData")

surv_obj <- input_figureF[["PRAD_all_OS_ASE"]]$surv_obj
km_fit <- input_figureF[["PRAD_all_OS_ASE"]]$km_fit
df <- input_figureF[["PRAD_all_OS_ASE"]]$df

plot_F <- ggsurvplot(km_fit, 
        data = df, risk.table = FALSE,
        pval = TRUE, conf.int = FALSE,
        ylim = c(0, 1), #xlim = c(0, 22), 
        xlab = "Time (years)", ylab = "Survival probability",
        title = "TCGA - PRAD\nP35",
        ggtheme = ggplot_theme_bw(), surv.median.line = "v", 
        legend.labs = c("High", "Low"), 
        legend.title = paste0("ASE iNMDeff"),
        legend = "top") 

plot_F_to_save <- ggsave_workaround(plot_F)

################################
############ PLOT G ############
################################

input_figureG <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig27/panel_A_H.RData")

surv_obj <- input_figureG[["ESCA_ac_all_OS_ASE"]]$surv_obj
km_fit <- input_figureG[["ESCA_ac_all_OS_ASE"]]$km_fit
df <- input_figureG[["ESCA_ac_all_OS_ASE"]]$df

plot_G <- ggsurvplot(km_fit, 
        data = df, risk.table = FALSE,
        pval = TRUE, conf.int = FALSE,
        ylim = c(0, 1), #xlim = c(0, 22), 
        xlab = "Time (years)", ylab = "Survival probability",
        title = "TCGA - ESCA_ac\nP30",
        ggtheme = ggplot_theme_bw(), surv.median.line = "v", 
        legend.labs = c("High", "Low"), 
        legend.title = paste0("ASE iNMDeff"),
        legend = "top") 

plot_G_to_save <- ggsave_workaround(plot_G)

################################
############ PLOT H ############
################################

input_figureH <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig27/panel_A_H.RData")

surv_obj <- input_figureH[["LUAD_all_OS_ASE"]]$surv_obj
km_fit <- input_figureH[["LUAD_all_OS_ASE"]]$km_fit
df <- input_figureH[["LUAD_all_OS_ASE"]]$df

plot_H <- ggsurvplot(km_fit, 
        data = df, risk.table = FALSE,
        pval = TRUE, conf.int = FALSE,
        ylim = c(0, 1), #xlim = c(0, 22), 
        xlab = "Time (years)", ylab = "Survival probability",
        title = "TCGA - LUAD\nP10",
        ggtheme = ggplot_theme_bw(), surv.median.line = "v", 
        legend.labs = c("High", "Low"), 
        legend.title = paste0("ASE iNMDeff"),
        legend = "top") #+ guides(color = guide_legend(ncol = 1))

plot_H_to_save <- ggsave_workaround(plot_H)

###############################################
############ PANEL UP --> PLOT A-D ############
###############################################

fig_up <- plot_grid(plot_B_to_save, plot_F_to_save, plot_A_to_save, plot_E_to_save, nrow = 1, ncol = 4,  
                        labels = c("A","B","C","D"), label_size = 20)

###################################################
############ PANEL BOTTOM --> PLOT E-F ############
###################################################

fig_bottom <- plot_grid(plot_C_to_save, plot_G_to_save, plot_D_to_save, plot_H_to_save, nrow = 1, ncol = 4,  
                        labels = c("E","F","G","H"), label_size = 20)

# fig_bottom <- ggarrange(plot_C_to_save,plot_G_to_save,plot_D_to_save, plot_H_to_save, labels = c("E","F","G","H"),
#                 font.label = list(size = 20), nrow = 1, ncol = 4, common.legend = TRUE, legend = "bottom")

########################################
############ FINAL SUPP FIG ############
########################################

SuppFig_final <- plot_grid(fig_up, fig_bottom, nrow = 2, ncol = 1)

final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/SuppFig/SuppFig27_complete.png"
ggsave(final_figure_path, SuppFig_final, width = 265, height = 160, units = "mm")
final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/SuppFig/SuppFig27_complete.pdf"
ggsave(final_figure_path, SuppFig_final, width = 265, height = 160, units = "mm")

