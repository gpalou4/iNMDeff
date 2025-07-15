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

input_figureA <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig29/panel_A_C.RData")

# Plot A1
surv_obj <- input_figureA[["SKCM_no_treatment_PFS_ASE"]]$surv_obj
km_fit <- input_figureA[["SKCM_no_treatment_PFS_ASE"]]$km_fit
df <- input_figureA[["SKCM_no_treatment_PFS_ASE"]]$df

plot_A <- ggsurvplot(km_fit, 
        data = df, risk.table = FALSE,
        pval = TRUE, conf.int = FALSE,
        ylim = c(0, 1), #xlim = c(0, 22), 
        xlab = "Time (years)", ylab = "PFS probability",
        title = "TCGA - SKCM \nw/o available treament (P50)",
        ggtheme = ggplot_theme_bw(), surv.median.line = "v", 
        legend.labs = c("High", "Low"), 
        legend.title = paste0("ASE iNMDeff"),
        legend = "top") 

plot_A_to_save <- ggsave_workaround(plot_A)

################################
############ PLOT B ############
################################

input_figureB <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig29/panel_A_C.RData")

surv_obj <- input_figureB[["SKCM_immunotherapy_PFS_ASE"]]$surv_obj
km_fit <- input_figureB[["SKCM_immunotherapy_PFS_ASE"]]$km_fit
df <- input_figureB[["SKCM_immunotherapy_PFS_ASE"]]$df

plot_B <- ggsurvplot(km_fit, 
        data = df, risk.table = FALSE,
        pval = TRUE, conf.int = FALSE,
        ylim = c(0, 1), #xlim = c(0, 22), 
        xlab = "Time (years)", ylab = "PFS probability",
        title = "TCGA - SKCM with\n at least immunotherapy (P50)",
        ggtheme = ggplot_theme_bw(), surv.median.line = "v", 
        legend.labs = c("High", "Low"), 
        legend.title = paste0("ASE iNMDeff"),
        legend = "top") 

plot_B_to_save <- ggsave_workaround(plot_B)

################################
############ PLOT C ############
################################

input_figureC <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig29/panel_A_C.RData")

surv_obj <- input_figureC[["SKCM_chemotherapy_PFS_ASE"]]$surv_obj
km_fit <- input_figureC[["SKCM_chemotherapy_PFS_ASE"]]$km_fit
df <- input_figureC[["SKCM_chemotherapy_PFS_ASE"]]$df

plot_C <- ggsurvplot(km_fit, 
        data = df, risk.table = FALSE,
        pval = TRUE, conf.int = FALSE,
        ylim = c(0, 1), #xlim = c(0, 22), 
        xlab = "Time (years)", ylab = "PFS probability",
        title = "TCGA - SKCM with chemotherapy, \nexcluding immunotherapy (P50)",
        ggtheme = ggplot_theme_bw(), surv.median.line = "v", 
        legend.labs = c("High", "Low"), 
        legend.title = paste0("ASE iNMDeff"),
        legend = "top") 

plot_C_to_save <- ggsave_workaround(plot_C)

################################################
############ PANEL UP --> PLOTS A-C ############
################################################

fig_up <- plot_grid(plot_A_to_save, plot_B_to_save, plot_C_to_save, nrow = 1, ncol = 3,  
                        labels = c("A","B","C"), label_size = 20)

################################
############ PLOT D ############
################################

input_figureD <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig29/panel_D_F.RData")

surv_obj <- input_figureD[["Liu_SKCM_PFS"]][["ASE_PTC_NMD_triggering_0.2"]]$surv_obj
km_fit <- input_figureD[["Liu_SKCM_PFS"]][["ASE_PTC_NMD_triggering_0.2"]]$km_fit
df <- input_figureD[["Liu_SKCM_PFS"]][["ASE_PTC_NMD_triggering_0.2"]]$df
# Correlation
cor.test(df$NMDeff,df$PFS, method = "pearson")

plot_D <- ggsurvplot(km_fit, 
        data = df, risk.table = FALSE,
        pval = TRUE, conf.int = FALSE,
        ylim = c(0, 1), #xlim = c(0, 22), 
        xlab = "Time (months)", ylab = "PFS probability",
        title = "Liu 2019 - SKCM\nwith immunotherapy (P40)",
        ggtheme = ggplot_theme_bw(), surv.median.line = "v", 
        legend.labs = c("High", "Low"), 
        legend.title = paste0("ASE iNMDeff"),
        legend = "top") 

plot_D_to_save <- ggsave_workaround(plot_D)

################################
############ PLOT E ############
################################

input_figureE <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig29/panel_D_F.RData")

surv_obj <- input_figureE[["Carrol_EAC_PFS"]][["ASE_PTC_NMD_triggering_0.2"]]$surv_obj
km_fit <- input_figureE[["Carrol_EAC_PFS"]][["ASE_PTC_NMD_triggering_0.2"]]$km_fit
df <- input_figureE[["Carrol_EAC_PFS"]][["ASE_PTC_NMD_triggering_0.2"]]$df
# Correlation
cor.test(df$NMDeff,df$PFS, method = "pearson")

plot_E <- ggsurvplot(km_fit, 
        data = df, risk.table = FALSE,
        pval = TRUE, conf.int = FALSE,
        ylim = c(0, 1), #xlim = c(0, 22), 
        xlab = "Time (months)", ylab = "PFS probability",
        title = "Carrol 2023 - EAC\nwith immunotherapy (P25)",
        ggtheme = ggplot_theme_bw(), surv.median.line = "v", 
        legend.labs = c("High", "Low"), 
        legend.title = paste0("ASE iNMDeff"),
        legend = "top") 

plot_E_to_save <- ggsave_workaround(plot_E)

################################
############ PLOT F ############
################################

input_figureF <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig29/panel_D_F.RData")

surv_obj <- input_figureF[["Motzer_RCC_PFS"]][["ASE_PTC_NMD_triggering_0.2"]]$surv_obj
km_fit <- input_figureF[["Motzer_RCC_PFS"]][["ASE_PTC_NMD_triggering_0.2"]]$km_fit
df <- input_figureF[["Motzer_RCC_PFS"]][["ASE_PTC_NMD_triggering_0.2"]]$df
# Correlation
cor.test(df$NMDeff,df$PFS_P, method = "pearson")

plot_F <- ggsurvplot(km_fit, 
        data = df, risk.table = FALSE,
        pval = TRUE, conf.int = FALSE,
        ylim = c(0, 1), #xlim = c(0, 22), 
        xlab = "Time (months)", ylab = "PFS probability",
        title = "Motzer 2020 - RCC\nwith immunotherapy (P45)",
        ggtheme = ggplot_theme_bw(), surv.median.line = "v", 
        legend.labs = c("High", "Low"), 
        legend.title = paste0("ASE iNMDeff"),
        legend = "top") 

plot_F_to_save <- ggsave_workaround(plot_F)

################################################
############ PANEL UP --> PLOTS D-F ############
################################################

fig_bottom <- plot_grid(plot_D_to_save, plot_E_to_save, plot_F_to_save, nrow = 1, ncol = 3,  
                        labels = c("D","E","F"), label_size = 20)

########################################
############ FINAL SUPP FIG ############
########################################

SuppFig_final <- plot_grid(fig_up, fig_bottom, nrow = 2, ncol = 1)

final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/SuppFig/SuppFig29_complete.png"
ggsave(final_figure_path, SuppFig_final, width = 260, height = 175, units = "mm")
final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/SuppFig/SuppFig29_complete.pdf"
ggsave(final_figure_path, SuppFig_final, width = 260, height = 175, units = "mm")
