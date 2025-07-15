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

source("/home/gpalou/projects/NMD/scripts/NMD_efficiency/Figures/ggplot_themes.R")

my_colors <- c(
  "Other" = "#F0027F",
  "Pan-GI" = "#FDC086",
  "Pan-kidney" = "#7FC97F",
  "Pan-nervous" = "#BEAED4",
  "Pan-reproductive" = "#FFFF99",
  "Pan-squamous" = "#386CB0"
)

################################
############ PLOT A ############
################################

input_figureA <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig9/panel_A.RData")

# Organ system
input_figureA$pan_organ_system <- "Other"
input_figureA[input_figureA$tissues %in% c("KIRC","KIRP"),"pan_organ_system"] <- "Pan-kidney"
input_figureA[grep("COAD|STAD|READ|ESCA",input_figureA$tissues),"pan_organ_system"] <- "Pan-GI"
input_figureA[grep("BRCA|OV|UCS|UCEC|CESC|TGCT|PRAD",input_figureA$tissues),"pan_organ_system"] <- "Pan-reproductive"
input_figureA[input_figureA$tissues %in% c("LGG","PCPG","GBM"),"pan_organ_system"] <- "Pan-nervous"
input_figureA[grep("HNSC|LUSC|BLCA",input_figureA$tissues),"pan_organ_system"] <- "Pan-squamous"
table(input_figureA$pan_organ_system)

# Plot
# Extract colors from the Accent palette
# magma_colors <- viridis(n = length(unique(input_figureC$pan_organ_system)), direction = -1, option = "viridis")
# accent_colors <- brewer.pal(n = length(unique(input_figureA$pan_organ_system)), name = "Accent")[length(unique(input_figureA$pan_organ_system)):1]  # Reverse the colors with [length(...):1]
# # Create a named vector of colors
# my_colors <- setNames(accent_colors, unique(input_figureA$pan_organ_system))
# Correct tissue p-values by FDR
tissues_p_value <- input_figureA[,c("NMD_method","tissues","p_value_above_NMDeff_median","p_value_below_NMDeff_median")]
tissues_p_value <- unique(tissues_p_value)
tissues_p_value <- tissues_p_value %>%
        rowwise() %>%
        mutate(p_value = if_else(p_value_above_NMDeff_median < p_value_below_NMDeff_median, 
                        p_value_above_NMDeff_median, p_value_below_NMDeff_median)) %>%
        group_by(NMD_method) %>%
        mutate(p_value_FDR_adjust = p.adjust(p_value, method = "fdr"))

# Add adjusted p-values to the dataframe
input_figureA <- merge(input_figureA,tissues_p_value[,c("NMD_method","tissues","p_value","p_value_FDR_adjust")], by = c("NMD_method","tissues"), all.x = TRUE)
input_figureA <- input_figureA
input_figureA <- input_figureA %>%
  group_by(tissues, NMD_method) %>%
  mutate(
    median_score = median(values, na.rm = TRUE),  # Calculate median for each group
    closest = which.min(abs(values - median_score)),  # Find closest row index in each group
    p_value_to_show = ifelse(row_number() == closest, p_value_FDR_adjust, NA)  # Set p_values to NA except for the closest row
  ) %>%
  ungroup() %>%
  select(-median_score, -closest)  # Optionally remove the helper columns

plot_A <- input_figureA %>%
    ggplot(aes(x = tissues, y = values, fill = pan_organ_system)) +
    geom_violin(draw_quantiles = TRUE, na.rm = TRUE) + coord_flip(ylim = c(-0.3,0.3)) +
    facet_wrap(. ~ NMD_method) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 0.5) +
    labs(y = "Tissue iNMDeff deviation", x = "", fill = "Pan-organ system", title = "TCGA") +
    ggplot_theme() +
    theme(legend.position = "top") +
    scale_fill_manual(values = my_colors) +
   geom_text(aes(label = ifelse(p_value_to_show < 0.05, "*", "")), 
                position = position_dodge(), size = 6, color = "black", hjust = 0.5, vjust = 0.4)

################################
############ PLOT B ############
################################

input_figureB <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig9/panel_B.RData")

# Organ system
input_figureB$pan_organ_system <- "Other"
input_figureB[input_figureB$tissues %in% c("KDNCTX","KDNMDL"),"pan_organ_system"] <- "Pan-kidney"
input_figureB[input_figureB$tissues %in% c("CLNSGM","CLNTRN","ESPGEJ","ESPMCS","ESPMSL","STMACH","SNTTRM"),"pan_organ_system"] <- "Pan-GI"
input_figureB[input_figureB$tissues %in% c("BREAST","UTERUS","TESTIS","PRSTTE","VAGINA","OVARY","UCS","CVXECT","CVSEND","FLLPNT"),"pan_organ_system"] <- "Pan-reproductive"
input_figureB[input_figureB$tissues %in% c("NERVET","BRNAMY","BRNACC","BRNCDT","BRNCHB","BRNCHA","BRNCTXA","BRNCTXB","BRNHPP","BRNHPT","BRNNCC","BRNPTM","BRNSPC","BRNSNG"),"pan_organ_system"] <- "Pan-nervous"
input_figureB[input_figureB$tissues %in% c("LUNG","BLDDER"),"pan_organ_system"] <- "Pan-squamous"
table(input_figureB$pan_organ_system)

# Correct tissue p-values by FDR
tissues_p_value <- input_figureB[,c("NMD_method","tissues","p_value_above_NMDeff_median","p_value_below_NMDeff_median")]
tissues_p_value <- unique(tissues_p_value)
tissues_p_value <- tissues_p_value %>%
        rowwise() %>%
        mutate(p_value = if_else(p_value_above_NMDeff_median < p_value_below_NMDeff_median, 
                        p_value_above_NMDeff_median, p_value_below_NMDeff_median)) %>%
        group_by(NMD_method) %>%
        mutate(p_value_FDR_adjust = p.adjust(p_value, method = "fdr"))
# Manual check for manuscript (here there are only 33 tissues, look at raw data...)
# df <- tissues_p_value[tissues_p_value$NMD_method == "ETG",]
# table(df$p_value_FDR_adjust < 0.05)
# Add adjusted p-values to the dataframe
input_figureB <- merge(input_figureB,tissues_p_value[,c("NMD_method","tissues","p_value","p_value_FDR_adjust")], by = c("NMD_method","tissues"), all.x = TRUE)
input_figureB <- input_figureB
input_figureB <- input_figureB %>%
  group_by(tissues, NMD_method) %>%
  mutate(
    median_score = median(values, na.rm = TRUE),  # Calculate median for each group
    closest = which.min(abs(values - median_score)),  # Find closest row index in each group
    p_value_to_show = ifelse(row_number() == closest, p_value_FDR_adjust, NA)  # Set p_values to NA except for the closest row
  ) %>%
  ungroup() %>%
  select(-median_score, -closest)  # Optionally remove the helper columns

plot_B <- input_figureB %>%
    ggplot(aes(x = tissues, y = values, fill = pan_organ_system)) +
    geom_violin(draw_quantiles = TRUE, na.rm = TRUE) + coord_flip(ylim = c(-0.3,0.3)) +
    facet_wrap(. ~ NMD_method) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 0.5) +
    labs(y = "Tissue iNMDeff deviation", x = "", fill = "Pan-organ system", title = "GTex") +
    ggplot_theme() +
    theme(legend.position = "top") +
    scale_fill_manual(values = my_colors) +
   geom_text(aes(label = ifelse(p_value_to_show < 0.05, "*", "")), 
                position = position_dodge(), size = 6, color = "black", hjust = 0.5, vjust = 0.4)

####################################################
############ PANEL UP --> PLOTS A-B ################
####################################################

fig_up <- ggarrange(plot_A, plot_B, labels = c("A","B"), font.label = list(size = 20),
                ncol=2, nrow=1, common.legend = TRUE, legend="bottom")

################################
############ PLOT C ############
################################

input_figureC <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig9/panel_C.RData")

# Create the boxplot
tmp_plot <- ggboxplot(input_figureC, x = "MSI_status_2", y = "endogenous_NMD_Consensus",
                      palette = "jco", color = "MSI_status_2", notch = TRUE,
                      facet.by = "cancer_type", short.panel.labs = TRUE)
# Add jitter with alpha using geom_jitter
tmp_plot <- tmp_plot + 
  geom_jitter(aes(color = MSI_status_2), width = 0.2, alpha = 0.5) 
tmp_plot <- tmp_plot + ggplot_theme_bw() + xlab("") +
            labs(y = "iNMDeff", x = "", title = "ETG") +
            theme(legend.position = "none",
              strip.background = element_blank()) + coord_cartesian(ylim = c(-2,2.75))
# Use only p.format as label. Remove method name.
plot_C <- tmp_plot + stat_compare_means(label = "p.format", method = "wilcox.test", 
            label.y = c(2.5), label.x = 1.5, hide.ns = TRUE, method.args = list(alternative = "greater"))

################################
############ PLOT D ############
################################

input_figureD <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig9/panel_D.RData")
# Create the boxplot
tmp_plot <- ggboxplot(input_figureD, x = "MSI_status_2", y = "ASE_PTC_NMD_triggering_0.2",
                      palette = "jco", color = "MSI_status_2", notch = TRUE,
                      facet.by = "cancer_type", short.panel.labs = TRUE)
# Add jitter with alpha using geom_jitter
tmp_plot <- tmp_plot + 
  geom_jitter(aes(color = MSI_status_2), width = 0.2, alpha = 0.5) 
tmp_plot <- tmp_plot + ggplot_theme_bw() + xlab("") +
            labs(y = "iNMDeff", x = "", title = "ASE") +
            theme(legend.position = "none",
                strip.background = element_blank()) + coord_cartesian(ylim = c(-2,2.5))
# Use only p.format as label. Remove method name.
plot_D <- tmp_plot + stat_compare_means(label = "p.format", method = "wilcox.test", 
            label.y = c(2), label.x = 1.5, hide.ns = TRUE, method.args = list(alternative = "greater"))

########################################################
############ PANEL BOTTOM --> PLOTS C-D ################
########################################################

fig_bottom <- ggarrange(plot_C, plot_D, labels = c("C","D"), font.label = list(size = 20),
                ncol=2, nrow=1, common.legend = TRUE, legend="none")

# final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/SuppFig/test.png"
# ggsave(final_figure_path, plot = fig_bottom, width = 225, height = 250, units = "mm")

########################################
############ FINAL SUPP FIG ############
########################################

SuppFig_final <- plot_grid(fig_up, fig_bottom, nrow = 2, ncol = 1, rel_heights = c(0.7,0.3)) + theme_classic()

final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/SuppFig/SuppFig9_complete.png"
ggsave(final_figure_path, SuppFig_final, width = 250, height = 300, units = "mm")
final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/SuppFig/SuppFig9_complete.pdf"
ggsave(final_figure_path, SuppFig_final, width = 250, height = 300, units = "mm")
