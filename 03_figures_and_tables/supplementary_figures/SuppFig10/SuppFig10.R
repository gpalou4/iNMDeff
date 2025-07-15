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

input_figureA <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig10/panel_A.RData")

# Organ system
input_figureA$pan_organ_system <- "Other"
input_figureA[input_figureA$tissues %in% c("KDNCTX","KDNMDL"),"pan_organ_system"] <- "Pan-kidney"
input_figureA[input_figureA$tissues %in% c("CLNSGM","CLNTRN","ESPGEJ","ESPMCS","ESPMSL","STMACH","SNTTRM"),"pan_organ_system"] <- "Pan-GI"
input_figureA[input_figureA$tissues %in% c("BREAST","UTERUS","TESTIS","PRSTTE","VAGINA","OVARY","UCS","CVXECT","CVSEND","FLLPNT"),"pan_organ_system"] <- "Pan-reproductive"
input_figureA[input_figureA$tissues %in% c("NERVET","BRNAMY","BRNACC","BRNCDT","BRNCHB","BRNCHA","BRNCTXA","BRNCTXB","BRNHPP","BRNHPT","BRNNCC","BRNPTM","BRNSPC","BRNSNG"),"pan_organ_system"] <- "Pan-nervous"
input_figureA[input_figureA$tissues %in% c("LUNG","BLDDER"),"pan_organ_system"] <- "Pan-squamous"
table(input_figureA$pan_organ_system)

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
  ungroup() #%>%
  # select(-median_score, -closest)  # Optionally remove the helper columns

plot_A <- input_figureA %>%
    ggplot(aes(x = tissues, y = values, fill = pan_organ_system)) +
    geom_violin(draw_quantiles = TRUE, na.rm = TRUE) + 
    coord_flip(ylim = c(-0.5,0.6)) +
    facet_wrap(. ~ NMD_method) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 0.5) +
    labs(y = "Tissue iNMDeff deviation", x = "", fill = "Pan-organ system", title = "GTex") +
    ggplot_theme() +
    theme(legend.position = "top") +
    scale_fill_manual(values = my_colors) +
    geom_text(aes(label = ifelse(p_value_to_show < 0.05, "*", "")), 
                position = position_dodge(), size = 6, color = "black", hjust = 0.5, vjust = 0.4) +
    guides(fill = guide_legend(nrow = 2))

####################################################
############ PANEL UP --> PLOTS A-B ################
####################################################

fig_up <- ggarrange(plot_A, labels = c(""), font.label = list(size = 20),
                ncol=1, nrow=1, common.legend = TRUE, legend="bottom")

########################################
############ FINAL SUPP FIG ############
########################################

SuppFig_final <- plot_grid(fig_up, nrow = 1, ncol = 1) + theme_classic()

final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/SuppFig/SuppFig10_complete.png"
ggsave(final_figure_path, SuppFig_final, width = 150, height = 175, units = "mm")
final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/SuppFig/SuppFig10_complete.pdf"
ggsave(final_figure_path, SuppFig_final, width = 150, height = 175, units = "mm")
