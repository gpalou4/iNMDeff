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

input_figureA <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig8/panel_A.RData")

# Organ system
input_figureA$pan_organ_system <- "Other"
input_figureA[input_figureA$GTEx_tissues %in% c("KDNCTX","KDNMDL"),"pan_organ_system"] <- "Pan-kidney"
input_figureA[input_figureA$GTEx_tissues %in% c("CLNSGM","CLNTRN","ESPGEJ","ESPMCS","ESPMSL","STMACH","SNTTRM"),"pan_organ_system"] <- "Pan-GI"
input_figureA[input_figureA$GTEx_tissues %in% c("BREAST","UTERUS","TESTIS","PRSTTE","VAGINA","OVARY","UCS","CVXECT","CVSEND","FLLPNT"),"pan_organ_system"] <- "Pan-reproductive"
input_figureA[input_figureA$GTEx_tissues %in% c("NERVET","BRNAMY","BRNACC","BRNCDT","BRNCHB","BRNCHA","BRNCTXA","BRNCTXB","BRNHPP","BRNHPT","BRNNCC","BRNPTM","BRNSPC","BRNSNG"),"pan_organ_system"] <- "Pan-nervous"
input_figureA[input_figureA$GTEx_tissues %in% c("LUNG","BLDDER"),"pan_organ_system"] <- "Pan-squamous"
table(input_figureA$pan_organ_system)

# Create a named vector of colors
# accent_colors <- brewer.pal(n = length(unique(input_figureA$pan_organ_system)), name = "Accent")[length(unique(input_figureA$pan_organ_system)):1]  # Reverse the colors with [length(...):1]
# my_colors <- setNames(accent_colors, unique(input_figureA$pan_organ_system))
my_colors <- c(
  "Other" = "#F0027F",
  "Pan-GI" = "#FDC086",
  "Pan-kidney" = "#7FC97F",
  "Pan-nervous" = "#BEAED4",
  "Pan-reproductive" = "#FFFF99",
  "Pan-squamous" = "#386CB0"
)

formula <- y ~ x

# Plot
plot_A_1 <- ggplot(data = input_figureA, 
                   mapping = aes(x = GTEx_ETG_iNMDeff_method_ranking, 
                                 y = GTEx_ASE_Teral_et_al_ranking)) +
    geom_point(aes(color = pan_organ_system), alpha = 1) +
    geom_text_repel(aes(label = GTEx_tissues), 
                    color = "black",  # Set color outside of aes()
                    size = 2.5, nudge_y = 0.05, max.overlaps = nrow(input_figureA)) +
    geom_smooth(method = "lm", formula = formula, se = FALSE, size = 1) +
    labs(x = "ETG iNMDeff GTex tissues ranking", 
         y = "Teran et al. tissues ranking", 
         title = "", 
         color = "Pan-organ system") +  # Add legend title for color
    ggplot_theme() +
    scale_color_manual(values = my_colors) + 
    theme(
        legend.position = "top"
    ) + 
    guides(color = guide_legend(override.aes = list(size = 5))) +  # Increase point size in the legend
    stat_cor(p.accuracy = 0.01, r.accuracy = 0.01, size = 4)

plot_A_2 <- ggplot(data = input_figureA, mapping = aes(x = GTEx_ASE_iNMDeff_method_ranking, 
        y = GTEx_ASE_Teral_et_al_ranking)) +
    geom_point(aes(color = pan_organ_system), alpha = 1) +
    geom_text_repel(aes(label = GTEx_tissues), 
                    color = "black",  # Set color outside of aes()
                    size = 2.5, nudge_y = 0.05, max.overlaps = nrow(input_figureA)) +
    geom_smooth(method = "lm", formula = formula, se = FALSE, size = 1) +
    scale_color_manual(values = my_colors) + 
    labs(x = "ASE iNMDeff GTex tissues ranking", y = "Teran et al. tissues ranking", title = "") +
    ggplot_theme() +
    theme(legend.position='top'
    #   plot.margin = unit(c(0, 2, 0, 2), "cm")
      ) + 
    guides(color = guide_legend(override.aes = list(size = 5))) +  # Increase point size in the legend
    stat_cor(p.accuracy = 0.01, r.accuracy = 0.01, size = 4)



######################################################
############ PANEL up --> PLOTS A ####################
######################################################

# fig_up <- plot_grid(plot_A, plot_B, nrow = 1, ncol = 2, label_size = 20)
fig_up <- ggarrange(plot_A_1,plot_A_2, labels = c("A","B"),
                font.label = list(size = 20),
                ncol=2, nrow=1, common.legend = TRUE)

########################################
############ FINAL SUPP FIG ############
########################################

SuppFig_final <- plot_grid(fig_up, nrow = 1, ncol = 1, rel_heights = c(0.5,0.5)) + theme_classic()

final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/SuppFig/SuppFig8_complete.png"
ggsave(final_figure_path, SuppFig_final, width = 225, height = 100, units = "mm")
final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/SuppFig/SuppFig8_complete.pdf"
ggsave(final_figure_path, SuppFig_final, width = 225, height = 100, units = "mm")
