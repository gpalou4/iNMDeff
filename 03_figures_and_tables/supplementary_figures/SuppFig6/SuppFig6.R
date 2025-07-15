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

input_figureA <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig6/panel_A.RData")
input_figureA$pan_organ_system <- "Other"
input_figureA[input_figureA$cancer_type_strat %in% c("KIRC","KIRP"),"pan_organ_system"] <- "Pan-kidney"
input_figureA[grep("COAD|STAD|READ|ESCA",input_figureA$cancer_type_strat),"pan_organ_system"] <- "Pan-GI"
input_figureA[grep("BRCA|OV|UCS|UCEC|CESC|TGCT|PRAD",input_figureA$cancer_type_strat),"pan_organ_system"] <- "Pan-reproductive"
input_figureA[input_figureA$cancer_type_strat %in% c("LGG","PCPG","GBM"),"pan_organ_system"] <- "Pan-nervous"
input_figureA[grep("HNSC|LUSC|BLCA",input_figureA$cancer_type_strat),"pan_organ_system"] <- "Pan-squamous"
table(input_figureA$pan_organ_system)

plot_A_tmp <- ggplot(data = input_figureA, mapping = aes(x = TCGA_ETG_iNMDeff_method_ranking, y = TCGA_ASE_iNMDeff_method_ranking)) +
    geom_point(aes(color = pan_organ_system), alpha = 1) +
    geom_text_repel(aes(label = cancer_type_strat), 
                    color = "black",  # Set color outside of aes()
                    size = 2.5, nudge_y = 0.05, max.overlaps = nrow(input_figureA)) +
    geom_smooth(method = "lm", formula = formula, se = FALSE, size = 1) +
    labs(x = "ETG iNMDeff TCGA tissues ranking", y = "ASE iNMDeff TCGA tissues ranking", title = "", color = "Pan-organ system") +
    ggplot_theme() +
    scale_color_manual(values = my_colors) + 
    theme(legend.position='top',
      plot.margin = unit(c(0, 2, 0, 2), "cm"),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10)
      ) + 
    guides(color = guide_legend(override.aes = list(size = 5))) +  # Increase point size in the legend
    stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, size = 4, label.y = 50)

empty_plot <- plot_grid(NULL, NULL, labels = "")

plot_A <- plot_grid(empty_plot, plot_A_tmp, empty_plot, nrow = 3, labels = c("A","",""), label_size = 20, 
                        rel_heights = c(0.04,1,0.04), vjust = 1, label_y = 1.5)

################################
############ PLOT B ############
################################

input_figureC <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig6/panel_B.RData")

# Change names of columns
cols <- c("ETG NMD Kar", "ETG NMD Col", "ETG NMD Tani", "ETG RndGn w/ NMD feat", "ETG NMD Court",
    "ETG NMD Ens", "ETG NMD All", "ETG NMD Cons", "ETG NMD SMG6", "ETG NMD SMG7", "ETG RndGn w/o NMD feat", "Leukocyte %", "Purity",
    "ETG # NMD targ", "ASE Syn", "ASE NMD-trig PTCs", "ASE NMD-ev PTCs", "ASE # PTCs", "TMB", "TIB", "TNB",
    "CNAburden", "SampleLibSize", "MSIscore","DaysToDeath","Age", "Sex")
colnames(input_figureC) <- cols
rownames(input_figureC) <- cols

input_figureC_filt <- input_figureC[,c("ETG NMD Cons","ASE NMD-trig PTCs")]

plot_B_1 <- ggcorrplot(as.matrix(input_figureC_filt),# method = "square", type = "upper",
                    pch = 10,
                    title = "",
                    legend.title = "Correlation",
                    tl.cex = 7,
                    tl.srt = 30,
                    lab = TRUE,
                    lab_size = 2.2,
                    pch.cex = 15,
                    digits = 1,
                    ggtheme = ggplot_theme(),
                    insig = "pch",
                    hc.order = FALSE,
                    show.diag = NULL)  + ggtitle ("TCGA")

plot_B_2 <- ggcorrplot(as.matrix(input_figureC), method = "square", type = "upper",
                    pch = 10,
                    title = "",
                    legend.title = "Correlation",
                    show.legend = FALSE,
                    tl.cex = 8,
                    tl.srt = 55,
                    # lab = TRUE,
                    lab_size = 2,
                    pch.cex = 15,
                    digits = 1,
                    ggtheme = ggplot_theme(),
                    insig = "pch",
                    hc.order = TRUE,
                    show.diag = NULL)  + ggtitle ("") #+
                    # theme(
                    #     plot.margin = unit(c(-1.5, -1, 0, -1), "cm")
                    #     ) #+ theme_classic()

plot_B <- plot_grid(plot_B_1, plot_B_2, nrow = 2, ncol = 1, labels = c("", ""), 
                    label_size = 20, vjust = 1, rel_heights = c(0.3,0.7))

# final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/SuppFig/test.png"
# ggsave(final_figure_path, plot = plot_A, width = 350, height = 350, units = "mm") 

################################
############ PLOT C ############
################################

input_figureC <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig6/panel_C.RData")

# Change names of columns
# cols <- c("ETG NMD Karousis", "ETG NMD Colombo", "ETG NMD Tani", "ETG RndGn w/ NMD feat.", "ETG NMD Courtney",
#     "ETG NMD Ensembl", "ETG NMD All", "ETG NMD Consensus", "ETG NMD SMG6", "ETG NMD SMG7", "ETG RndGn w/o NMD feat.",
#     "ETG num NMD targets", "ASE Synonymous", "ASE NMD-triggering PTCs", "ASE NMD-evading PTCs", "ASE num PTCs",
#     "Sample lib size", "Sex", "Age", "Death group")
cols <- c("ETG NMD Kar", "ETG NMD Col", "ETG NMD Tani", "ETG RndGn w/ NMD feat.", "ETG NMD Court",
    "ETG NMD Ens", "ETG NMD All", "ETG NMD Cons", "ETG NMD SMG6", "ETG NMD SMG7", "ETG RndGn w/o NMD feat.",
    "ETG # NMD targ", "ASE Syn", "ASE NMD-trig PTCs", "ASE NMD-ev PTCs", "ASE # PTCs",
    "SampleLibSize", "Sex", "Age", "DeathGroup")
colnames(input_figureC) <- cols
rownames(input_figureC) <- cols

input_figureC_filt <- input_figureC[,c("ETG NMD Cons","ASE NMD-trig PTCs")]

plot_C_1 <- ggcorrplot(as.matrix(input_figureC_filt),# method = "square", type = "upper",
                    pch = 10,
                    title = "",
                    legend.title = "Correlation",
                    tl.cex = 7,
                    tl.srt = 30,
                    lab = TRUE,
                    lab_size = 2.8,
                    pch.cex = 15,
                    digits = 1,
                    ggtheme = ggplot_theme(),
                    insig = "pch",
                    hc.order = FALSE,
                    show.diag = NULL)  + ggtitle ("GTex")

plot_C_2 <- ggcorrplot(as.matrix(input_figureC), method = "square", type = "upper",
                    pch = 10,
                    title = "",
                    show.legend = FALSE,
                    legend.title = "Correlation",
                    tl.cex = 8,
                    tl.srt = 50,
                    # lab = TRUE,
                    lab_size = 2,
                    pch.cex = 15,
                    digits = 1,
                    ggtheme = ggplot_theme(),
                    insig = "pch",
                    hc.order = TRUE,
                    show.diag = NULL) + ggtitle ("")# +
                    # theme(
                    #     plot.margin = unit(c(-1.5, -1, 0, -1), "cm")
                    #     ) #+ theme_classic()

# Example rows to highlight
# highlight_rows <- c("ETG NMD Consensus", "ASE NMD-triggering PTCs")

# # Create the plot
# plot_B_2 <- ggcorrplot(
#   as.matrix(input_figureC),
#   method = "square",
#   type = "upper",
#   pch = 10,
#   title = "",
#   legend.title = "Correlation",
#   tl.cex = 9,
#   tl.srt = 30,
#   lab_size = 2,
#   pch.cex = 15,
#   digits = 1,
#   ggtheme = ggplot_theme(),
#   insig = "pch",
#   hc.order = TRUE,
#   show.diag = NULL
# ) + ggtitle("GTex")

# # Add highlights for specific rows
# plot_B_2 <- plot_B_2 +
#   theme(axis.text.y = element_text(
#     size = ifelse(rownames(as.matrix(input_figureC)) %in% highlight_rows, 12, 9),
#     color = ifelse(rownames(as.matrix(input_figureC)) %in% highlight_rows, "red", "black")
#   ))

plot_C <- plot_grid(plot_C_1, plot_C_2, nrow = 2, ncol = 1, labels = c("", ""), 
                    label_size = 20, vjust = 1, rel_heights = c(0.3,0.7))

########################################
############ PANEL UP --> A ############
########################################

fig_up <- plot_grid(empty_plot,plot_A,empty_plot, nrow = 1, ncol = 3, label_size = 20, labels = c(""),
                        rel_widths = c(0.2,0.6,0.2))

# fig_up <- plot_grid(plot_A, nrow = 1, ncol = 1, labels = c("A"), 
#                     label_size = 20, vjust = 1)

##############################################
############ PANEL BOTTOM --> C-D ############
##############################################

fig_bottom <- plot_grid(plot_B, plot_C, nrow = 1, ncol = 2, labels = c("B","C"), 
                    label_size = 20, vjust = 1, rel_widths = c(0.5,0.5))

################################################
############ FINAL FIG --> SUPP FIG ############
################################################

SuppFig_final <- plot_grid(fig_up, fig_bottom, nrow = 2, rel_heights = c(0.35,0.65)) + ggplot_theme()

final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/SuppFig/SuppFig6_complete.png"
ggsave(final_figure_path, SuppFig_final, width = 275, height = 275, units = "mm")
