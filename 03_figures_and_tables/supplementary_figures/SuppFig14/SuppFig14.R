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
library(pheatmap)
library(ggplotify)

source("/home/gpalou/projects/NMD/scripts/NMD_efficiency/Figures/ggplot_themes.R")

################################
############ PLOT A ############
################################

input_figureA <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig14/panel_A.RData")

plot_A <- ggplot(data = input_figureA, aes(x = as.factor(PCs), y = correlation, fill = alpha)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_brewer(palette = "YlOrBr", direction = -1,guide = guide_legend(nrow = 1)) +
    labs(fill = "Alpha", y = "1% percentile of autocorrelations", 
        x = "Number of requested PCs", title = "Sparse-PCA parameters optimization") +
    geom_text(aes(label = eff_num_PCs), position = position_dodge(width = 0.9), 
                vjust = -0.25, size = 2, angle = 45) +
    ggplot_theme() + 
        theme(
            plot.margin = unit(c(0, 2, 0, 2), "cm"),
            legend.text = element_text(size = 9)
        )

####################################################
############ PANEL UP --> PLOTS A ##################
####################################################

fig_up <- plot_grid(plot_A, nrow = 1, ncol = 1, labels = c("A"))

################################
############ PLOT B ############
################################

input_figureB <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig14/panel_B.RData")
cols <- grep("PC",colnames(input_figureB$TCGA_CNV_PCA_var_df))

plot_B_tmp <- pheatmap(as.matrix((input_figureB$TCGA_CNV_PCA_var_df[,cols])),
        main = "Somatic pan-cancer CNA gain signatures",
        fontsize_col = 4,
        annotation_row = input_figureB$color_df,
        annotation_colors = input_figureB$color_list,
        show_colnames = TRUE,
        show_rownames = FALSE,
        cellwidth=4, 
        cellheight=0.009,
        cluster_rows = FALSE,
        na_col = "#DDDDDD",
        cluster_cols = FALSE,
        fontsize = 7)

# Convert to ggplot
plot_B <- as.ggplot(plot_B_tmp$gtable) + theme_classic() + 
        theme(
            axis.line = element_blank(),       # Remove axis lines
            axis.text.x = element_blank(),     # Remove X axis text
            axis.text.y = element_blank(),     # Remove Y axis text
            axis.ticks = element_blank(),      # Remove axis ticks
            axis.title.x = element_blank(),    # Remove X axis title
            axis.title.y = element_blank(),    # Remove Y axis title
            legend.text = element_text(size = 10),
            plot.margin = unit(c(1,0, 1, 0), "cm")
            )

################################
############ PLOT C ############
################################

input_figureC <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig14/panel_C.RData")
cols <- grep("PC",colnames(input_figureC$TCGA_CNV_PCA_var_df))

plot_C_tmp <- pheatmap(as.matrix((input_figureC$TCGA_CNV_PCA_var_df[,cols])),
        main = "Somatic pan-cancer CNA deletion signatures",
        fontsize_col = 4,
        annotation_row = input_figureC$color_df,
        annotation_colors = input_figureC$color_list,
        show_colnames = TRUE,
        show_rownames = FALSE,
        cellwidth=4, 
        cellheight=0.009,
        cluster_rows = FALSE,
        na_col = "#DDDDDD",
        cluster_cols = FALSE,
        fontsize = 7)

# Convert to ggplot
plot_C <- as.ggplot(plot_C_tmp$gtable) + theme_classic() + 
        theme(
            axis.line = element_blank(),       # Remove axis lines
            axis.text.x = element_blank(),     # Remove X axis text
            axis.text.y = element_blank(),     # Remove Y axis text
            axis.ticks = element_blank(),      # Remove axis ticks
            axis.title.x = element_blank(),    # Remove X axis title
            axis.title.y = element_blank(),    # Remove Y axis title
            legend.text = element_text(size = 10),
            plot.margin = unit(c(1,0, 1, 0), "cm")
            )

#####################################################
############ PANEL MID --> PLOTS B ##################
#####################################################

fig_mid <- plot_grid(plot_B, nrow = 1, ncol = 1, labels = c("B"))

########################################################
############ PANEL BOTTOM --> PLOTS C ##################
########################################################

fig_bottom <- plot_grid(plot_C, nrow = 1, ncol = 1, labels = c("C"))

# final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/SuppFig/SuppFig14_complete.png"
# ggsave(final_figure_path, fig_bottom, width = 300, height = 350, units = "mm")

########################################
############ FINAL SUPP FIG ############
########################################

SuppFig_final <- plot_grid(fig_up, fig_mid, fig_bottom, nrow = 3, ncol = 1) + theme_classic()

final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/SuppFig/SuppFig14_complete.png"
ggsave(final_figure_path, SuppFig_final, width = 200, height = 285, units = "mm")
final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/SuppFig/SuppFig14_complete.pdf"
ggsave(final_figure_path, SuppFig_final, width = 200, height = 285, units = "mm")
