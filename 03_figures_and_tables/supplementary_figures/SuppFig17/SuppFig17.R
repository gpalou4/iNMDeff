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

input_figureA <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig17/panel_A.RData")
input_figureA$dataset <- ifelse(input_figureA$dataset == "GTEx","GTex","TCGA")
# Plot
plot_A <- input_figureA %>% 
        filter(NMDeff_method == "endogenous_NMD_global_2_shared") %>%
            ggplot(aes(x = model, y = values, fill = type)) +
            facet_wrap(~ dataset) + scale_fill_brewer(palette = "Paired", direction = -1) +
            geom_bar(stat = "identity", position = "dodge", color = "black") +
            labs(x = "Model", y = "R-square", fill = "") +
            ggplot_theme() +
            theme(legend.position = "top")

################################
############ PLOT B ############
################################

input_figureB <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig17/panel_B.RData")
combinations_B <- combn(names(table(input_figureB$Chr_1q_CNV)), 2, simplify = FALSE)

plot_B <- ggplot(input_figureB, aes(x = Chr_1q_CNV, y = -(NMDeff), fill = Chr_1q_CNV)) +
                geom_violin() + scale_fill_brewer(palette = "Dark2", direction = 1) + 
                geom_boxplot(width=0.3, color="black", alpha=0.2) +
                scale_x_discrete(labels = c("Gain","Neutral","Deletion")) +
                labs( title = "Validation in cell lines",x = "CNA state - Chr 1q", y = "GTex-model ETG cNMDeff") +
                ggplot_theme() + coord_cartesian(ylim = c(-2,2)) +
                theme(legend.position = "none",
                        #plot.title = element_text(hjust = 0.90),
                        #plot.margin = unit(c(1, 0, 1, 0), "cm")
                        ) +
                stat_compare_means(comparisons = combinations_B, size = 4,
                                label.y = c(0.25,0.75,1.25),
                                label = "p.format", method = "wilcox.test", hide.ns = TRUE) +
                annotate("text",
                        #fontface = "bold",
                        x = 1:length(table(input_figureB$Chr_1q_CNV)),
                        y = aggregate( -(NMDeff) ~ Chr_1q_CNV, input_figureB, median)[ , 2],
                        label = table(input_figureB$Chr_1q_CNV),
                        col = "black",
                        vjust = - 0.5,
                        size = 4)

######################################################
############ PANEL UP --> PLOTS A-B ##################
######################################################

fig_up <- plot_grid(plot_A, plot_B, nrow = 1, ncol = 2, labels = c("A","B"))

########################################
############ FINAL SUPP FIG ############
########################################

SuppFig_final <- plot_grid(fig_up, nrow = 1, ncol = 1, rel_heights = c(0.5))

final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/SuppFig/SuppFig17_complete.png"
ggsave(final_figure_path, SuppFig_final, width = 180, height = 80, units = "mm")

