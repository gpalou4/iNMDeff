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

input_figureA <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig16/panel_A.RData")
input_figureA <- input_figureA[-which(input_figureA$gene_name %in% c("SNORA44","RN7SL371P","RNA5SP162")),]

# Plot
plot_A <- ggplot(data = input_figureA, aes(x = genome_location, y = values, group = ind, color = ind)) +
        geom_point(size = 1) + ggtitle(paste0("")) +
        geom_line() +  
        scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1,1.5,2)) +
        geom_hline(yintercept=c(0,1,2)) +
        labs(color = "CNA-PC86 groups", x = "Chr1 Genome location", y = "GISTIC CNA average state") +
        ggplot_theme() +
        scale_color_brewer(palette="Set2") +
        theme(legend.position = "top") +
        scale_x_discrete(guide=guide_axis(n.dodge = 3)) +
        guides(color = guide_legend(override.aes = list(size = 6), nrow = 1))

################################
############ PLOT B ############
################################

# Plot C
input_figureB <- readRDS("/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig16/panel_B.RData")
combinations_B <- combn(names(table(input_figureB$bins)), 2, simplify = FALSE)

plot_B <- input_figureB %>% 
        filter(NMD_method == "ASE") %>%
        group_by(bins) %>% dplyr::mutate(N=n()) %>%
        dplyr::mutate(N = ifelse(NMDeff == NMDeff[which.min(abs(NMDeff - median(NMDeff)))],paste0('',N),NA)) %>%
        dplyr::mutate(NMDeff = ifelse(is.na(N),NMDeff,NMDeff+0.25)) %>%
        ggplot(aes(x = factor(bins), y = NMDeff, fill = factor(bins),label = as.character(N))) +
                geom_violin() + coord_cartesian(ylim = c(-2,2)) + 
                geom_boxplot(width=0.3, color="black", alpha=0.2) +
                #facet_wrap( ~ NMD_method, scales = "free") + 
                labs(title = "", x = "CNA-PC3 groups", y = "ASE iNMDeff") +
                geom_text(size = 4) +
                scale_x_discrete(labels = c("High", "Mid", "Low")) +
                ggplot_theme() +
                scale_fill_brewer(palette = "Dark2") +
                theme(legend.position = "none",
                        plot.margin = unit(c(1, 0, 1, 0), "cm")) +
                stat_compare_means(comparisons = combinations_B, size = 4,
                                label.y = c(0.5,1,1.35),
                                label = "p.format", method = "wilcox.test", hide.ns = TRUE)

################################
############ PLOT C ############
################################

# Plot D1
input_figureC_1 <- readRDS("/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig16/panel_C.RData")
combinations_C_1 <- combn(names(table(input_figureC_1$bins)), 2, simplify = FALSE)

plot_C_1 <- input_figureC_1 %>% 
        filter(NMD_method == "Endogenous") %>%
        group_by(bins) %>% dplyr::mutate(N=n()) %>%
        dplyr::mutate(N = ifelse(NMDeff == NMDeff[which.min(abs(NMDeff - median(NMDeff)))],paste0('',N),NA)) %>%
        dplyr::mutate(NMDeff = ifelse(is.na(N),NMDeff,NMDeff+0.25)) %>%
        ggplot(aes(x = factor(bins), y = NMDeff, fill = factor(bins),label = as.character(N))) +
                geom_violin() + coord_cartesian(ylim = c(-2,2)) + 
                geom_boxplot(width=0.3, color="black", alpha=0.2) +
                #facet_wrap( ~ NMD_method, scales = "free") + 
                labs(title = "", x = "CNA-PC86 groups", y = "ETG iNMDeff") +
                geom_text(size = 4) +
                scale_x_discrete(labels = c("High", "Mid", "Low")) +
                ggplot_theme() +
                scale_fill_brewer(palette = "Dark2") +
                theme(legend.position = "none",
                         plot.margin = unit(c(0.5, 1, 0.5, 1), "cm")) +
                stat_compare_means(comparisons = combinations_C_1, size = 4,
                                label.y = c(-0.25,0.25,0.75),
                                label = "p.format", method = "wilcox.test", hide.ns = TRUE)

# Plot D2
input_figureC_2 <- readRDS("/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig16/panel_C.RData")
combinations_C_2 <- combn(names(table(input_figureC_2$bins)), 2, simplify = FALSE)

plot_C_2 <- input_figureC_2 %>% 
        filter(NMD_method == "ASE") %>%
        group_by(bins) %>% dplyr::mutate(N=n()) %>%
        dplyr::mutate(N = ifelse(NMDeff == NMDeff[which.min(abs(NMDeff - median(NMDeff)))],paste0('',N),NA)) %>%
        dplyr::mutate(NMDeff = ifelse(is.na(N),NMDeff,NMDeff+0.25)) %>%
        ggplot(aes(x = factor(bins), y = NMDeff, fill = factor(bins),label = as.character(N))) +
                geom_violin() + coord_cartesian(ylim = c(-2,2)) + 
                geom_boxplot(width=0.3, color="black", alpha=0.2) +
                #facet_wrap( ~ NMD_method, scales = "free") + 
                labs(title = "", x = "CNA-PC86 groups", y = "ASE iNMDeff") +
                geom_text(size = 4) +
                scale_x_discrete(labels = c("High", "Mid", "Low")) +
                ggplot_theme() +
                scale_fill_brewer(palette = "Dark2") +
                theme(legend.position = "none",
                        plot.margin = unit(c(0.5, 1, 0.5, 1), "cm")) +
                stat_compare_means(comparisons = combinations_C_2, size = 4,
                                label.y = c(0.5,1,1.35),
                                label = "p.format", method = "wilcox.test", hide.ns = TRUE)

plot_C <- plot_grid(plot_C_1, plot_C_2, nrow = 1, ncol = 2)

######################################################
############ PANEL UP --> PLOTS A-B ##################
######################################################

fig_up <- plot_grid(plot_A, plot_B, nrow = 1, ncol = 2, labels = c("A","B"))

######################################################
############ PANEL BOTTOM --> PLOTS C ################
######################################################

fig_bottom <- plot_grid(plot_C, nrow = 1, ncol = 1, labels = c("C"))

########################################
############ FINAL SUPP FIG ############
########################################

SuppFig_final <- plot_grid(fig_up, fig_bottom, nrow = 2, ncol = 1, rel_heights = c(0.55,0.45))

final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/SuppFig/SuppFig16_complete.png"
ggsave(final_figure_path, SuppFig_final, width = 185, height = 175, units = "mm")

