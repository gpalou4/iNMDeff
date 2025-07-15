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
library(ggh4x)
library(scales)

source("/home/gpalou/projects/NMD/scripts/NMD_efficiency/Figures/ggplot_themes.R")

################################
############ PLOT A ############
################################

input_figureA <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig12/panel_A_B.RData")

plot_A <- ggplot(input_figureA, aes(x = cancer_type, y = lambda_non_random, color = factor(dataset))) + 
        geom_point(size = 3) +
        labs(x = "", y = "Lambda", title = "Cancer and NMD genes", color = "Somatic mutations") +
        ylim(c(0,4)) +
        geom_hline(yintercept = 1.5, linetype = "dashed", color = "red", size = 2) +
        facet_grid(. ~ NMD_method) +
        ggplot_theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust=1),
                legend.title = element_text(size = 11),
                legend.text = element_text(colour = "black", size = 10),
            legend.position = "top")

####################################################
############ PANEL UP --> PLOTS A ##################
####################################################

fig_up <- plot_grid(plot_A, labels = c("A"), nrow = 1, ncol = 1)

################################
############ PLOT B ############
################################

input_figureB <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig12/panel_A_B.RData")

plot_B <- ggplot(input_figureB, aes(x = cancer_type, y = lambda_random, color = factor(dataset))) + 
        geom_point(size = 3) +
        labs(x = "", y = "Lambda", title = "Random genes", color = "Somatic mutations") +
        ylim(c(0,4)) +
        geom_hline(yintercept = 1.5, linetype = "dashed", color = "red", size = 2) +
        facet_grid(. ~ NMD_method) +
        ggplot_theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust=1),
                legend.title = element_text(size = 11),
                legend.text = element_text(colour = "black", size = 10),
            legend.position = "top")

########################################################
############ PANEL MIDDLE --> PLOTS B ##################
########################################################

fig_middle <- plot_grid(plot_B, labels = c("B"), nrow = 1, ncol = 1)

################################
############ PLOT C ############
################################

input_figureC <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig12/panel_C.RData")

plot_C <- input_figureC %>%
        filter(!cancer_type %in% c("PCPG","CESC")) %>%
        group_by(FDR_threshold_used, NMD_method_discovery, NMD_method_validation, dataset, genesets) %>%
        summarise(num_replicated_hits = ( n() / genes_tested ) * 100) %>%
        mutate(NMD_method_validation = paste0(NMD_method_validation," - val")) %>%
        mutate(NMD_method_discovery = paste0(NMD_method_discovery," - disc")) %>% data.frame() %>%
            ggplot(aes(x = factor(FDR_threshold_used), y = num_replicated_hits, fill = genesets)) + 
                geom_bar(stat='identity',position=position_dodge()) +
                facet_grid(NMD_method_validation + NMD_method_discovery ~ dataset ) +
                labs(title = "", x = "% FDR threshold", y = "% of Replicated Hits", fill = "") + 
                ggplot_theme_bw() + scale_x_discrete(labels=c("1","2","3","4","5","10")) +
                theme(legend.position = "top")

################################
############ PLOT D ############
################################

input_figureD <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig12/panel_D.RData")

plot_D <- input_figureD %>%
        filter(!cancer_type %in% c("PCPG","CESC")) %>%
        group_by(disc_val) %>%
            ggplot(aes(x = dataset, y = gene_cancer, fill = beta_coefficient)) +
            geom_tile(size = 2) +
            geom_text(aes(label = round(-beta_coefficient,2)), color = "black", size = 4) +
            ggplot_theme_bw() +
            facet_nested(. ~ facet_order) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1),
                    legend.position = "right") +
            labs(fill = "Effect size", x = "", y = "") +
            scale_fill_gradientn(colours = cols, 
                                    values = rescale(c(-2, -1, 0, 1, 2)),
                                    guide = guide_colorbar(barheight = 8, barwidth = 1.5),
                                    limits=c(-2, 2))

################################
############ PLOT E ############
################################

input_figureE <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig12/panel_E.RData")
input_figureE$NMD_method <- ifelse(input_figureE$NMD_method == "Endogenous","ETG","ASE")
# check the colours
cols <- brewer.pal(n = 5, name = "RdBu")

plot_E <- input_figureE %>%
            ggplot(aes(x = Gene_symbol, y = cancer_type, fill = beta_coefficient)) +
            geom_tile(size = 2) +
            geom_text(aes(label = round(-beta_coefficient,2)), color = "black", size = 4) +
            ggplot_theme_bw() + 
            facet_nested(. ~ NMD_method) +
            theme(panel.spacing = unit(3, "lines"),
                    legend.position = "right") +
            geom_text(aes(label = ifelse(p_value < 0.02, "*", "")), 
                               size = 10, hjust = -2, color = "black") +
            labs(fill = 'Effect size', x = "", y = "") +
            scale_fill_gradientn(colours = cols, na.value = 'white',
                                    guide = guide_colorbar(barheight = 8, barwidth = 1.5),
                                    values = rescale(c(-2, -1, 0, 1, 2)),
                                    limits=c(-2, 2))

########################################################
############ PANEL BOTTOM --> PLOTS C-E ##################
########################################################

plot_C_D <- plot_grid(plot_C, plot_D, nrow = 2, ncol = 1, labels = c("C","D"), rel_heights = c(0.6,0.4))
fig_bottom <- plot_grid(plot_C_D, plot_E, nrow = 1, ncol = 2, labels = c("","E"))

# final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/SuppFig/test.png"
# ggsave(final_figure_path, plot = fig_bottom, width = 400, height = 250, units = "mm")

########################################
############ FINAL SUPP FIG ############
########################################

SuppFig_final <- plot_grid(fig_up, fig_middle, fig_bottom, nrow = 3, ncol = 1, rel_heights = c(0.325,0.325,0.45))

final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/SuppFig/SuppFig12_complete.png"
ggsave(final_figure_path, SuppFig_final, width = 300, height = 350, units = "mm")


