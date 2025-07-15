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
library(ggh4x)

source("/home/gpalou/projects/NMD/scripts/NMD_efficiency/Figures/ggplot_themes.R")

ggplot_theme_bw <- function() {

  theme_bw() +
  theme(

    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.title.align = 0.5,
    legend.margin = margin(0,0,0,0),
    legend.box.margin = margin(-5,0,0,0),
    legend.text = element_text(colour = "black", size = 8),

    axis.text = element_text(colour = "black", size = 10),
    axis.title = element_text(colour = "black", size = 12, hjust = 0.5),

    plot.title = element_text(size = 14, hjust = 0.5),
     strip.background = element_rect(fill = "#f4f1f8e4", colour = "#000000"),

    strip.text = element_text(colour = "black", size = 12)
  )
}

################################
############ PLOT A ############
################################

input_figureA <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig23/panel_A.RData")

input_figureA$database_discovery <- gsub("GTEx","GTex",input_figureA$database_discovery)
input_figureA$database_validation <- gsub("GTEx","GTex",input_figureA$database_validation)


# Plot
plot_A <- input_figureA %>%
        filter(!matched_tissue %in% c("pantissue","pancancer")) %>%
        group_by(FDR_threshold_used, database_discovery, database_validation, randomization) %>%
        summarise(num_replicated_hits = n()) %>%
        mutate(database_validation = paste0(database_validation," - validation")) %>%
        mutate(database_discovery = paste0(database_discovery," - discovery")) %>% data.frame() %>%
            ggplot(aes(x = factor(FDR_threshold_used), y = num_replicated_hits, fill = randomization)) + 
                geom_bar(stat='identity',position=position_dodge()) +
                facet_grid(. ~ database_discovery + database_validation) +
                scale_fill_brewer(palette = "Set2", direction = -1, labels = c("Observed", "Randomization")) +
                labs(title = "", x = "% FDR threshold", y = "Replicated hits", fill = "") + 
                ggplot_theme_bw() +
                scale_x_discrete(labels=c("1","2","3","4","5","10")) +
                theme(legend.position = "top",
                    legend.text = element_text(size = 10),
                    plot.margin = unit(c(3, 0, 3, 0), "cm"),
                    legend.margin = margin(t = 0, r = 0, b = 0, l = 0)) #+
                # guides(fill = guide_legend(override.aes = list(size = 10)))

################################
############ PLOT B ############
################################

input_figureB <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig23/panel_B.RData")

# Order
order_vector <- input_figureB[order(input_figureB$matched_tissue),"gene_tissue"]
input_figureB$gene_tissue <- factor(input_figureB$gene_tissue, levels = unique(order_vector))

input_figureB$facet_order <- gsub("GTEx","GTex",input_figureB$facet_order)

plot_B <- input_figureB %>%
        mutate(NMD_method = if_else( NMD_method == "Endogenous","ETG","ASE")) %>%
        filter(database_discovery == "GTEx" & database_validation == "TCGA") %>%
        group_by(disc_val) %>%
            ggplot(aes(x = dataset,y = gene_tissue)) +
                geom_tile(aes(fill = n_carriers)) +
                geom_text(aes(label = n_carriers),color = "black") +
                facet_nested(. ~ NMD_method + facet_order) +
                ggplot_theme_bw() + 
                theme(axis.text.x = element_text(angle = 45, hjust = 1),
                        legend.position = "right") +
                scale_fill_gradientn(colours = c('grey','#F21A00')) +
                labs(title = "Replicated hits at a 5% FDR", x = "", y = "", fill = '# of individuals')

######################################################
############ PANEL UP --> PLOTS A-B ##################
######################################################

fig_up <- plot_grid(plot_A, plot_B, nrow = 1, ncol = 2, labels = c("A","B"), rel_widths = c(0.3,0.7))

########################################
############ FINAL SUPP FIG ############
########################################

SuppFig_final <- plot_grid(fig_up, nrow = 1, ncol = 1, rel_heights = c(0.5))

final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/SuppFig/SuppFig23_complete.png"
ggsave(final_figure_path, SuppFig_final, width = 300, height = 150, units = "mm")

