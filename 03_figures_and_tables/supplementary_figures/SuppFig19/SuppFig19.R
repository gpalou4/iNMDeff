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
library(ggforce)

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

    strip.text = element_text(colour = "black", size = 12)
  )
}

################################
############ PLOT A ############
################################

input_figureA <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig19/panel_A.RData")

# Data
input_figureA <- input_figureA[order(input_figureA$genome_location),]
CNV_color <- brewer.pal(n = 2, name = "Paired")
CNV_color <- CNV_color[-3]
# Outlier
input_figureA <- input_figureA[-which(input_figureA$gene_symbols == "SNORA44"),]

# Create a fake color only to output the legend
input_figureA$CNA_type <- "CNA_freq_1"
input_figureA$CNA_type[1:20] <- "CNA_freq_2"

plot_A1 <- ggplot(data = input_figureA, aes(x = as.numeric(as.character(genome_location)), y = CNV_freq_1, fill = factor(CNA_type) )) +
                # geom_bar(stat="identity", width = 1.2) +
                # geom_point()+
                geom_bar(aes(x = as.numeric(as.character(genome_location)), y = CNV_freq_1), stat = "identity", color = CNV_color[1]) +
                geom_bar(aes(x = as.numeric(as.character(genome_location)), y = CNV_freq_2), stat = "identity", color = CNV_color[2]) +
                labs(title = "", x = "", y = "CNA freq") +
                geom_hline(yintercept= 0, size = 0.5, color = "black") +
                ggplot_theme() +
                theme(axis.text.x = element_blank(),
                        axis.line.x = element_blank(),
                        axis.ticks.x = element_blank(),
                        plot.margin = unit(c(0, 0.5, -0.5, 0.25), "cm"),
                        legend.position = "top") +
                scale_fill_manual(values = c("CNA_freq_1" = CNV_color[1], "CNA_freq_2" = CNV_color[2]),
                                labels = c("CNA_freq_1" = "> 0", "CNA_freq_2" = ">= 1.5"),
                                  name = "GISTIC score")

# Plot A-2 (bottom)
# Colors and labels
coAmp_color <- brewer.pal(n = 8, name = "OrRd")
candidate_genes_NMD <- c("SMG5","RBM8A","SF3B4","INTS3")
NMD_related_genes <- c("SMG7")
genes_label <- c("SMG5","SMG7","RBM8A","SF3B4","INTS3","PMF1","BOLA1","CLK2","PRPF3","DENND4B","TPR","PI4KB","USF1","SETDB1","CRTC2","GORAB","ZNF687",
                "SCNM1","POGZ","YY1AP1","LYSMD1","DUSP12","USP21","GON4L","MSTO1","CHTOP","CHTOP","PPOX","VPS45","VPS72")
genes_label <- c("PMF1","GON4L","CRTC2","POGZ","PRPF3","VPS72")
genes_label <- c(candidate_genes_NMD,genes_label,NMD_related_genes)
input_figureA$label_color <- "#000000"
input_figureA[input_figureA$gene_symbols %in% candidate_genes_NMD,"label_color"] <- "#ff0000"
input_figureA[input_figureA$gene_symbols %in% NMD_related_genes,"label_color"] <- "#2015D0"
input_figureA[input_figureA$gene_symbols %in% genes_label & !input_figureA$gene_symbols %in% 
                        c(candidate_genes_NMD,NMD_related_genes),"label_color"] <- "#DCA620"
tmp <- subset(input_figureA, gene_symbols %in% genes_label)
label_color_char <- tmp[!tmp$gene_name %in% c("PMF1","GON4L"),"label_color"]
input_figureA$genome_location <- factor(input_figureA$genome_location, levels = unique(input_figureA$genome_location))
table(input_figureA$label_color)

plot_A2 <- ggplot(data = input_figureA, aes(x = as.numeric(as.character(genome_location)), y = eval(parse(text=paste0("NMDeff_ASE_1"))), fill = as.numeric(coCNV_freq_1),
                        colour = factor(label_color))) +
        geom_point(size = 1.25, shape = 21,  stroke = 0.4, alpha = 0.7) +
        geom_text_repel(data = subset(input_figureA, gene_symbols %in% genes_label[!genes_label %in% c("PMF1","GON4L")]),
                aes(label = gene_symbols), color = label_color_char, size = 2.5, arrow = arrow(length = unit(0.005, 'npc'), type = "closed"),
                        point.padding = 0.005, nudge_x = .01, nudge_y = .003, max.overlaps = 200) +
        geom_text_repel(data = subset(input_figureA, gene_symbols %in% c("PMF1","GON4L")), fontface = "bold",
                aes(label = gene_symbols), color = c("#DCA620","#DCA620"), size = 2.5, arrow = arrow(length = unit(0.05, 'npc'), type = "closed"),
                        nudge_x = .09, nudge_y = -.006, max.overlaps = 200) +
        labs(title = "", x = "Chr 1 genome location", y = "median ASE iNMDeff") +
        geom_hline(yintercept= 0, size = 0.5, color = "black") +
        # 1q21.1-23.1 bands
        geom_vline(xintercept = c(143955364,147242641), size = 0.5, color = "black", linetype = "dashed") + # q21.1
        geom_vline(xintercept = c(147258885,150574551), size = 0.5, color = "black", linetype = "dashed") + # q21.2
        geom_vline(xintercept = c(150600539,155078872), size = 0.5, color = "black", linetype = "dashed") + # q21.3
        geom_vline(xintercept = c(155127460,156579727), size = 0.5, color = "black", linetype = "dashed") + # q22
        geom_vline(xintercept = c(156594487,158999968), size = 0.5, color = "black", linetype = "dashed") + # q23.1
        scale_x_continuous(
                breaks = c(1.50e+08, 
                        mean(c(150600539,155078872)),
                        mean(c(155127460,156579727)),
                        mean(c(156594487,158999968)),
                        160555568),
                labels = c("q21.2", "q21.3", "q22", "q23.1", "q23.2")
        ) +
        scale_fill_gradient(name = paste0("Co-amp\n with SMG5"), 
                low = coAmp_color[3], high = coAmp_color[8],
                labels = scales::percent_format(scale = 100),
                guide = guide_colorbar(barheight = 1, barwidth = 6)) +
        scale_colour_manual(values = c("#2015D0" = "blue", "#DCA620" = "#DCA620", "#ff0000" = "#ff0000", "#000000" = "black"),
                        labels = c("#DCA620" = "Candidates", "#2015D0" = "NMD-related", "#ff0000" = "Candidates NMD", "#000000" = "Other"),
                                name = "") +
        scale_y_continuous(breaks = c(0,-0.05,-0.1,-0.15, -0.2,-0.25, -0.3)) +
        scale_x_break(breaks = c(121167646, 143955364))+
        ggplot_theme_bw() +
        guides(colour = guide_legend(override.aes = list(size = 4), nrow = 2)) +
        theme(
                plot.margin = unit(c(-0.5, 0.5, 0, 0), "cm"),
                legend.position = "bottom",  # Move legend to the bottom
                legend.text = element_text(size = 8),  # Adjust legend text size if needed
                legend.spacing.x = unit(0.001, "cm"),  # Reduce horizontal spacing
                legend.spacing.y = unit(0.001, "cm"),  # Reduce vertical spacing
                legend.box.margin = margin(1, 1, 1, 1),  # Minimize margin around legend box
                legend.box.spacing = unit(0.001, "cm")  # Minimize spacing between legend items
        )+
        facet_zoom(xlim = c(149887890-100, 161749758+100),
                        zoom.size = 0.5,
                        ylim = c(-0.15, -0.215),
                        show.area = TRUE,
                        horizontal = FALSE)

# Remove gene symbol labels and vertical dashes from the original panel
plot_A2 <- ggplot_build(plot_A2)
# plot_A2$data[[1]][plot_A2$data[[1]]$PANEL %in% c(1,2) & plot_A2$data[[1]]$x %in% c(160555568), 'alpha'] <- 0
for (i in c(2,3,4,5,6,7,8,9)) {
        plot_A2$data[[i]][plot_A2$data[[i]]$PANEL %in% c(1,2), 'alpha'] <- 0
}

plot_A2 <- ggplot_gtable(plot_A2)
plot_A <- plot_grid(plot_A1, plot_A2, nrow = 2, ncol = 1,  
                        rel_heights = c(0.3,0.7), labels = c("",""), label_size = 20)

####################################################
############ PANEL UP --> PLOTS A ##################
####################################################

fig_up <- plot_grid(plot_A, nrow = 1, ncol = 1, labels = c(""))

########################################
############ FINAL SUPP FIG ############
########################################

SuppFig_final <- plot_grid(fig_up, nrow = 1, ncol = 1, rel_heights = c(0.5))

final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/SuppFig/SuppFig19_complete.png"
ggsave(final_figure_path, SuppFig_final, width = 125, height = 150, units = "mm")
final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/SuppFig/SuppFig19_complete.pdf"
ggsave(final_figure_path, SuppFig_final, width = 125, height = 150, units = "mm")




















