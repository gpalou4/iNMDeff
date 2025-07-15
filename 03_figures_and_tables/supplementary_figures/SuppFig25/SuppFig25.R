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

input_figureA <- readRDS(file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig25/panel_A_B.RData")

# Remove GTEx non-European samples and missing PCs samples
outlier_path <- paste0("/g/strcombio/fsupek_cancer1/gpalou/GTEx/samples_metadata/ancestry_and_admixture_and_PCA_outliers.txt")
all_outlier_samples <- read.table(file = outlier_path)$V1
input_figureA <- input_figureA[!input_figureA$sample_short %in% all_outlier_samples,]
input_figureA$NMD_method <- ifelse(input_figureA$NMD_method == "NMD Consensus", "ETG","ASE")
input_figureA$dataset <- ifelse(input_figureA$dataset == "GTEx", "GTex","TCGA")

# Selected genes
selected_genes <- c("NUP153","KDM6B","TRAP1","CRTC1","PXDN","PDIA2","LAMC1","FIG4","COL19A1","CCDC114")
input_figureA_filt <- input_figureA %>%
                        filter(gene %in% selected_genes)
table(input_figureA_filt$gene)

# Add adjusted P-values from the SKAT-O
output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/rare_germline_variants/replication/nogenelist/between/RGWAS_replication_hits.txt")
RGWAS_replication_hits_final <- read.table(file = output_path, header = TRUE, sep = "\t")
RGWAS_replication_selected_hits <- RGWAS_replication_hits_final[RGWAS_replication_hits_final$variant %in% selected_genes & 
                        RGWAS_replication_hits_final$FDR_threshold_used == 2 &
                        RGWAS_replication_hits_final$dataset == "PTV_Missense_CADD15_0.1perc" &
                        RGWAS_replication_hits_final$randomization == "no",]
RGWAS_replication_selected_hits <- RGWAS_replication_selected_hits[,c("variant","database_validation","database_discovery","NMD_method","tissue","pvalue_FDR_adjusted")]
colnames(RGWAS_replication_selected_hits) <- c("gene","database_validation","database_discovery","NMD_method","tissue","pvalue_FDR_adjusted")
RGWAS_replication_selected_hits$NMD_method <- ifelse(RGWAS_replication_selected_hits$NMD_method == "Endogenous","ETG","ASE")
RGWAS_replication_selected_hits <- RGWAS_replication_selected_hits[order(RGWAS_replication_selected_hits$gene),]

# input_figureA_filt$gene <- factor(ifelse(input_figureA_filt$gene == "KDM6B","KDM6B (Thyroid-THCA)","PXDN (LGG-Brain_Substantia_Nigra)"))
genes <- c("CCDC114","COL19A1","CRTC1","FIG4","KDM6B")
plot_A_1 <- input_figureA_filt %>%
        filter(gene %in% genes) %>%
        ggplot(aes(x = germ_mut, y = NMD_efficiency, fill = factor(germ_mut, levels = c("WT","MUT")))) +
                geom_violin(position = position_dodge(width = 0.9)) +
                geom_jitter(aes(fill = germ_mut), 
                        position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0, dodge.width = 0.9),
                        alpha = 0.1, size = 1) +
                geom_boxplot(width = 0.2, position = position_dodge(width = 0.9), color="black", alpha=0.2) +
                coord_cartesian(ylim = c(-2,3)) +
                facet_nested( dataset ~ gene + NMD_method) +
                # new_scale("facet") +
                # facet_wrap(dataset ~ gene + NMD_method, nrow = 2) +
                geom_hline(yintercept = 0, linetype="dashed", color = "red", size=0.5) +
                labs(y = "iNMDeff", x = "", fill = "Rare germline pLoF variant") +
                ggplot_theme_bw() + scale_fill_brewer(palette = "Dark2") +
                theme(legend.position='top',
                        axis.text.x = element_text(size = 8),
                        axis.text.y =  element_text(size = 8),
                        strip.text = element_text(size = 9)) + #strip.background =  element_rect(colour="grey", fill="grey")
                stat_compare_means(data = subset(input_figureA_filt, gene %in% genes & ! ( (gene == "PXDN" & NMD_method == "ASE" 
                                        & dataset == "GTex") | (gene == "PDIA2" & NMD_method == "ASE" & dataset == "TCGA") )),
                                aes(group=germ_mut), size = 2,
                                label.y = c(2.4),
                                label.x = 1, na.rm = TRUE,
                                label = "p.format", method = "wilcox.test", hide.ns = TRUE)
genes <- c("LAMC1","NUP153","PDIA2","PXDN","TRAP1")
plot_A_2 <- input_figureA_filt %>%
        filter(gene %in% genes) %>%
        ggplot(aes(x = germ_mut, y = NMD_efficiency, fill = factor(germ_mut, levels = c("WT","MUT")))) +
                geom_violin(position = position_dodge(width = 0.9)) +
                geom_jitter(aes(fill = germ_mut), 
                        position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0, dodge.width = 0.9),
                        alpha = 0.1, size = 1) +
                geom_boxplot(width = 0.2, position = position_dodge(width = 0.9), color="black", alpha=0.2) +
                coord_cartesian(ylim = c(-2,3)) +
                facet_nested( dataset ~ gene + NMD_method) +
                # new_scale("facet") +
                # facet_wrap(dataset ~ gene + NMD_method, nrow = 2) +
                geom_hline(yintercept = 0, linetype="dashed", color = "red", size=0.5) +
                labs(y = "iNMDeff", x = "", fill = "Rare germline pLoF variant") +
                ggplot_theme_bw() + scale_fill_brewer(palette = "Dark2") +
                theme(legend.position='none',
                        axis.text.x = element_text(size = 8),
                        axis.text.y =  element_text(size = 8),
                        strip.text = element_text(size = 9)) + #strip.background =  element_rect(colour="grey", fill="grey")
                stat_compare_means(data = subset(input_figureA_filt, gene %in% genes & ! ( (gene == "PXDN" & NMD_method == "ASE" 
                                        & dataset == "GTex") | (gene == "PDIA2" & NMD_method == "ASE" & dataset == "TCGA") )),
                                aes(group=germ_mut), size = 2,
                                label.y = c(2.4),
                                label.x = 1, na.rm = TRUE,
                                label = "p.format", method = "wilcox.test", hide.ns = TRUE)
# # Add manually the p-values
# # To know the non-replicated ones, you have to go to the plots_all_genes and check the p.adjustment
# plot_C <- ggplot_build(plot_C)
# # RGWAS_replication_selected_hits[RGWAS_replication_selected_hits$database_discovery == "GTEx" &
# #                                 RGWAS_replication_selected_hits$database_validation == "TCGA",]
# # Should I put these adjusted p-values (based on GTEx discovery --> TCGA validation, or the raw ones?)
# pvals_fixed <- c("Wilcox, p = 0.46\nSKAT-O, p = 1.73-02", # 3.81e-05 (raw SKATO pval)
#                 "Wilcox, p = 0.47\nSKAT-O, ns", # 0.6 (raw SKATO pval), 1 (adjusted)
#                 "Wilcox, p = 0.13\nSKAT-O, p = 7.27-03", # 3.79e-05 (raw SKATO pval)
#                 "Wilcox, p = 0.22\nSKAT-O, p = 5.2-03", # 8.22e-05 (raw SKATO pval)
#                 "Wilcox, p = 0.055\nSKAT-O, ns",  # 0.19 (raw SKATO pval), 0.5 (adjusted)
#                 "Wilcox, p = 0.015\nSKAT-O, ns", # 5.36e-03 (raw SKATO pval)
#                 "Wilcox, p = 0.89\nSKAT-O, p = 1.84-02") # 4.18e-05 (raw SKATO pval)
# plot_C$data[[5]][,"label"] <- pvals_fixed
# plot_C <- ggplot_gtable(plot_C)

# final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/SuppFig/test.png"
# ggsave(final_figure_path, plot_C, width = 400, height = 100, units = "mm")

####################################################
############ PANEL UP --> PLOTS A ##################
####################################################

fig_up <- plot_grid(plot_A_1,plot_A_2, nrow = 2, ncol = 1, labels = c("",""), rel_widths = c(0.5))

########################################
############ FINAL SUPP FIG ############
########################################

# SuppFig_final <- plot_grid(fig_up, fig_bottom, nrow = 2, ncol = 1, rel_heights = c(0.4,0.6))
SuppFig_final <- plot_grid(fig_up, nrow = 1, ncol = 1, rel_heights = c(0.5))

final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/SuppFig/SuppFig25_complete.png"
ggsave(final_figure_path, SuppFig_final, width = 175, height = 150, units = "mm")
final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/SuppFig/SuppFig25_complete.pdf"
ggsave(final_figure_path, SuppFig_final, width = 175, height = 150, units = "mm")

