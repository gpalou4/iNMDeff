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
############ FIG 5A ############
################################

# Data
figureA_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig5/Fig5A.png"
# Read the PNG image using the png package
figureA_png <- readPNG(figureA_path)
# Convert the PNG image to a raster object
raster_grob <- rasterGrob(figureA_png, interpolate=TRUE)
# Create a ggplot object with the raster image
plot_A <- ggplot() + 
  annotation_custom(raster_grob, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  theme_void() +  # Remove axes and labels
  theme_classic()

################################
############ FIG 5B ############
################################

# Data
input_figureB_C <- readRDS("/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig5/Fig5B_C.RData")

# Plot
tmp_df <- input_figureB_C[,c("variant","FDR_threshold_used","database_discovery","database_validation","randomization","matched_tissue")]
tmp_df <- tmp_df[!duplicated(tmp_df),]
plot_B <- tmp_df %>%
        filter(database_discovery == "GTEx" & database_validation == "TCGA" & randomization == "no") %>%
        group_by(matched_tissue,FDR_threshold_used) %>%
        summarise(replicated_hits = n()) %>%
        arrange(replicated_hits) %>%
            ggplot(aes(x = matched_tissue, y = replicated_hits, fill = factor(FDR_threshold_used))) + 
                geom_bar( stat = 'identity', position = position_dodge() ) +
                coord_flip() +
                labs(title = "", x = "", y = "Replicated hits", fill = "% FDR threshold") +                 
                ggplot_theme() + scale_fill_brewer(palette = "OrRd", direction = -1) +
                theme(legend.position = "bottom", 
                        plot.margin = unit(c(2.5, 0.7, 2.5, 0), "cm"),
                        legend.margin = margin(t = 0, r = 0, b = 0, l = 0))

################################
############ FIG 5C ############
################################

# df <- input_figureB_C %>%
#         # filter(database_discovery %in% "GTEx" & database_validation %in% "TCGA") %>%
#         filter(database_discovery %in% "TCGA" & database_validation %in% "GTEx") %>%
#         dplyr::select(variant,dataset,NMD_method,database,tissue)

# Change stuff
input_figureB_C$gene_tissue <- paste0(input_figureB_C$variant," - ",input_figureB_C$matched_tissue)
input_figureB_C$database_translation <- paste0(input_figureB_C$database_discovery,"/",input_figureB_C$database_validation)
input_figureB_C$disc_val <- ifelse(input_figureB_C$disc_val == "discovery","disc","val")
input_figureB_C <- input_figureB_C %>%
            filter(FDR_threshold_used == 2 & randomization == "no") %>%
            mutate(facet_order = factor(paste0(database, " - ",disc_val),
                            levels = c("GTEx - disc", "TCGA - val",
                                        "TCGA - disc", "GTEx - val")))
input_figureB_C <- input_figureB_C %>%
                mutate(n_carriers = ifelse(variant == "LAMC1" & NMD_method == "Endogenous" & dataset == "PTV_Missense_CADD15_0.1perc", 4, n_carriers)) %>%
                mutate(n_carriers = ifelse(variant == "CRTC1" & NMD_method == "Endogenous" & dataset == "PTV_Missense_CADD15_0.1perc", 3, n_carriers))
# Dataset
input_figureB_C$dataset <- gsub("PTV","pLoF",input_figureB_C$dataset)
input_figureB_C$dataset <- gsub("_0.1perc","",input_figureB_C$dataset)

# Order
order_vector <- input_figureB_C[order(input_figureB_C$matched_tissue),"gene_tissue"]
input_figureB_C$gene_tissue <- factor(input_figureB_C$gene_tissue, levels = unique(order_vector))

input_figureB_C$facet_order <- gsub("GTEx","GTex",input_figureB_C$facet_order)

plot_C <- input_figureB_C %>%
        mutate(NMD_method = if_else( NMD_method == "Endogenous","ETG","ASE")) %>%
        filter(database_discovery == "GTEx" & database_validation == "TCGA") %>%
        group_by(disc_val) %>%
            ggplot(aes(x = dataset,y = gene_tissue)) +
                geom_tile(aes(fill = n_carriers)) +
                geom_text(aes(label = n_carriers),color = "black") +
                facet_nested(. ~ NMD_method + facet_order) +
                ggplot_theme_bw() + 
                theme(axis.text.x = element_text(angle = 45, hjust = 1),
                        legend.position = "top") +
                scale_fill_gradientn(colours = c('grey','#F21A00')) +
                labs(title = "Replicated hits at a 2% FDR", x = "", y = "", fill = '# of individuals')

###############################################
############ PANEL UP --> FIG 5A-B ############
###############################################

plot_A_B <- plot_grid(plot_A, plot_B, nrow = 1, ncol = 2, labels = c("A", "B"), 
                label_size = 20, rel_heights = c(0.7,0.3))
fig_up <- plot_grid( 
        plot_A_B,
        plot_C,
        nrow = 1, labels = c("", "C"), label_size = 20, 
        rel_widths = c(0.5,0.5))

# final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/test.png"
# ggsave(final_figure_path, fig_up, width = 250, height = 200, units = "mm")

################################
############ FIG 5D ############
################################

# # Data
# input_figureC <- readRDS("/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig5/Fig5C.RData")
# # Remove GTEx non-European samples and missing PCs samples
# outlier_path <- paste0("/g/strcombio/fsupek_cancer1/gpalou/GTEx/samples_metadata/ancestry_and_admixture_and_PCA_outliers.txt")
# all_outlier_samples <- read.table(file = outlier_path)$V1
# input_figureC <- input_figureC[!input_figureC$sample_short %in% all_outlier_samples,]
# input_figureC$NMD_method <- ifelse(input_figureC$NMD_method == "NMD Consensus", "ETG","ASE")

# # Selected genes
# #c("NUP153","KDM6B","TRAP1","CRTC1","PXDN","PDIA2","LAMC1","Fig5","COL19A1","CCDC114")
# # Fig5 amd CCDC114 are good only on ASE
# selected_genes <- c("KDM6B","PXDN")
# input_figureC_filt <- input_figureC %>%
#                         filter(gene %in% selected_genes)
# table(input_figureC_filt$gene)

# # Add adjusted P-values from the SKAT-O
# output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/rare_germline_variants/replication/nogenelist/between/RGWAS_replication_hits.txt")
# RGWAS_replication_hits_final <- read.table(file = output_path, header = TRUE, sep = "\t")
# RGWAS_replication_selected_hits <- RGWAS_replication_hits_final[RGWAS_replication_hits_final$variant %in% selected_genes & 
#                         RGWAS_replication_hits_final$FDR_threshold_used == 2 &
#                         RGWAS_replication_hits_final$dataset == "PTV_Missense_CADD15_0.1perc" &
#                         RGWAS_replication_hits_final$randomization == "no",]
# RGWAS_replication_selected_hits <- RGWAS_replication_selected_hits[,c("variant","database_validation","database_discovery","NMD_method","tissue","pvalue_FDR_adjusted")]
# colnames(RGWAS_replication_selected_hits) <- c("gene","database_validation","database_discovery","NMD_method","tissue","pvalue_FDR_adjusted")
# RGWAS_replication_selected_hits$NMD_method <- ifelse(RGWAS_replication_selected_hits$NMD_method == "Endogenous","ETG","ASE")
# RGWAS_replication_selected_hits <- RGWAS_replication_selected_hits[order(RGWAS_replication_selected_hits$gene),]
# # input_figureC_filt <- merge(input_figureC_filt,RGWAS_replication_selected_hits, all.x = TRUE)

# input_figureC_filt$gene <- factor(ifelse(input_figureC_filt$gene == "KDM6B","KDM6B (Thyroid-THCA)","PXDN (LGG-Brain_Substantia_Nigra)"))

# plot_C <- ggplot(data = input_figureC_filt, aes(x = germ_mut, y = NMD_efficiency, fill = factor(germ_mut, levels = c("WT","MUT")))) +
#         geom_violin(position = position_dodge(width = 0.9)) +
#         geom_jitter(aes(fill = germ_mut), 
#                     position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0, dodge.width = 0.9),
#                     alpha = 0.25, size = 1) +
#         geom_boxplot(width = 0.2, position = position_dodge(width = 0.9), color="black", alpha=0.2) +
#         coord_cartesian(ylim = c(-2,3)) +
#         facet_nested(dataset ~ gene + NMD_method ) +
#         geom_hline(yintercept = 0, linetype="dashed", color = "red", size=0.5) +
#         labs(y = "iNMDeff", x = "", fill = "Germline variant") +
#         ggplot_theme_bw() + scale_fill_brewer(palette = "Dark2") +
#         theme(legend.position='top') + #strip.background =  element_rect(colour="grey", fill="grey")
#         stat_compare_means(data = subset(input_figureC_filt, !(gene == "PXDN" & NMD_method == "ASE" 
#                                 & dataset == "GTEx")),
#                         aes(group=germ_mut), size = 5,
#                             label.y = c(2.4),
#                             label.x = 1, na.rm = TRUE,
#                             label = "p.format", method = "wilcox.test", hide.ns = TRUE)
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

# final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/test.png"
# ggsave(final_figure_path, plot_C, width = 250, height = 200, units = "mm")

#################################################
############ PANEL MIDDLE --> FIG 5C ############
#################################################

# fig_middle <- plot_grid(plot_C, nrow = 1, ncol = 1,  
#                        labels = c("C"), label_size = 20)

################################
############ FIG 5D ############
################################

# Data
input_figureD_1 <- readRDS("/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig5/Fig5D_1.RData")
input_figureD_1 <- input_figureD_1[,!colnames(input_figureD_1) %in% NA]
input_figureD_2 <- readRDS("/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig5/Fig5D_2.RData")
input_figureD_2 <- input_figureD_2[,!colnames(input_figureD_2) %in% NA]

plot_D_1 <- input_figureD_1 %>%
                filter(gene %in% "KDM6B") %>%
                filter(!tissue %in% c("ACC","LAML","HNSC","OV","LUSC")) %>%
                ggplot(aes(y = tissue, x = coefficient, color = factor(NMD_method))) +
                geom_point(size = 2, shape = 16) +
                facet_nested( ~ gene + NMD_method, scales = "free_x") +
                geom_errorbarh(aes(xmin = CI_2.5_values, xmax = CI_97.5_values), size = 1, height = 0.2) +
                geom_vline(xintercept = 0, color = "red", linetype = "dashed", cex = 0.8, alpha = 0.5) +
                labs(x = "Association effect size", y = "", fill = "", title = "TCGA") + scale_color_brewer(palette = "Paired") +
                ggplot_theme_bw() +
                theme(legend.position = "none") +
                geom_text(aes(label = ifelse(SKATO_p_value_FDR_adjusted < 0.05, "*", "")), 
                                        #position = position_dodge(width = 1), 
                                size = 8, hjust = 0.5, color = "black", nudge_y = 0.1) +
                geom_text(aes(label = ifelse(SKATO_p_value_FDR_adjusted > 0.05 & SKATO_p_value_FDR_adjusted < 0.2, "*", "")), 
                                        #position = position_dodge(width = 1), 
                                size = 8, hjust = 0.5, color = "red", nudge_y = 0.1)

plot_D_2 <- input_figureD_2 %>%
        filter(gene == "KDM6B") %>%
          ggplot(aes(y = tissue, x = coefficient, color = factor(NMD_method))) +
                geom_point(size = 2, shape = 16) +
                facet_nested( ~ gene + NMD_method, scales = "free_x") +
                geom_errorbarh(aes(xmin = CI_2.5_values, xmax = CI_97.5_values), size = 1, height = 0.2) +
                geom_vline(xintercept = 0, color = "red", linetype = "dashed", cex = 0.8, alpha = 0.5) +
                labs(x = "Association effect size", y = "", fill = "", title = "GTex") + scale_color_brewer(palette = "Paired") +
                ggplot_theme_bw() +
                theme(legend.position = "none") +
                geom_text(aes(label = ifelse(SKATO_p_value_FDR_adjusted < 0.05, "*", "")), 
                                #     position = position_dodge(width = 1),
                                size = 8, hjust = 0.5, color = "black", , nudge_y = 0.1) +
                geom_text(aes(label = ifelse(SKATO_p_value_FDR_adjusted > 0.05 & SKATO_p_value_FDR_adjusted < 0.2, "*", "")), 
                                        #position = position_dodge(width = 1), 
                                size = 8, hjust = 0.5, color = "red", nudge_y = 0.1)

plot_D <- plot_grid(plot_D_1, plot_D_2, nrow = 1, ncol = 2, rel_widths = c(0.4,0.6)) +
        theme(plot.margin = unit(c(0, 2, 0, 2), "cm"))

#################################################
############ PANEL BOTTOM --> FIG 5D ############
#################################################

fig_bottom <- plot_grid(plot_D, nrow = 1, ncol = 1,
                       labels = c("D"), label_size = 20)

#############################################
############ FINAL FIG --> FIG 5 ############
#############################################

Fig5_final <- plot_grid(fig_up, fig_bottom, nrow = 2, rel_heights = c(0.4,0.6))

final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/Fig5_complete.pdf"
ggsave(final_figure_path, Fig5_final, width = 300, height = 325, units = "mm", device = "pdf")
final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/Fig5_complete.png"
ggsave(final_figure_path, Fig5_final, width = 300, height = 325, units = "mm", device = "png")
final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/Fig5_complete.tiff"
ggsave(final_figure_path, Fig5_final, width = 300, height = 325, units = "mm", device = "tiff")





