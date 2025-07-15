modify_NMDeff_dataframe <- function(sample_NMDeff, dataset, scale = FALSE) {
  # Convert some columns to factors
  if (dataset == "TCGA") {
    factor_cols <- c("cancer_type","cancer_type_MSI","cancer_type_strat","cancer_subtype","LF_remove","purity_remove", "MSI_status",
                    "batch_portion","batch_plate","batch_center","batch_vial","TCGA_full_barcode")
    # Remove "TCGA" from the cancer type
    sample_NMDeff$cancer_type <- gsub("TCGA-","",sample_NMDeff$cancer_type)
    sample_NMDeff$cancer_type_strat <- gsub("TCGA-","",sample_NMDeff$cancer_type_strat)
    sample_NMDeff$cancer_type_MSI <- gsub("TCGA-","",sample_NMDeff$cancer_type_MSI)
  } else if (dataset == "GTEx") {
    factor_cols <- c("tissue","sample")
  }
  sample_NMDeff[factor_cols] <- lapply(sample_NMDeff[factor_cols], factor) 
  # Rename NMD genesets for the 3 methods
  all_NMD_genesets <- c(endogenous_NMD_genesets,ASE_NMD_genesets,"NMDeff_mean")
  # Endogenous
  colnames(sample_NMDeff)[colnames(sample_NMDeff) %in% c("endogenous_NMD_global_2_shared","endogenous_NMD_global_2_shared_randomized")] <- c("endogenous_NMD_Consensus","endogenous_NMD_Consensus_randomized")
  colnames(sample_NMDeff)[colnames(sample_NMDeff) %in% c("endogenous_NMD_global","endogenous_NMD_global_randomized")] <- c("endogenous_NMD_all","endogenous_NMD_all_randomized")
  # ASE
  colnames(sample_NMDeff)[colnames(sample_NMDeff) %in% c("ASE_stopgain_0.2","ASE_stopgain_0.2_randomized")] <- c("ASE_PTC_NMD_triggering_0.2","ASE_PTC_NMD_triggering_0.2_randomized")
  colnames(sample_NMDeff)[colnames(sample_NMDeff) %in% c("ASE_stopgain_0.01","ASE_stopgain_0.01_randomized")] <- c("ASE_PTC_NMD_triggering_0.01","ASE_PTC_NMD_triggering_0.01_randomized")
  colnames(sample_NMDeff)[colnames(sample_NMDeff) %in% c("ASE_stopgain_NMD_evading_0.2","ASE_stopgain_NMD_evading_0.2_randomized")] <- c("ASE_PTC_NMD_evading_0.2","ASE_PTC_NMD_evading_0.2_randomized")
  colnames(sample_NMDeff)[colnames(sample_NMDeff) %in% c("ASE_stopgain_NMD_evading_0.01","ASE_stopgain_NMD_evading_0.01_randomized")] <- c("ASE_PTC_NMD_evading_0.01","ASE_PTC_NMD_evading_0.01_randomized")

  # Scale NMD genesets for the three methods
  # Change the sign (coefficients are reversed), so higher values means high NMDeff
  sample_NMDeff[,all_NMD_genesets] <- -sample_NMDeff[,all_NMD_genesets]
  # Scale
  if (isTRUE(scale)) {
      sample_NMDeff[,all_NMD_genesets] <- scale(sample_NMDeff[,all_NMD_genesets])
  }
  # Filter samples with low PTC number in ASE
  sample_NMDeff[which(sample_NMDeff$ASE_num_PTCs_0.2 < 3),c("ASE_PTC_NMD_triggering_0.2","ASE_PTC_NMD_evading_0.2","NMDeff_mean")] <- NA
  sample_NMDeff[which(sample_NMDeff$ASE_num_PTCs_0.01 < 3),c("ASE_PTC_NMD_triggering_0.01","ASE_PTC_NMD_evading_0.01")] <- NA

  return(sample_NMDeff)
}

library(dplyr)
library(RColorBrewer)
library(ggplot2)

# 1.1) sample NMD efficiencies TCGA
endogenous_NMD_genesets <-  c("endogenous_NMD_Colombo","endogenous_NMD_Karousis","endogenous_NMD_Tani","endogenous_NMD_Courtney","endogenous_NMD_ensembl",
                      "endogenous_NMD_all","endogenous_NMD_Consensus","endogenous_SMG6","endogenous_SMG7",
                      "endogenous_non_NMD_neg_control","endogenous_non_NMD_neg_control_with_NMD_features")
ASE_NMD_genesets <- c("ASE_PTC_NMD_triggering_0.01","ASE_PTC_NMD_evading_0.01","ASE_synonymous_0.01",
                      "ASE_PTC_NMD_triggering_0.2","ASE_PTC_NMD_evading_0.2","ASE_synonymous_0.2")

# PTC // ASE // Endogenous
sample_NMD_efficiencies_TCGA_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/cancers/pancancer/NMD_efficiencies_TCGA.txt"
sample_NMD_efficiencies_TCGA <- read.table(file = sample_NMD_efficiencies_TCGA_path, header = TRUE, sep = "\t")
sample_NMD_efficiencies_TCGA <- modify_NMDeff_dataframe(sample_NMD_efficiencies_TCGA, dataset = "TCGA", scale = TRUE)

# 1.2) sample NMD efficiencies GTEx
sample_NMD_efficiencies_GTEx_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/GTEx/pantissue/NMD_efficiencies_GTEx.txt"
sample_NMD_efficiencies_GTEx <- read.table(file = sample_NMD_efficiencies_GTEx_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
sample_NMD_efficiencies_GTEx <- modify_NMDeff_dataframe(sample_NMD_efficiencies_GTEx, dataset = "GTEx", scale = TRUE)

# 2) Boxplot of NMDeff by bins, for Endogenous and ASE

NMDeff_by_bins_END_vs_ASE <- function(dataset) {

    if (dataset == "TCGA") {
        sample_NMD_efficiencies <- sample_NMD_efficiencies_TCGA  
    } else if (dataset == "GTEx") {
        sample_NMD_efficiencies <- sample_NMD_efficiencies_GTEx
    }

    sample_NMD_efficiencies <- sample_NMD_efficiencies %>%
        mutate(
            END_bins = cut(
                endogenous_NMD_Consensus, 
                breaks = quantile(endogenous_NMD_Consensus, probs = seq(0, 1, by = 0.1), na.rm = TRUE),
                labels = 1:10, 
                include.lowest = TRUE
            )
        ) %>%
        mutate(
            ASE_bins = cut(
            ASE_PTC_NMD_triggering_0.2, 
            breaks = quantile(ASE_PTC_NMD_triggering_0.2, probs = seq(0, 1, by = 0.1), na.rm = TRUE),
            labels = 1:10, 
            include.lowest = TRUE
            )
        )
    sample_NMD_efficiencies_filt <- sample_NMD_efficiencies[,c("endogenous_NMD_Consensus","ASE_PTC_NMD_triggering_0.2","END_bins","ASE_bins")]

    sample_NMD_efficiencies_stack <- stack(sample_NMD_efficiencies_filt)
    sample_NMD_efficiencies_stack$bins <- rep(sample_NMD_efficiencies_filt$ASE_bins,2)
    colnames(sample_NMD_efficiencies_stack) <- c("NMDeff","NMD_method","bins")
    sample_NMD_efficiencies_stack <- na.omit(sample_NMD_efficiencies_stack)
    
    plot_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/correlation_methods/sample_level/",dataset,"_NMDeff_bins_END_vs_ASE.png")
    p <- ggplot(data = sample_NMD_efficiencies_stack, aes(x = bins, y = NMDeff, fill = factor(NMD_method))) +
        #geom_violin(draw_quantiles = TRUE, na.rm = TRUE) + #coord_flip(ylim = ylim) +
        geom_boxplot(color="black", alpha=0.75) + ylim(c(-3,3)) +
        #geom_point(alpha=0.25)
        xlab("Percentiles (%), ordered by ASE method") + ylab("Individual NMD efficiency") +
        scale_x_discrete(labels = paste0(1:10,"0")) +
        geom_hline(yintercept = 0, linetype="dashed", color = "red", size=0.5) +
        theme_bw(base_size = 25) + 
        theme(plot.title = element_text(hjust = 0.5, size = 40),
                axis.title.x = element_text(color="black", size=45),
                axis.title.y = element_text(color="black", size=45),
                axis.text.x = element_text(color="black", size=40),
                axis.text.y = element_text(color="black", size=40),
                legend.position = "top",
                legend.text = element_text(color="black", size=40),
                legend.title = element_text(color="black", size=45)) +
        scale_fill_brewer(palette = "Accent", labels = c("ETG","ASE"), direction = -1) +
        guides(fill = guide_legend(override.aes = list(size = 15), title = "NMD methods"))
        # annotate("text",
        #         x = 1:length(table(sample_NMD_efficiencies_TCGA_stack$bins)),
        #         y = aggregate( NMDeff ~ bins, sample_NMD_efficiencies_TCGA_stack, median)[ , 2],
        #         label = table(sample_NMD_efficiencies_TCGA_stack$bins),
        #         col = "black",
        #         hjust = 1,
        #         vjust = -1,
        #         size = 10)
    png(plot_path, width = 4500, height = 4500, res = 300)
    print(p)
    dev.off()
    return(sample_NMD_efficiencies_stack)
}

TCGA_sample_NMD_efficiencies_stack <- NMDeff_by_bins_END_vs_ASE(dataset = "TCGA")
write.table(TCGA_sample_NMD_efficiencies_stack, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig3/SuppFig3E.txt", 
            sep = "\t", quote = TRUE, col.names = TRUE, row.names = FALSE)
# save R object
saveRDS(TCGA_sample_NMD_efficiencies_stack, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig3/SuppFig3E.RData")
GTEx_sample_NMD_efficiencies_stack <- NMDeff_by_bins_END_vs_ASE(dataset = "GTEx")
write.table(GTEx_sample_NMD_efficiencies_stack, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig3/SuppFig3F.txt", 
            sep = "\t", quote = TRUE, col.names = TRUE, row.names = FALSE)
# save R object
saveRDS(GTEx_sample_NMD_efficiencies_stack, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig3/SuppFig3F.RData")

# TCGA
df1 <- TCGA_sample_NMD_efficiencies_stack[TCGA_sample_NMD_efficiencies_stack$bins == 1,]
aggregate(NMDeff ~ NMD_method, data = df1, median)
df2 <- TCGA_sample_NMD_efficiencies_stack[TCGA_sample_NMD_efficiencies_stack$bins == 10,]
aggregate(NMDeff ~ NMD_method, data = df2, median)

df <- TCGA_sample_NMD_efficiencies_stack %>%
            filter(NMD_method == "endogenous_NMD_Consensus" & bins %in% c(1,10))

"Wilcoxon rank sum test is equivalent to the Mann-Whitney test, when pairing = FALSE"
wilcox.test(NMDeff ~ bins, data=df, paired = TRUE )

#GTEx

df1 <- GTEx_sample_NMD_efficiencies_stack[GTEx_sample_NMD_efficiencies_stack$bins == 1,]
aggregate(NMDeff ~ NMD_method, data = df1, median)
df2 <- GTEx_sample_NMD_efficiencies_stack[GTEx_sample_NMD_efficiencies_stack$bins == 10,]
aggregate(NMDeff ~ NMD_method, data = df2, median)

df <- GTEx_sample_NMD_efficiencies_stack %>%
            filter(NMD_method == "endogenous_NMD_Consensus" & bins %in% c(1,10))

"Wilcoxon rank sum test is equivalent to the Mann-Whitney test, when pairing = FALSE"
wilcox.test(NMDeff ~ bins, data=df, paired = FALSE )
