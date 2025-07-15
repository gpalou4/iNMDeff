library("tidyverse")
library("ggpubr")
library("ggplot2")
library("dplyr")

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

# 1) Data
endogenous_NMD_genesets <-  c("endogenous_NMD_Colombo","endogenous_NMD_Karousis","endogenous_NMD_Tani","endogenous_NMD_Courtney","endogenous_NMD_ensembl",
                      "endogenous_NMD_all","endogenous_NMD_Consensus","endogenous_SMG6","endogenous_SMG7",
                      "endogenous_non_NMD_neg_control","endogenous_non_NMD_neg_control_with_NMD_features")
ASE_NMD_genesets <- c("ASE_PTC_NMD_triggering_0.01","ASE_PTC_NMD_evading_0.01","ASE_synonymous_0.01",
                      "ASE_PTC_NMD_triggering_0.2","ASE_PTC_NMD_evading_0.2","ASE_synonymous_0.2")
# 1.1) sample NMD efficiencies TCGA
# PTC // ASE // Endogenous
sample_NMD_efficiencies_TCGA_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/cancers/pancancer/NMD_efficiencies_TCGA.txt"
sample_NMD_efficiencies_TCGA <- read.table(file = sample_NMD_efficiencies_TCGA_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
sample_NMD_efficiencies_TCGA <- modify_NMDeff_dataframe(sample_NMD_efficiencies_TCGA, dataset = "TCGA", scale = TRUE)

output_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/MSI"

# NMDeff Samples
for (NMD_method in c("ASE_0.2","endogenous")) {

    print(NMD_method)
    if (NMD_method == "ASE_0.2") {
        NMDeff <- "ASE_PTC_NMD_triggering_0.2"
        ylim <- c(-3,3)
        suppfig <- "C"
    } else if (NMD_method == "PTCs") {
        NMDeff <- "PTCs_stopgain"
    } else if (NMD_method == "endogenous") {
        NMDeff <- "endogenous_NMD_Consensus"
        ylim <- c(-3,3)
        suppfig <- "D"
    }

    sample_NMD_efficiencies_TCGA_filt <- sample_NMD_efficiencies_TCGA %>% filter(cancer_type %in% c("STAD","UCEC","COAD"))
    sample_NMD_efficiencies_TCGA_filt$MSI_status_2 <- "MSS"
    sample_NMD_efficiencies_TCGA_filt[grep("MSI",sample_NMD_efficiencies_TCGA_filt$cancer_type_strat),"MSI_status_2"] <- "MSI"
    sample_NMD_efficiencies_TCGA_filt[grep("POLE",sample_NMD_efficiencies_TCGA_filt$cancer_type_strat),"MSI_status_2"] <- "POLE"
    table(sample_NMD_efficiencies_TCGA_filt$MSI_status_2)
    sample_NMD_efficiencies_TCGA_filt <- sample_NMD_efficiencies_TCGA_filt %>%
                    filter(MSI_status_2 %in% c("MSI","MSS"))

    # Save
    write.table(sample_NMD_efficiencies_TCGA_filt, file = paste0("/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig5/SuppFig5",suppfig,".txt"), 
                    sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
    saveRDS(sample_NMD_efficiencies_TCGA_filt, paste0("/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig5/SuppFig5",suppfig,".RData"))

    p <- ggboxplot(sample_NMD_efficiencies_TCGA_filt, x = "MSI_status_2", y = NMDeff,
            color = "MSI_status_2", palette = "jco",
            add = "jitter",
            facet.by = "cancer_type", short.panel.labs = TRUE)
    p <- p + theme_bw() + xlab("") +
                # facet_wrap(.~cancer_type) +
                ylab("NMD efficiency") + ggtitle(gsub("0.01|0.2","",gsub("_"," ",gsub("endogenous_|ASE_","",NMDeff)))) +
                theme(plot.title = element_text(hjust = 0.5, size = 26),
                    axis.title.x = element_text(color="black", size=26, face="bold"),
                    axis.title.y = element_text(color="black", size=26, face="bold"),
                    axis.text.x = element_text(color="black", size=24),
                    axis.text.y = element_text(color="black", size=24), 
                    axis.text=element_text(size=18),
                    legend.title = element_blank(),
                    legend.key.size = unit(1, 'cm'),
                    legend.text = element_text(size=22),
                    legend.position = "top") #+ coord_cartesian(ylim = ylim)# +
    # Use only p.format as label. Remove method name.
    plot <- p + stat_compare_means(label = "p.format", method = "wilcox.test", hide.ns = TRUE, method.args = list(alternative =
          "greater"))

    png(paste0(output_path,"/MSI_vs_MSS_iNMDeff_",NMD_method,".png"), width = 4000, height = 3500, res = 300)
    print(plot)
    dev.off()

}
