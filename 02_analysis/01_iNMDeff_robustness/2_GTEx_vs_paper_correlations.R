library("ggplot2")
library("ggrepel")
library("ggpubr")
library("ggpmisc")

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

# 1.1) sample NMD efficiencies GTEx
endogenous_NMD_genesets <-  c("endogenous_NMD_Colombo","endogenous_NMD_Karousis","endogenous_NMD_Tani","endogenous_NMD_Courtney","endogenous_NMD_ensembl",
                      "endogenous_NMD_all","endogenous_NMD_Consensus","endogenous_SMG6","endogenous_SMG7",
                      "endogenous_non_NMD_neg_control","endogenous_non_NMD_neg_control_with_NMD_features")
ASE_NMD_genesets <- c("ASE_PTC_NMD_triggering_0.01","ASE_PTC_NMD_evading_0.01","ASE_synonymous_0.01",
                      "ASE_PTC_NMD_triggering_0.2","ASE_PTC_NMD_evading_0.2","ASE_synonymous_0.2")
NMD_genesets <- c("endogenous_NMD_Consensus","ASE_PTC_NMD_triggering_0.01","ASE_PTC_NMD_triggering_0.2")
sample_NMD_efficiencies_GTEx_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/GTEx/pantissue/NMD_efficiencies_GTEx.txt"
sample_NMD_efficiencies_GTEx <- read.table(file = sample_NMD_efficiencies_GTEx_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
sample_NMD_efficiencies_GTEx <- modify_NMDeff_dataframe(sample_NMD_efficiencies_GTEx, dataset = "GTEx", scale = TRUE)

# # Filter
# sample_NMD_efficiencies_GTEx[which(sample_NMD_efficiencies_GTEx$ASE_num_PTCs_0.2 < 2),c("ASE_stopgain_0.2","NMDeff_mean")] <- NA
# sample_NMD_efficiencies_GTEx[which(sample_NMD_efficiencies_GTEx$ASE_num_PTCs_0.01 < 2),c("ASE_stopgain_0.01")] <- NA
# # Scale
# sample_NMD_efficiencies_GTEx[,NMD_genesets] <- scale(sample_NMD_efficiencies_GTEx[,NMD_genesets])

# 1.2) Paper GTEx ranking from ASE method
GTEx_ASE_paper <- data.frame(tissues=c("BRNCHB","BRNCHA","SPLEEN","NERVET","PRSTTE","BRNSNG","OVARY","BRNCTXA","UTERUS","BRNHPP","CLNSGM","SNTTRM","SKINS","THYROID","ESPMSL","TESTIS","BRNAMY","BRNPTM","LUNG","SKINNS","BREAST","ESPGEJ","BRNNCC","PTTARY","ADPSBQ","ARTTBL","BRNCDT","BRNCTXB","HRTAA","SLVRYG","STMACH","VAGINA","WHLBLD","PNCREAS","ARTAORT","BRNHPT","HRTLV","CLNTRN","ESPMCS","KDNCTX","LCL","BRNACC","ADPVSC","ADRNLG","ARTCRN","LIVER","BRNSPC","MSCLSK","FIBRBLS"))
GTEx_ASE_paper$GTEx_ASE_paper_ranking <- rev(1:length(GTEx_ASE_paper$tissues))

# 1.3) Correlation for each NMD method
corr_output_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/correlation_methods"
formula <- y ~ x
ASE_VAF <- 0.2

GTEx_ranking_all <- list()

plots_list <- lapply(1:2, function(n) {

    if (n == 1) {
        NMD_method <- paste0("ASE_PTC_NMD_triggering_",ASE_VAF)
        NMD_method_char <- "ASE"
    } else if (n == 2) {
        NMD_method <- paste0("endogenous_NMD_Consensus")
        NMD_method_char <- "ETG"
    } else if (n == 3) {
        NMD_method <- paste0("PTCs_stopgain_NMD_triggering")
    } else if (n == 4) {
        NMD_method <- paste0("NMDeff_PCA")
    } else if (n == 5) {
        NMD_method <- paste0("NMDeff_mean")
    }

    # Ranking NMDeff from our NMD method
    GTEx_NMD_method <- aggregate(eval(parse(text = NMD_method)) ~ acronyms, data = sample_NMD_efficiencies_GTEx, median)
    GTEx_NMD_method_sort <- GTEx_NMD_method[order(GTEx_NMD_method[,2], decreasing = TRUE),]
    GTEx_NMD_method_sort[,paste0("GTEx_",NMD_method_char,"_iNMDeff_method_ranking")] <- 1:length(GTEx_NMD_method_sort$acronyms)
    # Merge
    GTEx_ranking <- merge(GTEx_NMD_method_sort,GTEx_ASE_paper, by.x = "acronyms", by.y = "tissues", all = TRUE)
    GTEx_ranking <- na.omit(GTEx_ranking)
    colnames(GTEx_ranking)[2] <- paste0(NMD_method_char,"_iNMDeff")
    # Save
    GTEx_ranking_all[[NMD_method_char]] <<- GTEx_ranking
    GTEx_ranking[,2] <- NULL

    # Plot
    p <- ggplot(data = GTEx_ranking, mapping = aes(x = eval(parse(text = paste0("GTEx_",NMD_method_char,"_iNMDeff_method_ranking"))), y = GTEx_ASE_paper_ranking)) +
        geom_point(alpha = 0.5) +
        geom_label_repel(aes(label=acronyms, color = "black"), size=3, nudge_y=0.05, max.overlaps = nrow(GTEx_ranking)) +
        geom_smooth(method = "lm", formula = formula, se = FALSE, size = 1) +
        #scale_color_brewer(palette = "Pastel1") + 
        # stat_fit_glance(method = 'cor.test',
        #                 method.args = list(x = GTEx_ranking$GTEx_NMD_method_ranking, y = GTEx_ranking$GTEx_ASE_paper_ranking, method = "spearman"),
        #                 geom = 'text',
        #                 aes(label = sprintf('r[s]~"="~%.2f~~italic(P)~"="~%.10f',
        #                     stat(estimate), stat(p.value))),
        #                 label.x = 10, label.y = c(40,41), size = 5, parse = TRUE, vstep = 5, hstep = 0.5, na.rm = FALSE) +
        xlab(paste0(gsub("0.01|0.2","",gsub("_"," ",gsub("endogenous_|ASE_","",NMD_method)))," tissue ranking")) + 
        ylab(paste0("Teral et al. tissue ranking")) +
        theme_classic() +
        theme(plot.title = element_text(hjust = 0.5, size = 20),
            axis.title.x = element_text(color="black", size=20, face="bold"),
            axis.title.y = element_text(color="black", size=20, face="bold"),
            legend.position='none')
    p + stat_cor(p.accuracy = 0.0000000001, r.accuracy = 0.01, size = 5)
})
png(paste0(corr_output_path,"/tissue_ranking/Teral_et_al/GTEx_vs_paper_correlation.png"), width = 6500, height = 2000, res = 300)
p <- cowplot::plot_grid(plotlist=plots_list, labels = "AUTO", align = "v", ncol = 2, nrow = 1)
print(p)
dev.off()

GTEx_ranking_iNMDeff <- do.call(cbind,GTEx_ranking_all)
GTEx_ranking_iNMDeff <- GTEx_ranking_iNMDeff[,c(1,2,3,4,6,7)]
colnames(GTEx_ranking_iNMDeff) <- c("GTEx_tissues","ASE_iNMDeff","GTEx_ASE_iNMDeff_method_ranking",
                                    "GTEx_ASE_Teral_et_al_ranking","ETG_iNMDeff","GTEx_ETG_iNMDeff_method_ranking")
GTEx_ranking_iNMDeff <- GTEx_ranking_iNMDeff[,c(1,2,5,3,6,4)]
# Save
write.table(GTEx_ranking_iNMDeff, file = paste0("/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig4/SuppFig4D.txt"), 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
saveRDS(GTEx_ranking_iNMDeff, paste0("/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig4/SuppFig4D.RData"))

