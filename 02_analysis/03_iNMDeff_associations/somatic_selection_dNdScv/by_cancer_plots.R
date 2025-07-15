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
 
# Libraries
library(dplyr)
#library(tidyverse)
library(ggpubr)
library(dndscv)
library(cowplot)
library(readxl)
library(stringr)

# 1) Data
# 1.1) Ind-NMDeff TCGA
endogenous_NMD_genesets <-  c("endogenous_NMD_Colombo","endogenous_NMD_Karousis","endogenous_NMD_Tani","endogenous_NMD_Courtney","endogenous_NMD_ensembl",
                      "endogenous_NMD_all","endogenous_NMD_Consensus","endogenous_SMG6","endogenous_SMG7",
                      "endogenous_non_NMD_neg_control","endogenous_non_NMD_neg_control_with_NMD_features")
ASE_NMD_genesets <- c("ASE_PTC_NMD_triggering_0.01","ASE_PTC_NMD_evading_0.01","ASE_synonymous_0.01",
                      "ASE_PTC_NMD_triggering_0.2","ASE_PTC_NMD_evading_0.2","ASE_synonymous_0.2")
# TCGA
sample_NMD_efficiencies_TCGA_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/cancers/pancancer/NMD_efficiencies_TCGA.txt"
sample_NMD_efficiencies_TCGA <- read.table(file = sample_NMD_efficiencies_TCGA_path, header = TRUE, sep = "\t")
sample_NMD_efficiencies_TCGA <- modify_NMDeff_dataframe(sample_NMD_efficiencies_TCGA, dataset = "TCGA", scale = FALSE)

# 1.3) Cancer genes (CGC) data
CGC <- read.csv(file = "/g/strcombio/fsupek_cancer1/gpalou/COSMIC/cancer_gene_census_updated.tsv", 
                    header = TRUE, stringsAsFactors = FALSE)
CGC_filt <- CGC[which(!CGC$Role.in.Cancer %in% "non_CGC_NMD"),]
# Gene type
CGC_filt$Role.in.Cancer <- ifelse(CGC_filt$Role.in.Cancer == "oncogene, fusion","oncogene",CGC_filt$Role.in.Cancer)
CGC_filt$Role.in.Cancer <- ifelse(CGC_filt$Role.in.Cancer == "oncogene, TSG","both",CGC_filt$Role.in.Cancer)
CGC_filt[CGC_filt$Gene.Symbol %in% "TP53","Role.in.Cancer"] <- "TSG"
CGC_filt$Role.in.Cancer <- ifelse(CGC_filt$Role.in.Cancer == "oncogene, TSG, fusion","both",CGC_filt$Role.in.Cancer)
CGC_filt$Role.in.Cancer <- ifelse(CGC_filt$Role.in.Cancer == "TSG, fusion","TSG",CGC_filt$Role.in.Cancer)
# Remove oncogenes that are not Dominant 
CGC_filt <- CGC_filt[-which(CGC_filt$Role.in.Cancer == "oncogene" & CGC_filt$Molecular.Genetics != "Dom"),]
# Manual removal of leukemia-specific cancer genes
CGC_filt <- CGC_filt[-grep("leukaemia|T-ALL",CGC_filt$Tumour.Types.Somatic.),]
table(CGC_filt$Role.in.Cancer)

# 1.4) Mutpanning Cancer Genes
Mutpanning <- read.csv(file = "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/dNdScv/MutPanningGeneTumorPairs.csv", 
                    header = TRUE, stringsAsFactors = FALSE)
Mutpanning <- Mutpanning[order(Mutpanning$Q.value..false.discovery.rate.),]

# 1.5) Solimini TSG STOP genes
STOP_genes <- read_excel("/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/dNdScv/Solimini_NL_et_al_2016.xls", 
                sheet = "Table S7", skip = 1)
STOP_genes <- na.omit(STOP_genes[order(STOP_genes[,"Average Log2 Ratio"],decreasing = TRUE),])
STOP_genes_top <- unique(STOP_genes[,c("Gene Symbol(s)")]) %>% pull()
STOP_genes_top_200 <- STOP_genes_top[1:200]
STOP_genes_top_200 <- gsub(" ","",unlist(strsplit(STOP_genes_top_200,",")))
STOP_genes_top_all <- gsub(" ","",unlist(strsplit(STOP_genes_top,",")))

# 1.6) Solimini Oncogenes GO genes
OG_GO <- read_excel("/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/dNdScv/Solimini_NL_et_al_2016.xls", 
                sheet = "Table S11", skip = 1)
OG_GO <- data.frame(OG_GO)
OG_GO <- OG_GO$Gene.Symbol

# 1.7) Selection of cancer genes
# We will use genes from CGC to do the dNdScv test. 
# Then we will filter for interested genes (Mutpanning, Solimin, etc) in the plots
# 1.7.1) Oncogenes: from CGC
oncogenes_CGC <- as.character(na.omit(CGC_filt[CGC_filt$Role.in.Cancer == "oncogene","Gene.Symbol"]))
# 1.7.2) TSG, from CGC
TSG_CGC <- as.character(na.omit(CGC_filt[CGC_filt$Role.in.Cancer == "TSG","Gene.Symbol"]))
# 1.7.3) Final set of cancer genes
OG_TSG <- unique(c(oncogenes_CGC,TSG_CGC))
# 1.8 Selection of cancer genes for the analysis
# 1.8.1) TSG
# Match CGG with STOP 
matched_genes <- c()
for (gene in TSG_CGC) {
  res <- grep(paste0(gene,"$"),STOP_genes_top)
  if (length(res)==0) {
    next
  } else {
    matched_genes <- c(matched_genes,gene)
  }
}
TSG_STOP_CGC <- matched_genes
# TSG Mutpanning genes with an incidence >=3 cancer types
Mutpanning_freq_genes <- data.frame(table(Mutpanning$Gene))
Mutpanning_freq_genes <- Mutpanning_freq_genes[order(Mutpanning_freq_genes$Freq),]
Mutpanning_freq_genes <- Mutpanning_freq_genes[Mutpanning_freq_genes$Freq >= 4,]
colnames(Mutpanning_freq_genes) <- c("Gene","frequency")
# And also q-value < 0.01
Mutpanning_filt <- Mutpanning[Mutpanning$Gene %in% Mutpanning_freq_genes$Gene,]
filter <- Mutpanning_filt[,"Q.value..false.discovery.rate."] < 0.01
mutpanning_freq_genes <- unique(Mutpanning_filt[filter,"Gene"])
# Intersect with TSG CGC
TSG_mutpanning_CGC <- intersect(TSG_CGC,mutpanning_freq_genes)
# Add STOP
TSG_STOP_CGC_and_mutpanning_CGC <- unique(c(TSG_STOP_CGC,TSG_mutpanning_CGC))
# TSG_STOP_CGC_and_mutpanning_CGC[!TSG_STOP_CGC_and_mutpanning_CGC %in% dNdScv_all_res$gene_name]
# 1.8.2) Oncogenes
# OG Mutpanning genes with an incidence >=3 cancer types
Mutpanning_freq_genes <- data.frame(table(Mutpanning$Gene))
Mutpanning_freq_genes <- Mutpanning_freq_genes[order(Mutpanning_freq_genes$Freq),]
Mutpanning_freq_genes <- Mutpanning_freq_genes[Mutpanning_freq_genes$Freq >= 1,]
colnames(Mutpanning_freq_genes) <- c("Gene","frequency")
# And also q-value < 0.01
Mutpanning_filt <- Mutpanning[Mutpanning$Gene %in% Mutpanning_freq_genes$Gene,]
filter <- Mutpanning_filt[,"Q.value..false.discovery.rate."] < 0.01
mutpanning_freq_genes <- unique(Mutpanning_filt[filter,"Gene"])
# Intersect with OG CGC
OG_mutpanning_CGC <- intersect(oncogenes_CGC,mutpanning_freq_genes)
# Intersect GO genes with CGC
OG_GO_CGC <- intersect(OG_GO,oncogenes_CGC)
# Merge
OG_GO_and_mutpanning_CGC <- unique(c(OG_GO_CGC,OG_mutpanning_CGC))
# OG_GO_and_mutpanning_CGC[!OG_GO_and_mutpanning_CGC %in% dNdScv_all_res$gene_name]

# 2) Plots
cancers <- na.omit(unique(sample_NMD_efficiencies_TCGA$cancer_type_strat))
cancer_ylims <- data.frame(cancer = cancers, nonsense_ylim = 200, missense_ylim = 80)
cancer_ylims <- cancer_ylims %>% 
      mutate(missense_ylim = ifelse(cancer == "BRCA_LumA",40,missense_ylim)) %>%
      mutate(missense_ylim = ifelse(cancer == "COAD_MSI",10,missense_ylim)) %>%
      mutate(nonsense_ylim = ifelse(cancer == "COAD_MSI",40,nonsense_ylim)) %>%
      mutate(nonsense_ylim = ifelse(cancer == "ESCA_ac",500,nonsense_ylim)) %>% # Check
      mutate(missense_ylim = ifelse(cancer == "GBM",20,missense_ylim)) %>%
      mutate(nonsense_ylim = ifelse(cancer == "KICH",500,nonsense_ylim)) %>% # Check
      mutate(missense_ylim = ifelse(cancer == "KIRP",10,missense_ylim)) %>%
      mutate(nonsense_ylim = ifelse(cancer == "KIRP",35,nonsense_ylim)) %>%
      mutate(missense_ylim = ifelse(cancer == "LUAD",40,missense_ylim)) %>%
      mutate(nonsense_ylim = ifelse(cancer == "LUAD",90,nonsense_ylim)) %>%
      mutate(nonsense_ylim = ifelse(cancer == "LUSC",70,nonsense_ylim)) %>%
      mutate(missense_ylim = ifelse(cancer == "LUSC",23,missense_ylim)) %>%
      mutate(nonsense_ylim = ifelse(cancer == "OV",15,nonsense_ylim)) %>% # CHECK
      mutate(nonsense_ylim = ifelse(cancer == "PAAD",15,nonsense_ylim)) %>% # Check
      mutate(nonsense_ylim = ifelse(cancer == "SKCM",25,nonsense_ylim)) %>% # Check
      mutate(nonsense_ylim = ifelse(cancer == "STAD",500,nonsense_ylim)) %>%
      mutate(nonsense_ylim = ifelse(cancer == "STAD_MSI",100,nonsense_ylim)) %>%
      mutate(missense_ylim = ifelse(cancer == "STAD_MSI",20,missense_ylim)) %>%
      mutate(missense_ylim = ifelse(cancer == "THCA",10,missense_ylim)) %>%
      mutate(nonsense_ylim = ifelse(cancer == "THCA",25,nonsense_ylim)) %>%
      mutate(missense_ylim = ifelse(cancer == "UCEC",25,missense_ylim)) %>%
      mutate(nonsense_ylim = ifelse(cancer == "UCEC",125,nonsense_ylim)) %>%
      mutate(missense_ylim = ifelse(cancer == "UCEC_MSI",60,missense_ylim)) %>%
      mutate(nonsense_ylim = ifelse(cancer == "UCEC_MSI",125,nonsense_ylim)) %>%
      mutate(nonsense_ylim = ifelse(cancer == "UCEC_POLE",25,nonsense_ylim)) %>%
      mutate(missense_ylim = ifelse(cancer == "UCEC_POLE",7,missense_ylim))

sig_cancers <- c("COAD_MSI","STAD_MSI","UCEC_POLE","GBM","LUAD","LUSC")
all_list_plots <- c()
# for (cancer in cancers) {
for (cancer in sig_cancers) {

  if (cancer %in% c("LAML","MESO")) { next }
  input <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/selection/dNdScv_all_res_",cancer,".txt")
  dNdScv_all_res <- read.table(file = input, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  # Re-classify cancer genes
  dNdScv_all_res$gene_type <- ifelse(dNdScv_all_res$gene_type == "TSG",NA,dNdScv_all_res$gene_type)
  table(dNdScv_all_res$gene_type)
  gene_filter <- TSG_STOP_CGC_and_mutpanning_CGC
  dNdScv_all_res[which(is.na(dNdScv_all_res$gene_type) & dNdScv_all_res$gene_name %in% c(gene_filter)),"gene_type"] <- "TSG"
  table(dNdScv_all_res$gene_type)
  dNdScv_all_res$gene_type <- ifelse(dNdScv_all_res$gene_type == "oncogene",NA,dNdScv_all_res$gene_type)
  table(dNdScv_all_res$gene_type)
  gene_filter <- OG_GO_and_mutpanning_CGC
  dNdScv_all_res[which(is.na(dNdScv_all_res$gene_type) & dNdScv_all_res$gene_name %in% gene_filter),"gene_type"] <- "oncogene"
  table(dNdScv_all_res$gene_type)

  # Change
  dNdScv_all_res$mutation_type <- ifelse(dNdScv_all_res$mutation_type == "wmis_cv", "Missense",dNdScv_all_res$mutation_type)
  dNdScv_all_res$mutation_type <- ifelse(dNdScv_all_res$mutation_type == "wnon_cv", "Nonsense",dNdScv_all_res$mutation_type)
  dNdScv_all_res$mutation_type <- ifelse(dNdScv_all_res$mutation_type == "wind_cv", "Indel",dNdScv_all_res$mutation_type)

  # Plot
  list_plots <- list() 
  combinations <- combn(names(table(dNdScv_all_res$NMDeff)), 2, simplify = FALSE)
  # for (NMD_method_char in c("END","ASE")) {
  for (NMD_method_char in c("END")) {
    for (mutation_type_char in c("Missense","Nonsense")) {
      if ( mutation_type_char == "Nonsense" ) {
        if (NMD_method_char == "ASE") { 
              ylim <- c(0,200)
        } else if (NMD_method_char == "END") {
              # ylim <- c(0,200)
              ylim <- c(0,cancer_ylims[cancer_ylims$cancer %in% cancer,"nonsense_ylim"])
        } 
      } else if ( mutation_type_char == "Missense" ) {
        if (NMD_method_char == "ASE") { 
          ylim <- c(0,80)
        } else if (NMD_method_char == "END") {
          # ylim <- c(0,80)
          ylim <- c(0,cancer_ylims[cancer_ylims$cancer %in% cancer,"missense_ylim"])
        }
      }

      # df <- dNdScv_all_res %>%
      #     filter(gene_type == "TSG") %>%
      #     filter(percentile %in% c(50)) %>%
      #     filter(NMD_method %in% NMD_method_char) %>%
      #     filter(mutation_type %in% mutation_type_char)
      # if(nrow(df)==0){list_plots[[length(list_plots) + 1]] <- ggplot(); next}

      p <- dNdScv_all_res %>%
          filter(gene_type == "TSG") %>%
          filter(percentile %in% c(50)) %>%
          filter(NMD_method %in% NMD_method_char) %>%
          filter(mutation_type %in% mutation_type_char) %>%
            ggplot(aes(x = factor(NMDeff), y = dNdScv_ratio, fill = factor(NMDeff))) +#,label = as.character(N))) +
                    geom_violin() + scale_fill_brewer(palette = "Pastel1") + xlab("") + ylab("dNdScv ratio") +
                    geom_jitter(aes(fill = factor(NMDeff)), 
                      position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0, dodge.width = 0.9),
                      alpha = 0.1, size = 2) + xlab("") +
                    labs(fill = "ETG iNMDeff") + guides(fill = guide_legend(override.aes = list(size = 12))) +
                    ggtitle(paste0(mutation_type_char)) +
                    geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
                    facet_grid(. ~ NMD_variant_type, scales = "free_y") + coord_cartesian(ylim = ylim) + 
                    geom_boxplot(width=0.3, color="black", alpha=0.2, position = position_dodge(width = 0.9)) +
                    theme_bw() +
                  stat_compare_means(
                                    aes(group = NMDeff),
                                    method.args = list(alternative = "greater"),
                                    comparisons = combinations, 
                                    size = 5, paired = FALSE,
                                    label.y = c(1.25),
                                    label = "p.format", method = "wilcox.test", hide.ns = FALSE)
        list_plots[[length(list_plots) + 1]] <- p
      }
  }

#   wilcox_test <- function(dNdScv_all_res, mutation_type_char, gene_type_char, NMD_variant_type_char,
#                           NMD_method_char, percentile_char) {
#     dNdScv_all_res_filt <- dNdScv_all_res %>% 
#               filter(mutation_type == mutation_type_char & gene_type == gene_type_char &
#                     percentile == percentile_char & NMD_method == NMD_method_char & 
#                     NMD_variant_type == NMD_variant_type_char)
#     # Subset data for 'High' and 'Low' NMDeff
#     high_group <- subset(dNdScv_all_res_filt, NMDeff == "High")
#     colnames(high_group)[colnames(high_group) %in% c("dNdScv_ratio")] <- "dNdScv_ratio_high"
#     low_group <- subset(dNdScv_all_res_filt, NMDeff == "Low")
#     colnames(low_group)[colnames(low_group) %in% c("dNdScv_ratio")] <- "dNdScv_ratio_low"
#     # Merge subsets by 'gene_name'
#     paired_data <- merge(high_group[,c("gene_name","dNdScv_ratio_high")], low_group[,c("gene_name","dNdScv_ratio_low")], by = c("gene_name"))
#     if (nrow(paired_data) == 0 ) {return(1)}       
#     wilcox_res <- wilcox.test(paired_data$dNdScv_ratio_high, paired_data$dNdScv_ratio_low, paired = TRUE, 
#                               method = "wilcox.test", alternative = "greater")
#     return(wilcox_res$p.value)
#   }

#  error <- FALSE
#     tryCatch( { 
#     # Wilcoxon tests with pairing data
#     pvalue1 <- wilcox_test(dNdScv_all_res, mutation_type_char = "Missense", gene_type_char = "TSG",
#                 NMD_variant_type_char = "NMD-evading", NMD_method_char = "END", percentile_char = 50)
#     pvalue2 <- wilcox_test(dNdScv_all_res, mutation_type_char = "Missense", gene_type_char = "TSG",
#                 NMD_variant_type_char = "NMD-triggering", NMD_method_char = "END", percentile_char = 50)
#     legend_A <- get_legend(list_plots[[1]])
#     plot_1 <- ggplot_build(list_plots[[1]] + guides(fill = FALSE))
#     plot_1$data[[5]][,"annotation"] <- as.character(plot_1$data[[5]][,"annotation"])
#     plot_1$data[[5]][,"annotation"][1:3] <- rep(round(pvalue1,2),3)
#     plot_1$data[[5]][,"annotation"][4:6] <- rep(round(pvalue2,2),3)
#     plot_1$data[[5]][,"annotation"] <- as.factor(plot_1$data[[5]][,"annotation"])
#     plot_1 <- ggplot_gtable(plot_1)
#     plot_1 <- ggplotify::as.ggplot(plot_1)
#     pvalue3 <- wilcox_test(dNdScv_all_res, mutation_type_char = "Nonsense", gene_type_char = "TSG",
#                 NMD_variant_type_char = "NMD-evading", NMD_method_char = "END", percentile_char = 50)
#     pvalue4 <- wilcox_test(dNdScv_all_res, mutation_type_char = "Nonsense", gene_type_char = "TSG",
#                 NMD_variant_type_char = "NMD-triggering", NMD_method_char = "END", percentile_char = 50)
#     # plot_2 <- ggplot_build(list_plots[[2]])
#     plot_2 <- ggplot_build(list_plots[[2]] + guides(fill = FALSE))
#     plot_2$data[[5]][,"annotation"] <- as.character(plot_2$data[[5]][,"annotation"])
#     plot_2$data[[5]][,"annotation"][1:3] <- rep(round(pvalue3,2),3)
#     plot_2$data[[5]][,"annotation"][4:6] <- rep(round(pvalue4,2),3)
#     plot_2$data[[5]][,"annotation"] <- as.factor(plot_2$data[[5]][,"annotation"])
#     plot_2 <- ggplot_gtable(plot_2)
#     plot_2 <- ggplotify::as.ggplot(plot_2)

#   },error = function(e) {error <<- TRUE})

#   if (isTRUE(error)) {
#     final_fig <- plot_grid(list_plots[[1]],list_plots[[2]], ncol = 2, nrow = 1, rel_widths = c(0.5,0.5))
#   } else {
#     final_fig <- ggarrange(plot_1,plot_2, labels = c("",""), legend.grob = legend_A,
#                     font.label = list(size = 20), #widths = c(0.45,0.28,0.28),
#                     ncol=2, nrow=1, common.legend = TRUE, legend = "right")
#     final_fig <- annotate_figure(final_fig, top = text_grob("Tumor Supressor Genes", 
#                                 size = 16, face = "plain")) #+ ggplot_theme_bw()
#   }

  final_fig <- ggarrange(list_plots[[1]],list_plots[[2]], ncol = 2, nrow = 1,
                  font.label = list(size = 20), legend = "none")
  final_fig <- annotate_figure(final_fig, top = text_grob(cancer,  
                            size = 16, face = "bold"), ) + theme_classic()
  all_list_plots[[length(all_list_plots) + 1]] <- final_fig
  # final_fig <- plot_grid(list_plots[[1]],list_plots[[2]], ncol = 2, nrow = 1, rel_widths = c(0.5,0.5))
  # final_figure_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/selection/by_cancer/dNdScv_res_",cancer,".png")
  # ggsave(final_figure_path, p, width = 250, height = 100, units = "mm") 

}

final_figure <- plot_grid(all_list_plots[[1]],all_list_plots[[2]],all_list_plots[[3]], 
                          all_list_plots[[4]],all_list_plots[[5]],all_list_plots[[6]], 
                        ncol = 2, nrow = 3, rel_widths = c(0.5,0.5))
final_figure_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/selection/by_cancer/dNdScv_res_sig_cancers.png")
ggsave(final_figure_path, final_figure, width = 400, height = 250, units = "mm") 


# FDR correction of wilxocon tests
Mann_whitney_tests <- data.frame(TCGA_cancer = cancers, miss_NMD_evading = NA, miss_NMD_triggering = NA,
                                nonsense_NMD_evading = NA, nonsense_NMD_triggering = NA)

for (cancer in cancers) {

  if (cancer %in% c("LAML","MESO")) { next }
  input <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/selection/dNdScv_all_res_",cancer,".txt")
  dNdScv_all_res <- read.table(file = input, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  # Re-classify cancer genes
  dNdScv_all_res$gene_type <- ifelse(dNdScv_all_res$gene_type == "TSG",NA,dNdScv_all_res$gene_type)
  table(dNdScv_all_res$gene_type)
  gene_filter <- TSG_STOP_CGC_and_mutpanning_CGC
  dNdScv_all_res[which(is.na(dNdScv_all_res$gene_type) & dNdScv_all_res$gene_name %in% c(gene_filter)),"gene_type"] <- "TSG"
  table(dNdScv_all_res$gene_type)
  dNdScv_all_res$gene_type <- ifelse(dNdScv_all_res$gene_type == "oncogene",NA,dNdScv_all_res$gene_type)
  table(dNdScv_all_res$gene_type)
  gene_filter <- OG_GO_and_mutpanning_CGC
  dNdScv_all_res[which(is.na(dNdScv_all_res$gene_type) & dNdScv_all_res$gene_name %in% gene_filter),"gene_type"] <- "oncogene"
  table(dNdScv_all_res$gene_type)

  # Change
  dNdScv_all_res$mutation_type <- ifelse(dNdScv_all_res$mutation_type == "wmis_cv", "Missense",dNdScv_all_res$mutation_type)
  dNdScv_all_res$mutation_type <- ifelse(dNdScv_all_res$mutation_type == "wnon_cv", "Nonsense",dNdScv_all_res$mutation_type)
  dNdScv_all_res$mutation_type <- ifelse(dNdScv_all_res$mutation_type == "wind_cv", "Indel",dNdScv_all_res$mutation_type)

  Mann_whitney_test <- function(dNdScv_all_res, mutation_type_char, gene_type_char, NMD_variant_type_char,
                          NMD_method_char, percentile_char) {
    dNdScv_all_res_filt <- dNdScv_all_res %>% 
              filter(mutation_type == mutation_type_char & gene_type == gene_type_char &
                    percentile == percentile_char & NMD_method == NMD_method_char & 
                    NMD_variant_type == NMD_variant_type_char) 

    sample_size <- sum(table(dNdScv_all_res_filt$NMDeff)  <= 5)
    
    if ( nrow(dNdScv_all_res_filt) <= 1 | sample_size %in% c(1,2) | length(table(dNdScv_all_res_filt$NMDeff)) <= 1) {return(NA)}   
    mann_whitney_res <- wilcox.test(dNdScv_ratio ~ NMDeff, data = dNdScv_all_res_filt, paired = FALSE, 
                              method = "wilcox.test", alternative = "greater")
    return(mann_whitney_res$p.value)
  }

  pvalue1 <- Mann_whitney_test(dNdScv_all_res, mutation_type_char = "Missense", gene_type_char = "TSG",
            NMD_variant_type_char = "NMD-evading", NMD_method_char = "END", percentile_char = 50)
  pvalue2 <- Mann_whitney_test(dNdScv_all_res, mutation_type_char = "Missense", gene_type_char = "TSG",
            NMD_variant_type_char = "NMD-triggering", NMD_method_char = "END", percentile_char = 50)
  pvalue3 <- Mann_whitney_test(dNdScv_all_res, mutation_type_char = "Nonsense", gene_type_char = "TSG",
            NMD_variant_type_char = "NMD-evading", NMD_method_char = "END", percentile_char = 50)
  pvalue4 <- Mann_whitney_test(dNdScv_all_res, mutation_type_char = "Nonsense", gene_type_char = "TSG",
            NMD_variant_type_char = "NMD-triggering", NMD_method_char = "END", percentile_char = 50)

  Mann_whitney_tests[Mann_whitney_tests$TCGA_cancer == cancer,-1] <- c(pvalue1, pvalue2, pvalue3, pvalue4)

}

outliers <- c("")
Mann_whitney_tests <- Mann_whitney_tests[!is.na(Mann_whitney_tests$nonsense_NMD_triggering),]

Mann_whitney_tests$nonsense_NMD_triggering_FDR <- p.adjust(Mann_whitney_tests$nonsense_NMD_triggering, method = "fdr")
Mann_whitney_tests <- Mann_whitney_tests[order(Mann_whitney_tests$nonsense_NMD_triggering_FDR,Mann_whitney_tests$nonsense_NMD_triggering),]
Mann_whitney_tests[,c("TCGA_cancer","nonsense_NMD_triggering","nonsense_NMD_triggering_FDR")]

Mann_whitney_tests$nonsense_NMD_evading_FDR <- p.adjust(Mann_whitney_tests$nonsense_NMD_evading, method = "fdr")
Mann_whitney_tests <- Mann_whitney_tests[order(Mann_whitney_tests$nonsense_NMD_evading_FDR),]
Mann_whitney_tests[,c("TCGA_cancer","nonsense_NMD_evading","nonsense_NMD_evading_FDR")]


Mann_whitney_tests$miss_NMD_triggering_FDR <- p.adjust(Mann_whitney_tests$miss_NMD_triggering, method = "fdr")
Mann_whitney_tests <- Mann_whitney_tests[order(Mann_whitney_tests$miss_NMD_triggering_FDR),]
Mann_whitney_tests[,c("TCGA_cancer","miss_NMD_triggering","miss_NMD_triggering_FDR")]

Mann_whitney_tests$miss_NMD_evading_FDR <- p.adjust(Mann_whitney_tests$miss_NMD_evading, method = "fdr")
Mann_whitney_tests <- Mann_whitney_tests[order(Mann_whitney_tests$miss_NMD_evading_FDR),]
Mann_whitney_tests




#### Genes

dndscv_all_cancers <- c()
for (cancer in as.character(cancers)) {
# for (cancer in sig_cancers) {

  if (cancer %in% c("LAML","MESO")) { next }
  input <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/selection/dNdScv_all_res_",cancer,".txt")
  dNdScv_all_res <- read.table(file = input, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

  dNdScv_all_res_filt <- dNdScv_all_res %>% 
      filter(gene_type == "TSG" & percentile == 50 & NMD_variant_type == "NMD-triggering" & 
              NMD_method == "END" & mutation_type == "wnon_cv")
  if (nrow(dNdScv_all_res_filt) == 0) {next}
  dNdScv_all_res_filt$cancer <- cancer
  if (length(dndscv_all_cancers) == 0) {
    dndscv_all_cancers <- dNdScv_all_res_filt
  } else {
    # Shared columns
    missing_cols1 <- setdiff(colnames(dndscv_all_cancers), colnames(dNdScv_all_res_filt))
    missing_cols2 <- setdiff(colnames(dNdScv_all_res_filt), colnames(dndscv_all_cancers))
    dNdScv_all_res_filt[,missing_cols1] <- NA
    dndscv_all_cancers[,missing_cols2] <- NA
    dndscv_all_cancers <- rbind(dndscv_all_cancers,dNdScv_all_res_filt)
  }

}


cancers_remove <- names(table(dndscv_all_cancers$cancer)[table(dndscv_all_cancers$cancer) == 1])


gene <- "PTEN"
# Split the data by gene_name and then calculate the difference
diff_res_df <- dndscv_all_cancers %>%
  filter(gene_name == gene) %>%
  filter(n_non >= 1) %>%
  group_by(cancer) %>%
    summarize(
      dNdScv_diff = dNdScv_ratio[NMDeff == "High"] - dNdScv_ratio[NMDeff == "Low"],
      dNdScv_mean = mean(dNdScv_ratio),
    ) %>% arrange(desc(dNdScv_diff))
data.frame(diff_res_df)
cancer_order <- diff_res_df$cancer

res_df <- dndscv_all_cancers %>%
  filter(cancer %in% cancer_order) %>%
  filter(gene_name == gene) %>%
  filter(n_non >= 1) %>%
  arrange(cancer,NMDeff) %>%
  select(cancer,dNdScv_ratio,n_syn,n_non,NMDeff) %>%
  group_by(cancer) %>%
  mutate(dNdScv_mean = mean(dNdScv_ratio)) %>% as.data.frame()

# res_df$cancer <- factor(res_df$cancer, levels = unique(res_df[order(res_df$dNdScv_mean),"cancer"]))
res_df$cancer <- factor(res_df$cancer, levels = cancer_order)
# dNdScv_diff = dNdScv_ratio[NMDeff == "High"] - dNdScv_ratio[NMDeff == "Low"],


plot_B <- ggplot(data = res_df, aes(x = cancer, 
                    y = dNdScv_ratio, fill = factor(NMDeff))) +
        geom_bar(stat = "identity", position=position_dodge(width=0.9)) + #coord_flip() +
        labs(title = "Nonsense NMD-triggering in TSG", x = gene, y = "dNdS", fill = "ETG iNMDeff") + 
        theme_bw() + #ggplot_theme() +
        # scale_fill_brewer(palette = "Set2", labels = c("pos_sel" = " >1 (Positive selection)", "neg_sel" = " <1 (Negative selection)"), direction = -1) +
        theme( legend.position = "bottom",
                legend.text = element_text(colour = "black", size = 10),
                plot.margin = unit(c(1,0.5,1,0.5), "cm"),
                axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 9)
              ) +
        geom_text(aes(label = paste0(n_non,"-",n_syn)), position = position_dodge(width=0.9), 
            size = 1.5, color = "black", hjust = 0.5, vjust = -0.5)

final_figure_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/selection/by_cancer/dNdScv_res_",gene,".png")
ggsave(final_figure_path, plot_B, width = 210, height = 105, units = "mm")




##### PTEN
gene <- "PTEN"
df <- diff_res_df[,c("cancer","NMDeff","dNdScv_ratio","n_non","n_syn","ptrunc_cv")]
df <- df %>%
arrange(cancer,NMDeff)

plot_B <- ggplot(data = df, aes(x = cancer, 
                    y = dNdScv_ratio, fill = factor(NMDeff))) +
        geom_bar(stat = "identity", position=position_dodge(width=0.9)) + #coord_flip() +
        labs(title = "Nonsense NMD-triggering in TSG", x = gene, y = "dNdS", fill = "ETG iNMDeff") + 
        theme_bw() + #ggplot_theme() +
        # scale_fill_brewer(palette = "Set2", labels = c("pos_sel" = " >1 (Positive selection)", "neg_sel" = " <1 (Negative selection)"), direction = -1) +
        theme( legend.position = "bottom",
                legend.text = element_text(colour = "black", size = 10),
                plot.margin = unit(c(1,0.5,1,0.5), "cm"),
                axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 9)
              ) +
        geom_text(aes(label = paste0(n_non,"-",n_syn)), position = position_dodge(width=0.9), 
            size = 1.5, color = "black", hjust = 0.5, vjust = -0.5)

final_figure_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/selection/by_cancer/dNdScv_res_",gene,"_all.png")
ggsave(final_figure_path, plot_B, width = 210, height = 105, units = "mm")
