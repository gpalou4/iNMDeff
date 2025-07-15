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

dNdScv_function <- function(NMD_method, percentile, SNP_table, variants_NMDeff_TCGA_df, geneList) {

  # 1) split variants in NMD-triggering and NMD-evading
  # variants_NMDeff_TCGA_df_filt <- variants_NMDeff_TCGA_df #%>%
        # filter(VAF >= 0.1) %>%
        # filter(somatic_CNV_SNV == "no") %>%
        # filter(median_TPM_exp_transcript >= 3) %>%
        # filter(LOEUF_decile != 0) %>%
        # filter(PSG == "no")# %>%
        # filter(coeff_var <= 0.25)
  # SNVs
  NMD_evading_variants <- variants_NMDeff_TCGA_df %>%
                filter((X55_nt_last_exon == "NMD-evading" | TSS_PTC_dist < 250 | PTC_CDS_exon_length > 1000)) 
  NMD_triggering_variants <- variants_NMDeff_TCGA_df %>%
                filter(! (transcript_CDS_exon_num == 1 & splice_site_3UTR == "no") ) %>%
                filter((X55_nt_last_exon == "NMD-triggering" & TSS_PTC_dist >= 250)) %>%
                filter(PTC_CDS_exon_length <= 500)      
  # NMD_evading_variants <- variants_NMDeff_TCGA_df %>%
  #               filter((X55_nt_last_exon == "NMD-evading" | TSS_PTC_dist < 150 | PTC_CDS_exon_length > 407)) 
  # NMD_triggering_variants <- variants_NMDeff_TCGA_df %>%
  #               filter(! (transcript_CDS_exon_num == 1 & splice_site_3UTR == "no") ) %>%
  #               filter((X55_nt_last_exon == "NMD-triggering" & TSS_PTC_dist >= 150)) %>%
  #               filter(PTC_CDS_exon_length <= 407)               

  # variants_NMDeff_TCGA_df %>%  
  #   filter(stopgain == "nonsense") %>%
  #   filter((X55_nt_last_exon == "NMD-evading" | TSS_PTC_dist <= 125)) %>%
  #   filter(VAF >= 0.2) %>%
  #   filter(germline_SNV == "no") %>%
  #   filter(somatic_CNV_SNV == "no") %>%
  #   summarize(iNMDeff = median(NMD_efficiency_TPM, na.rm =TRUE))
  
  # variants_NMDeff_TCGA_df %>%  
  #   filter(stopgain == "nonsense") %>%
  #   filter((X55_nt_last_exon == "NMD-triggering" & TSS_PTC_dist > 250)) %>%
  #   filter(germline_SNV == "no") %>%
  #   filter(VAF >= 0.2) %>%
  #   summarize(iNMDeff = median(NMD_efficiency_TPM, na.rm =TRUE))

  # 1) Split samples between NMDeff high and lows
  higher_quantile <- paste0(as.character(100-percentile),"%")
  lower_quantile <- paste0(as.character(percentile),"%")

  # Exclude MSI samples
  # variants_NMDeff_TCGA_filt <- variants_NMDeff_TCGA[-which(variants_NMDeff_TCGA$MSI_status %in% c("MSI-H")),]

  samples_NMDeff <- unique(variants_NMDeff_TCGA_df[,c("TCGA_barcode","endogenous_NMD_Consensus","ASE_PTC_NMD_triggering_0.2")])

  variants_NMDeff_TCGA_df[,c("END_iNMDeff_group","ASE_iNMDeff_group")] <- NULL
  END_quantiles <- quantile(samples_NMDeff$endogenous_NMD_Consensus, seq(0,1,0.05), na.rm = TRUE)
  NMDeffHighSamples <- samples_NMDeff[which(samples_NMDeff$endogenous_NMD_Consensus >= END_quantiles[higher_quantile]),"TCGA_barcode"]
  NMDeffLowSamples <- samples_NMDeff[which(samples_NMDeff$endogenous_NMD_Consensus < END_quantiles[lower_quantile]),"TCGA_barcode"]
  variants_NMDeff_TCGA_df[which(variants_NMDeff_TCGA_df$TCGA_barcode %in% NMDeffHighSamples),"END_iNMDeff_group"] <- "NMDeff-high"
  variants_NMDeff_TCGA_df[which(variants_NMDeff_TCGA_df$TCGA_barcode %in% NMDeffLowSamples),"END_iNMDeff_group"] <- "NMDeff-low"
  table(variants_NMDeff_TCGA_df$END_iNMDeff_group)

  ASE_quantiles <- quantile(samples_NMDeff$ASE_PTC_NMD_triggering_0.2, seq(0,1,0.05), na.rm = TRUE)
  NMDeffHighSamples <- samples_NMDeff[which(samples_NMDeff$ASE_PTC_NMD_triggering_0.2 >= ASE_quantiles[higher_quantile]),"TCGA_barcode"]
  NMDeffLowSamples <- samples_NMDeff[which(samples_NMDeff$ASE_PTC_NMD_triggering_0.2 < ASE_quantiles[lower_quantile]),"TCGA_barcode"]
  variants_NMDeff_TCGA_df[which(variants_NMDeff_TCGA_df$TCGA_barcode %in% NMDeffHighSamples),"ASE_iNMDeff_group"] <- "NMDeff-high"
  variants_NMDeff_TCGA_df[which(variants_NMDeff_TCGA_df$TCGA_barcode %in% NMDeffLowSamples),"ASE_iNMDeff_group"] <- "NMDeff-low"
  table(variants_NMDeff_TCGA_df$ASE_iNMDeff_group)

  # Maybe generate our own GencodeV22 !!! buildref()
  path <- "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/dNdScv/"
  # buildref(cdsfile = paste0(path,"ensembl_v80_biomart.txt"), #excludechrs="MT",
  #         genomefile = "/g/strcombio/fsupek_cancer1/gpalou/human_genome/TCGA/GRCh38.d1.vd1.fa", 
  #         outfile = paste0(path,"RefCDS_GRCh38_ensembl_v80.rda"))

  dNdScn_res_all <- c()

  for (NMD_variant_type in c("NMD-triggering","NMD-evading")) {
    print(NMD_variant_type)
    if (NMD_variant_type == "NMD-triggering") {
      SNP_table <- NMD_triggering_variants
    } else if (NMD_variant_type == "NMD-evading") {
      SNP_table <- NMD_evading_variants
    }
    SNP_table <- SNP_table[,c("TCGA_barcode","chr","start_pos","Ref","Alt")]
    # Remove duplicates
    SNP_table <- SNP_table[!duplicated(SNP_table),]
    colnames(SNP_table) <- c("sampleID","chr","pos","ref","mut")
    # SNP_table$chr <- gsub("chr","",SNP_table$chr)

    # dNdScv 
    # Covariates for the model
    load(paste0(path,"covariates_hg19_hg38_epigenome_pcawg.rda")) # Loads the covs object

    SNP_table_NMDeff_high <- SNP_table[SNP_table$sampleID %in% unique(variants_NMDeff_TCGA_df[variants_NMDeff_TCGA_df[,paste0(NMD_method,"_iNMDeff_group")] == "NMDeff-high","TCGA_barcode"]),]
    SNP_table_NMDeff_low <- SNP_table[SNP_table$sampleID %in% unique(variants_NMDeff_TCGA_df[variants_NMDeff_TCGA_df[,paste0(NMD_method,"_iNMDeff_group")] == "NMDeff-low","TCGA_barcode"]),]
    dndsout_NMDeff_high <- dndscv(SNP_table_NMDeff_high, gene_list=geneList, max_muts_per_gene_per_sample = Inf, 
                      max_coding_muts_per_sample = Inf, sm = "192r_3w",
                      # refdb = paste0(path,"RefCDS_human_GRCh38_GencodeV18_recommended.rda"))
                      refdb = paste0(path,"RefCDS_GRCh38_ensembl_v80.rda"), cv = covs)
    dndsout_NMDeff_low <- dndscv(SNP_table_NMDeff_low, gene_list=geneList, max_muts_per_gene_per_sample = Inf, 
                      max_coding_muts_per_sample = Inf, sm = "192r_3w",
                      # refdb = paste0(path,"RefCDS_human_GRCh38_GencodeV18_recommended.rda"))
                      refdb = paste0(path,"RefCDS_GRCh38_ensembl_v80.rda"), cv = covs)
    # Results
    sel_cv_NMDeffHigh <- dndsout_NMDeff_high$sel_cv
    sel_cv_NMDeffHigh$NMDeff <- "High"
    sel_cv_NMDeffLow <- dndsout_NMDeff_low$sel_cv
    sel_cv_NMDeffLow$NMDeff <- "Low"
    print(dim(sel_cv_NMDeffHigh))
    print(dim(sel_cv_NMDeffLow))
    # Merge
    shared_cols <- intersect(colnames(sel_cv_NMDeffHigh),colnames(sel_cv_NMDeffLow))
    sel_cv_all_res <- rbind(sel_cv_NMDeffHigh[,shared_cols],sel_cv_NMDeffLow[,shared_cols])
    if ("wind_cv" %in% shared_cols) {
      sel_cv_stacked <- stack(sel_cv_all_res[,c("wmis_cv","wnon_cv","wind_cv")])
      n <- 3
    } else {
      sel_cv_stacked <- stack(sel_cv_all_res[,c("wmis_cv","wnon_cv")])
      n <- 2
    }
    keep_cols <- shared_cols[!shared_cols %in% c("wmis_cv","wnon_cv","wind_cv","wspl_cv")]
    for (col in keep_cols) {
      sel_cv_stacked[,col] <- rep(sel_cv_all_res[,col],n)
    }
    sel_cv_stacked <- sel_cv_stacked[order(sel_cv_stacked$gene_name),]
    # sel_cv_stacked$NMDeff <- rep(sel_cv_all_res$NMDeff,n)
    # sel_cv_stacked$gene_name <- rep(sel_cv_all_res$gene_name,n)
    # Add TSG/OG info
    sel_cv_stacked <- merge(sel_cv_stacked,CGC_filt[,c("Gene.Symbol","Role.in.Cancer")], by.x = "gene_name", by.y  = "Gene.Symbol", all.x = TRUE)
    table(sel_cv_stacked$Role.in.Cancer)
    colnames(sel_cv_stacked)[colnames(sel_cv_stacked) %in% c("values","ind","Role.in.Cancer")] <- c("dNdScv_ratio","mutation_type","gene_type")
    # colnames(sel_cv_stacked) <- c("gene_name","dNdScv_ratio","mutation_type","NMDeff","gene_type")
    # Remove 0 values
    sel_cv_stacked <- sel_cv_stacked[sel_cv_stacked$dNdScv_ratio != 0,]
    #sort(unique(sel_cv_stacked[which(sel_cv_stacked$Role.in.Cancer == "oncogene, TSG, fusion"),"gene_name"]))
    #aggregate(dNdScv_ratio ~ gene_type + mutation_type, data = sel_cv_stacked, median)
    sel_cv_stacked$NMD_variant_type <- NMD_variant_type
    if (length(dNdScn_res_all) == 0) {
      dNdScn_res_all <- sel_cv_stacked
    } else {
      # Shared columns
      missing_cols <- setdiff(colnames(dNdScn_res_all), colnames(sel_cv_stacked))
      sel_cv_stacked[,missing_cols] <- NA
      dNdScn_res_all <- rbind(dNdScn_res_all,sel_cv_stacked)
    }
    # Boxplot
    # Apply the function to each unique gene_type in the data
    # unique_gene_types <- na.omit(unique(sel_cv_stacked$gene_type))
    # purrr::walk(unique_gene_types, save_plot_for_gene_type, NMD_method = NMD_method, percentile = percentile, data = sel_cv_stacked)
  }
  return(dNdScn_res_all)
}

# Define a function to create and save the plot for each gene_type
save_plot_for_gene_type <- function(gene_type, NMD_method, percentile, data) {
  # Filter data for the specific gene_type
  plot_data <- data %>% filter(gene_type == !!gene_type)
  combinations <- combn(names(table(plot_data$NMDeff)), 2, simplify = FALSE)
  # Create the plot
  p <- ggplot(plot_data , aes(x = factor(NMDeff), y = dNdScv_ratio, fill = factor(NMDeff))) +#,label = as.character(N))) +
              geom_violin() + scale_fill_brewer(palette = "Dark2") + xlab("") + ylab("dNdScv ratio") +
              geom_jitter(aes(fill = factor(NMDeff)), 
                position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0, dodge.width = 0.9),
                alpha = 0.35, size = 2) +
              #ggtitle(paste0("iNMDeff for the CNA-PC",PC," groups")) +
              #geom_text(fontface = "bold", size = 10) +
              #scale_x_discrete(labels = c("High", "Mid", "Low")) +
              geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
              facet_wrap( ~ mutation_type, scales = "free") + coord_cartesian(ylim = c(0,10)) + 
              geom_boxplot(width=0.3, color="black", alpha=0.2) +
              theme_bw(base_size = 30) +
              theme(plot.title = element_text(hjust = 0.5, size = 55, margin = margin(t = 20, b = 20)),,
                    axis.text.x = element_text(size = 50, hjust = 0.5, vjust = 0.5),
                    axis.title.x = element_text(size = 55),
                    axis.text.y = element_text(size = 50),
                    axis.title.y = element_text(size = 55),
                    strip.text = element_text(size = 50),
                    legend.position = "none") +
            stat_compare_means(comparisons = combinations, size = 12,
                              label.y = c(1.25),
                              label = "p.format", method = "wilcox.test", hide.ns = TRUE)
  # Define the file name based on the gene_type
  file_name <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/selection/",NMD_method,"_dNdScv_boxplot_", gene_type,"_perc_",percentile,".png")
  # Save the plot to the file
  ggsave(filename = file_name, plot = p, width = 7500, height = 5000, units = "px")
}

# Libraries
library(dplyr)
# library(tidyverse)
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

# 1.2) PTC-NMDeff
# PTC
PTCs_NMD_efficiencies_TCGA_path <- "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/PTCs_dataset/somatic_PTCs_all_TCGA_seq.txt"
PTCs_NMD_efficiencies_TCGA <- read.table(file = PTCs_NMD_efficiencies_TCGA_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# ASE
# PTCs_ASE_NMD_efficiencies_TCGA_path <- "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/PTCs_dataset/germline_PTCs_ASE_all_TCGA_confident.txt"
# PTCs_ASE_NMD_efficiencies_TCGA <- read.table(file = PTCs_ASE_NMD_efficiencies_TCGA_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# Add sample metadata
# PTCs_ASE_NMD_efficiencies_TCGA <- merge(PTCs_ASE_NMD_efficiencies_TCGA,sample_NMD_efficiencies_TCGA, by.x = "TCGA_barcode", by.y = "sample", all.x = TRUE)
# PTCs_NMD_efficiencies_TCGA <- merge(PTCs_NMD_efficiencies_TCGA,sample_NMD_efficiencies_TCGA, by.x = "TCGA_barcode", by.y = "sample", all.x = TRUE)

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

# 1.5) Missense/syn NMDeff
MisSyn_NMD_efficiencies_TCGA_path <- "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/PTCs_dataset/somatic_MisSyn_all_TCGA.txt"
MisSyn_NMD_efficiencies_TCGA <- read.table(file = MisSyn_NMD_efficiencies_TCGA_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# ENSG00000134982 (APC) in COAD
# df <- MisSyn_NMD_efficiencies_TCGA %>%
#     filter(variant_type == "synonymous_SNV" & gene_id == "ENSG00000134982" & TCGA_cancer == "TCGA-COAD")
# df <- PTCs_NMD_efficiencies_TCGA %>%
#     filter(stopgain == "nonsense" & gene_id == "ENSG00000134982" & TCGA_cancer == "TCGA-COAD")

# 1.6) Solimini TSG STOP genes
STOP_genes <- read_excel("/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/dNdScv/Solimini_NL_et_al_2016.xls", 
                sheet = "Table S7", skip = 1)
STOP_genes <- na.omit(STOP_genes[order(STOP_genes[,"Average Log2 Ratio"],decreasing = TRUE),])
STOP_genes_top <- unique(STOP_genes[,c("Gene Symbol(s)")]) %>% pull()
STOP_genes_top_200 <- STOP_genes_top[1:200]
STOP_genes_top_200 <- gsub(" ","",unlist(strsplit(STOP_genes_top_200,",")))
STOP_genes_top_all <- gsub(" ","",unlist(strsplit(STOP_genes_top,",")))

# Manual check
# for (gene in TSG_STOP_CGC) {
#   res <- grep(paste0(gene,"$"),STOP_genes[,c("Gene Symbol(s)")][[1]])
#   if (length(res)==0) {
#     next
#   } else {
#     mean_log2_ratio <- mean(STOP_genes[res,"Average Log2 Ratio"][[1]])
#     print(gene)
#     print(mean_log2_ratio)
#   }
# }

# 1.7) Solimini Oncogenes GO genes
OG_GO <- read_excel("/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/dNdScv/Solimini_NL_et_al_2016.xls", 
                sheet = "Table S11", skip = 1)
OG_GO <- data.frame(OG_GO)
OG_GO <- OG_GO$Gene.Symbol
#X..Cell.Lines.in.which.shRNA.is.Lethal

# 1.8) Selection of cancer genes
# We will use genes from CGC to do the dNdScv test. 
# Then we will filter for interested genes (Mutpanning, Solimin, etc) in the plots
# 1.8.1) Oncogenes: from CGC
oncogenes_CGC <- as.character(na.omit(CGC_filt[CGC_filt$Role.in.Cancer == "oncogene","Gene.Symbol"]))
# 1.8.2) TSG, from CGC
TSG_CGC <- as.character(na.omit(CGC_filt[CGC_filt$Role.in.Cancer == "TSG","Gene.Symbol"]))
# 1.8.3) Final set of cancer genes
OG_TSG <- unique(c(oncogenes_CGC,TSG_CGC))

# 1.9) Add info
# Add missense variants from the samples where we have PTC
PTC_samples <- unique(PTCs_NMD_efficiencies_TCGA$TCGA_barcode)
PTC_genes <- unique(PTCs_NMD_efficiencies_TCGA$gene_id)
MisSyn_NMD_efficiencies_TCGA_filt <- MisSyn_NMD_efficiencies_TCGA %>%
            filter(TCGA_barcode %in% PTC_samples) %>%
            filter(gene_id %in% PTC_genes) #%>%
            #filter(variant_type == "synonymous_SNV")
rem_cols <- colnames(PTCs_NMD_efficiencies_TCGA)[!colnames(PTCs_NMD_efficiencies_TCGA) %in% colnames(MisSyn_NMD_efficiencies_TCGA_filt)]
variants_NMDeff_TCGA <- rbind(PTCs_NMD_efficiencies_TCGA[,!colnames(PTCs_NMD_efficiencies_TCGA) %in% rem_cols],MisSyn_NMD_efficiencies_TCGA_filt)
# sample NMD eff
variants_NMDeff_TCGA <- merge(variants_NMDeff_TCGA,sample_NMD_efficiencies_TCGA[,c("sample","endogenous_NMD_Consensus","ASE_PTC_NMD_triggering_0.2","cancer_type_strat","MSI_status")], by.x = "TCGA_barcode", by.y  = "sample", all.x = TRUE)

# e <- gsub("TCGA-","",variants_NMDeff_TCGA$TCGA_cancer)
# table(e == variants_NMDeff_TCGA$cancer_type_strat)

# 2) dNdScv

# 192 rates (used as default)
data("submod_192r_3w", package="dndscv")
colnames(substmodel) = c("syn","mis","non","spl")
head(substmodel)

# geneList (TSG & OG)
excludedGenes <- c("HMGN2P46", "IGH", "IGK", "IGL", "MALAT1",
        "MLLT4","CARS", "CRLF2", "FGFR1OP", "H3F3A", "H3F3B", "HIST1H3B", "HIST1H4I", "CASC5", "LHFP", "MDS2", "WHSC1", "WHSC1L1", 
        "P2RY8", "SEPT5", "SEPT6", "SEPT9", "KIAA1598", "C2orf44",
        "CDKN2A", "HRASLS5", "KIAA0947", "MKL1", "SLMO2", "SF3B14", "SNRPD2P1", "SNRPEP2", "NHP2L1", "FAM46C", "TRA", "TRB", "TRD", "YAE1D1", 
        "LOC100291317",  "DHFR", "ZNF663", "LOC645822",  "USP32", "C4orf14", "KIAA1324L", "ORC4L", "EGFP", "BAI2", "C9orf6",  "ROCK1P1", "C14orf109",  "LOC645453",  
        "CCDC74B", "CTSL1", "CCBL2", "CSDA", "LOC729222",  "PPFIBP1", "C3orf23", "KIAA2022", "DPH3P1", "C7orf41", "C13orf34", "CPS1IT", "CFLP1", "KLHDC5", "C14orf37")
genes_vector <- c("BPIL2", "AGXT2L1", "C5orf27", "LOC100289899", "LOC344382", "LOC100132780", "LOC654085", 
                  "LOC730107", "LOC100329108", "LOC100329109", "CCDC76", "KGFLP1", "C2orf29", "OR8U8", 
                  "FAM108B1", "TMEM56", "B3GNT1", "DCDC5", "C16orf75", "LOC653354", "C1orf212", "C14orf135", 
                  "LOC644563", "C10orf28", "LOC442454", "LOC100129942", "C13orf27", "JUB", "KGFLP2", "MRE11A", 
                  "LOC729080", "LOC641746", "LOC255025", "SRSF2IP", "LOC100131261", "LOC257358", "LOC100129138", 
                  "SQRDL", "CXorf61", "LOC100128010", "C1orf114", "C7orf58", "KIAA0141", "LOC402175", "ANKRD36BP2", 
                  "IMPAD1", "CXorf1", "C2orf47", "C19orf43", "RNF138P1", "C15orf24", "KIAA1161", "LOC220980", 
                  "C8orf4", "RBMY2FP", "C16orf63", "LOC642585", "C5orf41", "MST1P9", "DDX26B", "LOC100286977", 
                  "LOC100292152", "LOC100290750", "ATP5J", "C15orf37", "TSGA14", "GATS", "GPR110", "LOC643387", 
                  "AQPEP", "HNRPLL", "MUDENG", "ATP5E", "ATP5EP2", "SMCR7L", "LOC647349", "AGXT2L2", "ZAK", 
                  "METT5D1", "PTPLA", "KIAA1045", "C20orf197", "LOC150568", "APITD1", "C19orf62", "C14orf43", 
                  "FAM175A", "PAK7", "ZNF286B", "MLL", "LUCIFERASE", "LOC652726", "LOC400987", "LOC100134365", "CCDC99", "C1")
genes_vector2 <- c("C10orf79", "C4orf23", "ORC6L", "LOC731275", "LOC100288069", "LOC100294335", "LOC644397", 
                  "LOC728758", "RAD51L1", "TMEM105", "LOC284294", "C14orf45", "FKSG83", "LOC643486", "FAM173B", 
                  "C11orf66", "MARCH6", "HDHD1", "LOC441666", "KIAA1147", "SFRS15", "STON1-GTF2A1L", "IL8", 
                  "FAM63A", "GPR125", "FAM183B", "FAM55A", "LOC401805", "C3orf15", "C6orf170", "ZNF645", 
                  "LOC286512", "FAM26D", "C1orf49", "LOC286135", "LOC392452", "LOC285033", "LOC100133189", 
                  "LOC644496", "CTSLL2", "C6orf57", "C6orf165", "LOC100128889", "LOC100133775", "C7orf28B", 
                  "MARCH8", "FLJ10088", "GMCL1L", "KIAA1949", "EMX2OS", "AMZ2P1", "GBA3", "FAM22F", "FAM22A", 
                  "FAM22B", "FAM22E", "FAM22G", "C12orf51", "C15orf42", "PTENP1", "C19orf42", "C18orf10", 
                  "LOC100294340", "LOC646043", "CCRL1", "LOC283335", "PTTG3P", "LOC285741", "SGK223", "C9orf53", 
                  "LOC100133641", "DOPEY1", "LOC730167", "LOC340094", "LOC254312", "LOC646879", "FAM84A", "IL3RA", 
                  "HN1", "PMS2L11", "LOC100132832", "PMS2L5", "LOC728931", "ADRBK2", "SNHG3", "C10orf96", 
                  "LOC728323", "LOC389787", "CCDC72", "LOC729973", "LOC729255", "C17orf68", "LOC283788", 
                  "LOC653653", "HYALP1", "BRD7P3")
genes_vector3 <- c("WDR85", "C1orf163", "LRRC16A", "HIST1H1A", "TMEM8A", "LOC374491", "SLC24A6", 
                  "C14orf184", "LOC100133779", "LOC100288560", "LOC100294048", "PSIMCT-1", 
                  "FAM41AY1", "FAM41AY2", "LOC100291873", "LOC441554", "LOC728657", "LOC728725", 
                  "LOC728813", "LOC399815", "C15orf5", "LOC728211", "LOC727947", "KIAA0182", 
                  "LOC728147", "LOC728003", "LOC1720", "LOC100289134", "ZNF812", "C20orf117", 
                  "LOC400927", "LOC100129726", "PMS2L1", "LOC158376", "LOC442293", "C1orf144", 
                  "LOC100133760", "NCRNA00282", "LOC100288623", "COPG2IT1", "ARHGAP11B", 
                  "C14orf23", "LOC154822", "HCG4P6", "C20orf54", "TCL6", "LOC730323", "LOC100288713", 
                  "LOC100294407", "TPTE2P2", "TPTE2P3", "SNORA63", "LOC340970", "LOC642425", 
                  "TRIM53", "LOC221442", "C11orf30", "LOC100132992", "LOC440180", "LOC100287654", 
                  "LOC100288102", "LOC100294393", "C3orf48", "SNORD36C", "LOC339809", "KIAA0427", 
                  "LOC100292420", "LOC100133317", "KIAA0355", "KBTBD10", "GDEP", "LOC643699", 
                  "LOC727909", "LOC728047", "LOC728080", "LOC653061", "LOC653720", "LOC100288229", 
                  "LOC649489", "LOC1518", "LOC646548", "LOC152225", "NBR2", "PDXDC2P", "CSDAP1", "TTTY7B", "TTTY7")
genes_vector4 <- c("LOC100131149", "ANKRD32", "LOC728416", "LOC728627", "C11orf93", "LOC100132167", 
                   "LOC100132161", "LOC643630", "LOC100132644", "LOC100306975", "ZNF876P", "KIAA0494", 
                   "LOC100127922", "PMS2L4", "GATSL1", "GATSL2", "FAM157B", "AKAP2", "STAG3L2", "IGLL5", 
                   "HIST1H1C", "C10orf113", "LILRA3", "ACPP")

geneList <- OG_TSG[OG_TSG %in% c(excludedGenes,genes_vector,genes_vector2,genes_vector3,genes_vector4) == FALSE]

# geneList <- OG_TSG # With my new ensembl gene V80 I can use all genes

###########################
## Run dNdScv PAN-CANCER ##
###########################

dNdScv_all_res <- c()
for (percentile in c(5,10,20,30,40,50)) {
# for (percentile in c(5)) {
  for (NMD_method in c("END","ASE")) {
  # for (NMD_method in c("END")) {
    print(percentile)
    print(NMD_method)
    dNdScv_res <- dNdScv_function(NMD_method = NMD_method, percentile = percentile, 
                  variants_NMDeff_TCGA_df = variants_NMDeff_TCGA, geneList = geneList)
    dNdScv_res$percentile <- percentile
    dNdScv_res$NMD_method <- NMD_method
    print(dim(dNdScv_all_res))
    if (length(dNdScv_all_res) == 0 ) {
      dNdScv_all_res <- dNdScv_res
    } else {
      dNdScv_all_res <- rbind(dNdScv_all_res,dNdScv_res)
    }
  }
}

# Save
# output <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/selection/dNdScv_all_res_pancancer.txt"
output <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/selection/dNdScv_all_res_pancancer_subset_NMD_triggering.txt"
write.table(dNdScv_all_res, file = output, sep = "\t", quote = FALSE,
                      col.names = TRUE, row.names = FALSE)
# dNdScv_all_res <- read.table(file = output, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# df <- dNdScv_all_res %>%
#   filter(NMD_variant_type == "NMD-evading") %>%
#   filter(NMD_method == "END") %>%
#   filter(mutation_type == "wnon_cv") %>%
#   # filter(NMDeff == "Low") %>%
#   # filter(dNdScv_ratio >= 1.5) %>%
#   filter(gene_type == "TSG") %>%
#   arrange(gene_name,desc(dNdScv_ratio)) %>%
#   select(gene_name,dNdScv_ratio,n_syn,n_mis,n_non,NMDeff)

# head(df,100)

# low_higher_genes <- df %>%
#   dplyr::group_by(gene_name) %>%
#   dplyr::summarise(
#     high_ratio = dNdScv_ratio[NMDeff == "High"],
#     low_ratio = dNdScv_ratio[NMDeff == "Low"]
#   ) %>%
#   filter(low_ratio > high_ratio)
# low_higher_genes

# List of genes
# genes <- c("ACKR3", "ACVR1B", "AJUBA", "AKT3", "ALK", "ANAPC1", "APC", "ARHGAP35", "ARID1A", "ARID2",
#            "ATF1", "ATF1", "ATM", "B2M", "BAZ1A", "BCL10", "BCLAF1", "BCLAF1", "BCLAF1", "BCOR", "BMP5",
#            "BRAF", "BRCA1", "BRCA1", "BRD7", "BRIP1", "CARD11", "CASP8", "CASP8", "CDH10", "CDH10", "CDH17",
#            "CDH17", "CDK12", "CIITA", "CPEB3", "CR1", "CREBBP", "CRTC1", "CTCF", "CTCF", "CTCF", "CTNND1",
#            "CUL3", "CUL3", "CYLD", "CYLD", "CYSLTR2", "DDR2", "DDX3X", "DDX6", "DNMT3A", "EBF1", "ELF3",
#            "EPHA7", "ERBB2", "ETV1", "ETV1", "ETV6", "FANCC", "FANCF", "FANCG", "FBXW7", "FBXW7", "FGFR2",
#            "FGFR3", "FGFR4", "FGFR4", "FLCN", "FLT3", "FLT3", "FLT3", "FLT3", "FLT4", "FLT4", "FUBP1", "GATA3",
#            "GLI1", "GRM3", "GRM3", "GRM3", "H3F3B", "HIF1A", "HOOK3", "IDH1", "KAT7", "KCNT2", "KDM5C",
#            "KDM6A", "KIT", "KMT2D", "KMT2D", "KMT2D", "LRIG3", "MAP2K4", "MAP2K4", "MAP2K7", "MBD6", "MBD6",
#            "MEN1", "MEN1", "MGA", "MLH1", "MSH6", "MUC16", "MYB", "MYC", "NFATC2", "NOTCH1", "NR4A3", "NRG1",
#            "NRG1", "NRG1", "NSD1", "NUTM1", "PCBP1", "PCCA", "PDCD1LG2", "PDGFB", "PHOX2B", "PIK3R1", "PMS2",
#            "POLD1", "POLG", "POU2AF1", "PRDM1", "PRF1", "PRF1", "PRKACA", "PTCH1", "PTEN", "PTPN11", "PTPN13",
#            "PTPRD", "RBBP6", "RBM10", "RFWD3", "RHOA", "RHOH", "RNF43", "RNF43", "SDHA", "SETD1B", "SETD2",
#            "SFPQ", "SGK1", "SIRPA", "SIX2", "SMARCA4", "SMARCB1", "SMARCD1", "SMARCD1", "SOX2", "SPEN",
#            "SPEN", "SPEN", "SPEN", "STK11", "SUZ12", "TAOK1", "TAOK1", "TBL1XR1", "TBX3", "TCF12", "TET2",
#            "THRAP3", "THRAP3", "THRAP3", "THRAP3", "TRAF3", "TTK", "UBR5", "USP6", "WBP1", "WNK2", "WRN",
#            "WRN", "WRN", "YLPM1", "YLPM1", "YLPM1", "ZC3H13", "ZFP36L1", "ZNF208", "ZNF292", "ZNF292",
#            "ZNF331", "ZNF331", "ZNF331", "ZNF479")

# # Extract unique genes
# unique_genes <- unique(genes)

# # Print the unique vector
# unique_genes

# df %>%
#   filter(gene_name %in% unique_genes)

######################
### BY CANCER TYPE ###
######################

# There are some samples with NA in cancer_type_strat. Some of them we already knew (BRCA, ESCA, BLCA) from our strat classification
# Some others seems they are not found in our NMDeff table I think because of this:
# print("VCF germline file not found, skipping...")
# lost_sample <- "TCGA-HC-7736"
# PTCs_NMD_efficiencies_TCGA[PTCs_NMD_efficiencies_TCGA$TCGA_barcode %in% lost_sample,]
# MisSyn_NMD_efficiencies_TCGA_filt[MisSyn_NMD_efficiencies_TCGA_filt$TCGA_barcode %in% lost_sample,]
# sample_NMD_efficiencies_TCGA[sample_NMD_efficiencies_TCGA$sample %in% lost_sample,]

cancers <- na.omit(as.character(unique(variants_NMDeff_TCGA$cancer_type_strat)))
for (cancer in cancers) {
  print(cancer)
  variants_NMDeff_TCGA_filt <- variants_NMDeff_TCGA %>% filter(cancer_type_strat == cancer)
  dNdScv_all_res_cancers <- c()
  # for (percentile in c(5,10,20,30,40,50)) {
  for (percentile in c(50)) {
    for (NMD_method in c("END","ASE")) {
      print(percentile)
      print(NMD_method)
      dNdScv_res <- dNdScv_function(NMD_method = NMD_method, percentile = percentile, 
                    variants_NMDeff_TCGA_df = variants_NMDeff_TCGA_filt, geneList = geneList)
      dNdScv_res$percentile <- percentile
      dNdScv_res$NMD_method <- NMD_method
      print(dim(dNdScv_all_res_cancers))
      if (length(dNdScv_all_res_cancers) == 0 ) {
        dNdScv_all_res_cancers <- dNdScv_res
      } else {
        # Shared columns
        missing_cols <- setdiff(colnames(dNdScv_all_res_cancers), colnames(dNdScv_res))
        dNdScv_res[,missing_cols] <- NA
        dNdScv_all_res_cancers <- rbind(dNdScv_all_res_cancers,dNdScv_res)
      }
    }
  }
  output <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/selection/dNdScv_all_res_",cancer,".txt")
  write.table(dNdScv_all_res_cancers, file = output, sep = "\t", quote = FALSE,
                        col.names = TRUE, row.names = FALSE)
}

# 3) Selection of cancer genes for the analysis

# 3.1) TSG

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

# TSG Mutpanning genes with an incidence >=4 cancer types
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
TSG_STOP_CGC_and_mutpanning_CGC[!TSG_STOP_CGC_and_mutpanning_CGC %in% dNdScv_all_res$gene_name]

# 3.2) Oncogenes
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
OG_GO_and_mutpanning_CGC[!OG_GO_and_mutpanning_CGC %in% dNdScv_all_res$gene_name]

# 3.3) Re-classify cancer genes
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

# Save
#TSG ETG
output <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/selection/dNdScv_all_res.txt"
write.table(dNdScv_all_res, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig5/fig5A.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
# saveRDS(dNdScv_all_res, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig5/fig5A.RData")
saveRDS(dNdScv_all_res, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig5/fig5A.RData")
#OG ETG
write.table(dNdScv_all_res, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig19/SuppFig19C.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
saveRDS(dNdScv_all_res, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig19/SuppFig19C.RData")
# dNdScv_all_res <- read.table(file =  "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig5/fig5A.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
#TSG ASE
write.table(dNdScv_all_res, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig19/SuppFig19A.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
saveRDS(dNdScv_all_res, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig19/SuppFig19A.RData")
#OGs ASE
write.table(dNdScv_all_res, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig19/SuppFig19D.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
saveRDS(dNdScv_all_res, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig19/SuppFig19D.RData")

# 4) Plots

# NMD triggering
NMD_variant_type_char <- "NMD-triggering"
for (NMD_method_char in c("ASE","END")) {
  list_plots <- list()
  combinations <- combn(names(table(dNdScv_all_res$NMDeff)), 2, simplify = FALSE)
  # for (NMD_variant_type_char in c("NMD-triggering","NMD-evading")) {
    for (gene_type_char in c("TSG","oncogene")) {
      for (mutation_type_char in c("wmis_cv","wnon_cv","wind_cv")) {
        # Change ylim
        # if ( mutation_type_char == "wind_cv" ) {
        #   if (NMD_method_char == "ASE") {
        #     if (gene_type_char == "oncogene") {
        #       ylim <- c(0,45)
        #     } else {
        #       ylim <- c(0,120)
        #     }
        #   } else if (NMD_method_char == "END") {
        #     ylim <- c(0,10)
        #   }
        # } else if (NMD_method_char == "END") {
        #     if (mutation_type_char == "wnon_cv") {
        #       ylim <- c(0,4)
        #     } else if (mutation_type_char == "wmis_cv") {
        #       ylim <- c(0,2.5)
        #     } else {
        #       ylim <- c(0,5)
        #   }      
        # } else {
        #   ylim <- c(0,5)
        # }
        if ( mutation_type_char == "wind_cv" ) {
          if (NMD_method_char == "ASE") {
             if (NMD_variant_type_char == "NMD-triggering") {
                ylim <- c(0,7)
             } else if (NMD_variant_type_char == "NMD-evading") {
                ylim <- c(0,30)
             }  
          } else if (NMD_method_char == "END") {
             if (NMD_variant_type_char == "NMD-triggering") {
                ylim <- c(0,5)
             } else if (NMD_variant_type_char == "NMD-evading") {
                ylim <- c(0,20)
             }          }
        } else if ( mutation_type_char == "wnon_cv" ) {
          if (NMD_method_char == "ASE") { 
             if (NMD_variant_type_char == "NMD-triggering") {
                ylim <- c(0,4)
             } else if (NMD_variant_type_char == "NMD-evading") {
                ylim <- c(0,10)
             }
          } else if (NMD_method_char == "END") {
             if (NMD_variant_type_char == "NMD-triggering") {
                ylim <- c(0,4)
             } else if (NMD_variant_type_char == "NMD-evading") {
                ylim <- c(0,6)
             }
          }
        } else if ( mutation_type_char == "wmis_cv" ) {
          if (NMD_method_char == "ASE") { 
             if (NMD_variant_type_char == "NMD-triggering") {
                ylim <- c(0,3)
             } else if (NMD_variant_type_char == "NMD-evading") {
                ylim <- c(0,5)
             }          } else if (NMD_method_char == "END") {
             ylim <- c(0,2.5)
          }
        }

        p <- dNdScv_all_res %>%
            filter(NMD_variant_type %in% NMD_variant_type_char) %>%
            filter(percentile %in% c(20,30,40,50)) %>%
            filter(NMD_method %in% NMD_method_char) %>%
            filter(gene_type %in% gene_type_char) %>%
            filter(mutation_type %in% mutation_type_char) %>%
              ggplot(aes(x = factor(percentile), y = dNdScv_ratio, fill = factor(NMDeff))) +#,label = as.character(N))) +
                      geom_violin() + scale_fill_brewer(palette = "Dark2") + xlab("") + ylab("dNdScv ratio") +
                      geom_jitter(aes(fill = factor(NMDeff)), 
                        position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0, dodge.width = 0.9),
                        alpha = 0.1, size = 2) + xlab("Percentiles") +
                      labs( fill = "iNMDeff") + guides(fill = guide_legend(override.aes = list(size = 12))) +
                      ggtitle(paste0(gene_type_char)) +
                      geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
                      facet_wrap(. ~ mutation_type, scales = "free_y") + coord_cartesian(ylim = ylim) + 
                      geom_boxplot(width=0.5, color="black", alpha=0.2, position = position_dodge(width = 0.9)) +
                      #theme_bw(base_size = 30) +
                      theme_bw()+
                      theme(plot.title = element_text(hjust = 0.5, size = 18, margin = margin(t = 20, b = 20)),,
                            axis.text.x = element_text(size = 16, hjust = 0.5, vjust = 0.5),
                            axis.title.x = element_text(size = 18),
                            axis.text.y = element_text(size = 16),
                            axis.title.y = element_text(size = 16),
                            strip.text = element_text(size = 18),
                            legend.title = element_text(size = 18),
                            legend.text = element_text(size = 16),
                            legend.position = "top") #+
                    # stat_compare_means(comparisons = combinations, size = 12,
                    #                   label.y = c(1.25),
                    #                   label = "p.format", method = "wilcox.test", hide.ns = TRUE)
          list_plots[[length(list_plots) + 1]] <- p
      }
    }
    # Define the file name based on the gene_type
    file_name <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/selection/",NMD_method_char,"_",NMD_variant_type_char,"_dNdScv_boxplot_perc.png")
    # Save the plot to the file
    plot <- plot_grid(plotlist = list_plots, nrow = 2, ncol = 3)
          #labels = c("a", "b", "c", "d"), label_size = 18, rel_widths = c(0.6,0.6,1,1))
    ggsave(filename = file_name, plot = plot, width = 450, height = 300, units = "mm", limitsize = FALSE)
  # }
}

NMD_variant_type_char <- "NMD-evading"
for (NMD_method_char in c("ASE","END")) {
  list_plots <- list()
  combinations <- combn(names(table(dNdScv_all_res$NMDeff)), 2, simplify = FALSE)
  # for (NMD_variant_type_char in c("NMD-triggering","NMD-evading")) {
    for (gene_type_char in c("TSG","oncogene")) {
      for (mutation_type_char in c("wmis_cv","wnon_cv","wind_cv")) {
        # Change ylim
        # if ( mutation_type_char == "wind_cv" ) {
        #   if (NMD_method_char == "ASE") {
        #     if (gene_type_char == "oncogene") {
        #       ylim <- c(0,45)
        #     } else {
        #       ylim <- c(0,120)
        #     }
        #   } else {
        #     ylim <- c(0,25)
        #   }
        # } else if (NMD_method_char == "END") {
        #     if (mutation_type_char == "wnon_cv") {
        #       ylim <- c(0,15)
        #     } else if (mutation_type_char == "wmis_cv") {
        #       ylim <- c(0,2.5)
        #     } else {
        #       ylim <- c(0,5)
        #   }      
        # } else if ( mutation_type_char == "wnon_cv" ) {
        #     if (NMD_method_char == "ASE") { 
        #       ylim <- c(0,25)
        #     }
        # } else {
        #   ylim <- c(0,18)
        # }
        if ( mutation_type_char == "wind_cv" ) {
          if (NMD_method_char == "ASE") {
             if (NMD_variant_type_char == "NMD-triggering") {
                ylim <- c(0,7)
             } else if (NMD_variant_type_char == "NMD-evading") {
                ylim <- c(0,30)
             }  
          } else if (NMD_method_char == "END") {
             if (NMD_variant_type_char == "NMD-triggering") {
                ylim <- c(0,5)
             } else if (NMD_variant_type_char == "NMD-evading") {
                ylim <- c(0,20)
             }          }
        } else if ( mutation_type_char == "wnon_cv" ) {
          if (NMD_method_char == "ASE") { 
             if (NMD_variant_type_char == "NMD-triggering") {
                ylim <- c(0,4)
             } else if (NMD_variant_type_char == "NMD-evading") {
                ylim <- c(0,10)
             }
          } else if (NMD_method_char == "END") {
             if (NMD_variant_type_char == "NMD-triggering") {
                ylim <- c(0,4)
             } else if (NMD_variant_type_char == "NMD-evading") {
                ylim <- c(0,6)
             }
          }
        } else if ( mutation_type_char == "wmis_cv" ) {
          if (NMD_method_char == "ASE") { 
             if (NMD_variant_type_char == "NMD-triggering") {
                ylim <- c(0,3)
             } else if (NMD_variant_type_char == "NMD-evading") {
                ylim <- c(0,5)
             }          } else if (NMD_method_char == "END") {
             ylim <- c(0,2.5)
          }
        }
        p <- dNdScv_all_res %>%
            filter(NMD_variant_type %in% NMD_variant_type_char) %>%
            filter(percentile %in% c(20,30,40,50)) %>%
            filter(NMD_method %in% NMD_method_char) %>%
            filter(gene_type %in% gene_type_char) %>%
            filter(mutation_type %in% mutation_type_char) %>%
              ggplot(aes(x = factor(percentile), y = dNdScv_ratio, fill = factor(NMDeff))) +#,label = as.character(N))) +
                      geom_violin() + scale_fill_brewer(palette = "Dark2") + xlab("") + ylab("dNdScv ratio") +
                      geom_jitter(aes(fill = factor(NMDeff)), 
                        position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0, dodge.width = 0.9),
                        alpha = 0.1, size = 2) + xlab("Percentiles") +
                      labs( fill = "iNMDeff") + guides(fill = guide_legend(override.aes = list(size = 12))) +
                      ggtitle(paste0(gene_type_char)) +
                      geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
                      facet_wrap(. ~ mutation_type, scales = "free_y") + coord_cartesian(ylim = ylim) + 
                      geom_boxplot(width=0.5, color="black", alpha=0.2, position = position_dodge(width = 0.9)) +
                      #theme_bw(base_size = 30) +
                      theme_bw()+
                      theme(plot.title = element_text(hjust = 0.5, size = 18, margin = margin(t = 20, b = 20)),,
                            axis.text.x = element_text(size = 16, hjust = 0.5, vjust = 0.5),
                            axis.title.x = element_text(size = 18),
                            axis.text.y = element_text(size = 16),
                            axis.title.y = element_text(size = 16),
                            strip.text = element_text(size = 18),
                            legend.title = element_text(size = 18),
                            legend.text = element_text(size = 16),
                            legend.position = "top") #+
                    # stat_compare_means(comparisons = combinations, size = 12,
                    #                   label.y = c(1.25),
                    #                   label = "p.format", method = "wilcox.test", hide.ns = TRUE)
          list_plots[[length(list_plots) + 1]] <- p
      }
    }
    # Define the file name based on the gene_type
    file_name <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/selection/",NMD_method_char,"_",NMD_variant_type_char,"_dNdScv_boxplot_perc.png")
    # Save the plot to the file
    plot <- plot_grid(plotlist = list_plots, nrow = 2, ncol = 3)
          #labels = c("a", "b", "c", "d"), label_size = 18, rel_widths = c(0.6,0.6,1,1))
    ggsave(filename = file_name, plot = plot, width = 450, height = 300, units = "mm", limitsize = FALSE)
  # }
}

# 3) Repeating Rik's analysis 2016 Nat Genetics
# Selection in OG and TSG difference between NMD-triggering vs NMD-evading regions
# Positive selection in NMD-triggering in TSG and negative selection in OG


# 5) By-gene analysis
##### ETG #####
df <- dNdScv_all_res %>% 
      filter(mutation_type == "wnon_cv") %>%
      filter(NMD_variant_type == "NMD-triggering") %>%
      filter(percentile == 50) %>%
      filter(NMD_method == "END") %>%
      filter(gene_type == "TSG")
head(df)

# Split the data by gene_name and then calculate the difference
df_genes <- df %>%
  group_by(gene_name) %>%
  summarize(
    dNdScv_diff = dNdScv_ratio[NMDeff == "High"] - dNdScv_ratio[NMDeff == "Low"],
    dNdScv_mean = mean(dNdScv_ratio),
    n_non_High = n_non[NMDeff == "High"],
    n_non_Low = n_non[NMDeff == "Low"]
  ) %>% arrange(desc(dNdScv_diff))
data.frame(df_genes)

df_genes$dNdScv_group <- NA
df_genes$dNdScv_group <- ifelse(df_genes$dNdScv_mean > 1, "pos_sel","neg_sel")
table(df_genes$dNdScv_group)
# Save
write.table(df_genes, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig5/fig5B.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
saveRDS(df_genes, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig5/fig5B.RData")

##### ASE #####
df <- dNdScv_all_res %>% 
      filter(mutation_type == "wnon_cv") %>%
      filter(NMD_variant_type == "NMD-triggering") %>%
      filter(percentile == 50) %>%
      filter(NMD_method == "ASE") %>%
      filter(gene_type == "TSG")
head(df)

# Split the data by gene_name and then calculate the difference
df_genes <- df %>%
  group_by(gene_name) %>%
  summarize(
    dNdScv_diff = dNdScv_ratio[NMDeff == "High"] - dNdScv_ratio[NMDeff == "Low"],
    dNdScv_mean = mean(dNdScv_ratio),
    n_non_High = n_non[NMDeff == "High"],
    n_non_Low = n_non[NMDeff == "Low"]
  ) %>% arrange(desc(dNdScv_diff))
data.frame(df_genes)

df_genes$dNdScv_group <- NA
df_genes$dNdScv_group <- ifelse(df_genes$dNdScv_mean > 1, "pos_sel","neg_sel")
table(df_genes$dNdScv_group)
# Save
write.table(df_genes, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig19/SuppFig19B.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
saveRDS(df_genes, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig19/SuppFig19B.RData")



