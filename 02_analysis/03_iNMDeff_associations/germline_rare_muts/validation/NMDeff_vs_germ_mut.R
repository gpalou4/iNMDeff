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
  colnames(sample_NMDeff)[colnames(sample_NMDeff) %in% c("endogenous_NMD_global_2_shared","endogenous_NMD_global_2_shared_randomized")] <- c("NMD Consensus","NMD Consensus randomized")
  colnames(sample_NMDeff)[colnames(sample_NMDeff) %in% c("endogenous_NMD_global","endogenous_NMD_global_randomized")] <- c("NMD All","NMD All randomized")
  colnames(sample_NMDeff)[colnames(sample_NMDeff) %in% c("endogenous_NMD_Karousis","endogenous_NMD_Karousis_randomized")] <- c("NMD Karousis","NMD Karousis randomized")
  colnames(sample_NMDeff)[colnames(sample_NMDeff) %in% c("endogenous_NMD_Colombo","endogenous_NMD_Colombo_randomized")] <- c("NMD Colombo","NMD Colombo randomized")
  colnames(sample_NMDeff)[colnames(sample_NMDeff) %in% c("endogenous_non_NMD_neg_control_with_NMD_features","endogenous_non_NMD_neg_control_with_NMD_features_randomized")] <- c("RandomGenes with NMD features","RandomGenes with NMD features - randomized")
  colnames(sample_NMDeff)[colnames(sample_NMDeff) %in% c("endogenous_non_NMD_neg_control","endogenous_non_NMD_neg_control_randomized")] <- c("RandomGenes without NMD features","RandomGenes without NMD features - randomized")
  colnames(sample_NMDeff)[colnames(sample_NMDeff) %in% c("endogenous_NMD_Courtney","endogenous_NMD_Courtney_randomized")] <- c("NMD Courtney","NMD Courtney randomized")
  colnames(sample_NMDeff)[colnames(sample_NMDeff) %in% c("endogenous_NMD_ensembl","endogenous_NMD_ensembl_randomized")] <- c("NMD Ensembl","NMD Ensembl randomized")
  colnames(sample_NMDeff)[colnames(sample_NMDeff) %in% c("endogenous_SMG6","endogenous_SMG6_randomized")] <- c("NMD SMG6","NMD SMG6 randomized")
  colnames(sample_NMDeff)[colnames(sample_NMDeff) %in% c("endogenous_SMG7","endogenous_SMG7_randomized")] <- c("NMD SMG7","NMD SMG7 randomized")
  colnames(sample_NMDeff)[colnames(sample_NMDeff) %in% c("endogenous_NMD_Tani","endogenous_NMD_Tani_randomized")] <- c("NMD Tani","NMD Tani randomized")
  # ASE
  colnames(sample_NMDeff)[colnames(sample_NMDeff) %in% c("ASE_stopgain_0.2","ASE_stopgain_0.2_randomized")] <- c("PTC NMD-triggering 0.2","PTC NMD-triggering 0.2 randomized")
  colnames(sample_NMDeff)[colnames(sample_NMDeff) %in% c("ASE_stopgain_0.01","ASE_stopgain_0.01_randomized")] <- c("PTC NMD-triggering 0.01","PTC NMD-triggering 0.01 randomized")
  colnames(sample_NMDeff)[colnames(sample_NMDeff) %in% c("ASE_stopgain_NMD_evading_0.2","ASE_stopgain_NMD_evading_0.2_randomized")] <- c("PTC NMD-evading 0.2","PTC NMD-evading 0.2 randomized")
  colnames(sample_NMDeff)[colnames(sample_NMDeff) %in% c("ASE_stopgain_NMD_evading_0.01","ASE_stopgain_NMD_evading_0.01_randomized")] <- c("PTC NMD-evading 0.01","PTC NMD-evading 0.01 randomized")
  colnames(sample_NMDeff)[colnames(sample_NMDeff) %in% c("ASE_synonymous_0.2","ASE_synonymous_0.2_randomized")] <- c("Synonymous 0.2","Synonymous 0.2 randomized")
  colnames(sample_NMDeff)[colnames(sample_NMDeff) %in% c("ASE_synonymous_0.01","ASE_synonymous_0.01_randomized")] <- c("Synonymous 0.01","Synonymous 0.01 randomized")
  # Scale NMD genesets for the three methods
  # Change the sign (coefficients are reversed), so higher values means high NMDeff
  sample_NMDeff[,all_NMD_genesets] <- -sample_NMDeff[,all_NMD_genesets]
  # Scale
  if (isTRUE(scale)) {
      sample_NMDeff[,all_NMD_genesets] <- scale(sample_NMDeff[,all_NMD_genesets])
  }
  # Filter samples with low PTC number in ASE
  sample_NMDeff[which(sample_NMDeff$ASE_num_PTCs_0.2 < 3),c("PTC NMD-triggering 0.2","PTC NMD-evading 0.2","NMDeff_mean")] <- NA
  sample_NMDeff[which(sample_NMDeff$ASE_num_PTCs_0.01 < 3),c("PTC NMD-triggering 0.01","PTC NMD-evading 0.01")] <- NA

  return(sample_NMDeff)
}

library(ggplot2)
library(dplyr)
library(ggpubr)

# 1) Data
endogenous_NMD_genesets <-  c("NMD Colombo","NMD Karousis","NMD Tani","NMD Courtney","NMD Ensembl",
                      "NMD All","NMD Consensus","NMD SMG6","NMD SMG7",
                      "RandomGenes without NMD features","RandomGenes with NMD features")
ASE_NMD_genesets <- c("PTC NMD-triggering 0.01","PTC NMD-evading 0.01","Synonymous 0.01",
                      "PTC NMD-triggering 0.2","PTC NMD-evading 0.2","Synonymous 0.2")

# 1.1) sample NMD efficiencies GTEx
sample_NMD_efficiencies_TCGA_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/cancers/pancancer/NMD_efficiencies_TCGA.txt"
sample_NMD_efficiencies_TCGA <- read.table(file = sample_NMD_efficiencies_TCGA_path, header = TRUE, sep = "\t")
sample_NMD_efficiencies_TCGA <- modify_NMDeff_dataframe(sample_NMD_efficiencies_TCGA, dataset = "TCGA", scale = TRUE)
# Change some columns
filter <- colnames(sample_NMD_efficiencies_TCGA) %in% c("sample","sex")
colnames(sample_NMD_efficiencies_TCGA)[filter] <- c("sample_short","sex")

# 1.2) sample NMD efficiencies GTEx
sample_NMD_efficiencies_GTEx_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/GTEx/pantissue/NMD_efficiencies_GTEx.txt"
sample_NMD_efficiencies_GTEx <- read.table(file = sample_NMD_efficiencies_GTEx_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
sample_NMD_efficiencies_GTEx <- modify_NMDeff_dataframe(sample_NMD_efficiencies_GTEx, dataset = "GTEx", scale = TRUE)
# Change some columns
filter <- colnames(sample_NMD_efficiencies_GTEx) %in% c("sample","sex")
colnames(sample_NMD_efficiencies_GTEx)[filter] <- c("sample_short","sex")

# 1.3) GTEx Rare germline variants
rare_germline_variants_name <- "PTV_Missense_CADD15_0.1perc"
input_path <- paste0("/g/strcombio/fsupek_cancer1/gpalou/GTEx/WES_VCF/rare_germline_variants/GTEx_germline_input_variants_",rare_germline_variants_name,".txt")
GTEx_variants_gnomad_allinfo <- read.csv(file = input_path, head=T,sep ="\t",stringsAsFactors = F)
# Change columns
GTEx_variants_gnomad_allinfo$Gene.refGene <- GTEx_variants_gnomad_allinfo$gene_name
filter <- colnames(GTEx_variants_gnomad_allinfo) %in% c("GTEx_sample")
colnames(GTEx_variants_gnomad_allinfo)[filter] <- c("sample_short")

# 1.4) TCGA Rare germline variants
rare_germline_variants_name <- "PTV_Missense_CADD15_0.1perc"
input_path <- paste('/g/strcombio/fsupek_cancer1/gpalou/Mischan/TCGA_germline_input_variants_',rare_germline_variants_name,'.txt',sep='')
TCGA_variants_gnomad_allinfo <- read.csv(file = input_path, head=T,sep ="\t",stringsAsFactors = F)

# 2) Show NMDeff for individuals with mutations vs WT in 2-3 selected genes in both TCGA and GTEx
GTEx_tissues <- c("Brain_Nucleus_accumbens_basal_ganglia","Thyroid",
            "Brain_Nucleus_accumbens_basal_ganglia","Brain_Nucleus_accumbens_basal_ganglia|Brain_Hypothalamus",
            "Brain_Substantia_nigra","Thyroid","Brain_Nucleus_accumbens_basal_ganglia","Thyroid",
            "Brain_Nucleus_accumbens_basal_ganglia","Thyroid")
genes <- c("NUP153","KDM6B","TRAP1","CRTC1","PXDN","PDIA2","LAMC1","FIG4","COL19A1","CCDC114")
TCGA_cancers <- c("LGG","THCA","LGG","LGG","LGG","THCA","LGG","THCA","GBM","THCA")
hits_and_tissues <- data.frame(gene = genes, tissue = GTEx_tissues, cancer = TCGA_cancers)

# 3.1) GTEx
GTEx_df_stack_all <- c()
for (gene in genes) {
    print(gene)
    # Filter variants in the gene
    rare_variants_filter <- GTEx_variants_gnomad_allinfo %>% 
                    filter(Gene.refGene %in% gene)
    # Filter NMDeff samples for the tissue where we found the signiticant association
    tissue_hit <- as.character(hits_and_tissues[hits_and_tissues$gene %in% gene,"tissue"])
    print(tissue_hit)
    # sample_NMD_efficiencies_GTEx_filt <- sample_NMD_efficiencies_GTEx
    sample_NMD_efficiencies_GTEx_filt <- sample_NMD_efficiencies_GTEx[grep(tissue_hit, sample_NMD_efficiencies_GTEx$tissue),]
    # Create variable of MUT vs WT samples
    sample_NMD_efficiencies_GTEx_filt$germ_mut <- "WT"
    filter <- sample_NMD_efficiencies_GTEx_filt$sample_short %in% rare_variants_filter$sample_short
    sample_NMD_efficiencies_GTEx_filt[filter,"germ_mut"] <- "MUT"
    table(sample_NMD_efficiencies_GTEx_filt$germ_mut)
    df_stack <- stack(sample_NMD_efficiencies_GTEx_filt[,c("NMD Consensus","PTC NMD-triggering 0.2")])
    df_stack$germ_mut <- rep(sample_NMD_efficiencies_GTEx_filt$germ_mut,2)
    df_stack$tissue <- rep(sample_NMD_efficiencies_GTEx_filt$tissue,2)
    df_stack$sample <- rep(sample_NMD_efficiencies_GTEx_filt$sample_short,2)
    colnames(df_stack) <- c("NMD_efficiency","NMD_method","germ_mut","tissue","sample_short")
    df_stack$gene <- gene
    # Save
    if (length(GTEx_df_stack_all) == 0) {
        GTEx_df_stack_all <- df_stack
    } else {
        GTEx_df_stack_all <- rbind(GTEx_df_stack_all,df_stack)
    }
}
table(GTEx_df_stack_all$gene,GTEx_df_stack_all$germ_mut)

# NMDeff plot
p <- ggplot(data = GTEx_df_stack_all, aes(x = germ_mut, y = NMD_efficiency, fill = factor(germ_mut, levels = c("WT","MUT")))) +
        geom_violin(position = position_dodge(width = 0.9)) +
        geom_boxplot(width=0.3, position = position_dodge(width = 0.9), color="black", alpha=0.2) +
        # coord_cartesian(ylim = c(-2,2)) + 
        coord_flip(ylim = c(-2,2)) +
        facet_wrap(gene ~ NMD_method) +
        geom_hline(yintercept = 0, linetype="dashed", color = "red", size=0.5) +
        ylab("NMD efficiency") + ggtitle("") + xlab("") +
        theme_bw() + scale_fill_brewer(palette = "Dark2") +
        theme(plot.title = element_text(hjust = 0.5, size = 20),
            axis.title.x = element_text(color="black", size=25),
            strip.text = element_text(size = 30),
            axis.title.y = element_text(color="black", size=20),
            axis.text.x = element_text(color="black", size=15),# angle = 90, hjust = 1),
            axis.text.y = element_text(color="black", size=22),
            legend.position='top', legend.text = element_text(size = 25)) +
            guides(fill = guide_legend(title = "", override.aes = list(size = 12))) +
        stat_compare_means(aes(group=germ_mut), size = 7,
                            label.y = c(0.75,1,1.25),
                            label = "p.format", method = "wilcox.test", hide.ns = TRUE)

final_figure_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/rare_germline_variants/replication/nogenelist/effect_size/GTEx_hits_NMDeff.png")
ggsave(final_figure_path, p, width = 1050, height = 1100, units = "mm", limitsize = FALSE)

GTEx_df_stack_all$dataset <- "GTEx"
df_stack_all <- GTEx_df_stack_all

# 3.2) TCGA

TCGA_df_stack_all <- c()
for (gene in genes) {
    print(gene)
    # Filter variants in the gene
    rare_variants_filter <- TCGA_variants_gnomad_allinfo %>% 
                    filter(Gene.refGene %in% gene)
    # Filter NMDeff samples for the tissue where we found the signiticant association
    tissue_hit <- as.character(hits_and_tissues[hits_and_tissues$gene %in% gene,"cancer"])
    print(tissue_hit)
    # sample_NMD_efficiencies_TCGA_filt <- sample_NMD_efficiencies_TCGA
    sample_NMD_efficiencies_TCGA_filt <- sample_NMD_efficiencies_TCGA[grep(tissue_hit, sample_NMD_efficiencies_TCGA$cancer_type),]
    # Create variable of MUT vs WT samples
    sample_NMD_efficiencies_TCGA_filt$germ_mut <- "WT"
    filter <- sample_NMD_efficiencies_TCGA_filt$sample_short %in% rare_variants_filter$sample_short
    sample_NMD_efficiencies_TCGA_filt[filter,"germ_mut"] <- "MUT"
    table(sample_NMD_efficiencies_TCGA_filt$germ_mut)
    df_stack <- stack(sample_NMD_efficiencies_TCGA_filt[,c("NMD Consensus","PTC NMD-triggering 0.2")])
    df_stack$germ_mut <- rep(sample_NMD_efficiencies_TCGA_filt$germ_mut,2)
    df_stack$tissue <- rep(sample_NMD_efficiencies_TCGA_filt$cancer_type,2)
    df_stack$sample <- rep(sample_NMD_efficiencies_TCGA_filt$sample_short,2)
    colnames(df_stack) <- c("NMD_efficiency","NMD_method","germ_mut","tissue","sample_short")
    df_stack$gene <- gene
    # Save
    if (length(TCGA_df_stack_all) == 0) {
        TCGA_df_stack_all <- df_stack
    } else {
        TCGA_df_stack_all <- rbind(TCGA_df_stack_all,df_stack)
    }
}
table(TCGA_df_stack_all$gene,TCGA_df_stack_all$germ_mut)
TCGA_df_stack_all$dataset <- "TCGA"

p <- ggplot(data = TCGA_df_stack_all, aes(x = germ_mut, y = NMD_efficiency, fill = factor(germ_mut, levels = c("WT","MUT")))) +
        geom_violin(position = position_dodge(width = 0.9)) +
        geom_boxplot(width=0.3, position = position_dodge(width = 0.9), color="black", alpha=0.2) +
        coord_cartesian(ylim = c(-2,2)) +
        facet_wrap(gene ~ NMD_method) +
        geom_hline(yintercept = 0, linetype="dashed", color = "red", size=0.5) +
        ylab("NMD efficiency") + ggtitle("") + xlab("") +
        theme_bw() + scale_fill_brewer(palette = "Dark2") +
        theme(plot.title = element_text(hjust = 0.5, size = 20),
            axis.title.x = element_text(color="black", size=25),
            strip.text = element_text(size = 30),
            axis.title.y = element_text(color="black", size=20),
            axis.text.x = element_text(color="black", size=15),# angle = 90, hjust = 1),
            axis.text.y = element_text(color="black", size=22),
            legend.position='top', legend.text = element_text(size = 25)) +
            guides(fill = guide_legend(title = "", override.aes = list(size = 12))) +
        stat_compare_means(aes(group=germ_mut), size = 7,
                            label.y = c(0.75,1,1.25),
                            label = "p.format", method = "wilcox.test", hide.ns = TRUE)

final_figure_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/rare_germline_variants/replication/nogenelist/effect_size/TCGA_hits_NMDeff.png")
ggsave(final_figure_path, p, width = 1050, height = 1100, units = "mm")

df_stack_all <- rbind(GTEx_df_stack_all,TCGA_df_stack_all)

# Save
write.table(df_stack_all, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig4/fig4D.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
saveRDS(df_stack_all, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig4/fig4D.RData")
# Save
write.table(df_stack_all, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig18/SuppFig18A_B.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
saveRDS(df_stack_all, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig18/SuppFig18A_B.RData")
# 4) Consistensy of effect size across tissues


