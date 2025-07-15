
library(ggplot2)
library(dplyr)

# 1) Data
# 1.1) sample NMD efficiencies

# TCGA NMDeff
sample_NMD_efficiencies_TCGA_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/cancers/pancancer/NMD_efficiencies_TCGA.txt"
sample_NMD_efficiencies_TCGA <- read.table(file = sample_NMD_efficiencies_TCGA_path, header = TRUE, sep = "\t")
# Convert some columns to factors
factor_cols <- c("cancer_type","cancer_type_strat","cancer_subtype","LF_remove","purity_remove", "MSI_status",
                "batch_portion","batch_plate","batch_center","batch_vial","TCGA_full_barcode")
sample_NMD_efficiencies_TCGA[factor_cols] <- lapply(sample_NMD_efficiencies_TCGA[factor_cols], factor) 
# Remove samples with <3 PTCs
if (NMD_method %in% c("ASE")) {
    cols_rm <- paste0(NMD_method,"_stopgain_",VAF)
    sample_NMD_efficiencies_TCGA[which(sample_NMD_efficiencies_TCGA[,paste0(NMD_method,"_num_PTCs_",VAF)] < 3),cols_rm] <- NA
}
input_phenotypes <- sample_NMD_efficiencies_TCGA
# Change some columns
colnames(input_phenotypes)[1] <- "sample_short"
input_phenotypes$cancer_type <- gsub("TCGA-","",input_phenotypes$cancer_type)

# 1.2) Rare germline variants

rare_germline_variants_name <- ""

# TCGA
# 1.2) Rare germline variants
variants_gnomad_allinfo <- read.csv(file= paste('/g/strcombio/fsupek_cancer1/gpalou/Mischan/TCGA_germline_input_variants_',rare_germline_variants_name,'.txt',sep=''),
                                    head=T,sep ="\t",stringsAsFactors = F)
# GTEx

#

# 2) Enrichment analysis
# For each NMDeff threshold
# For each NMDeff method
# For each dataset

database <- "TCGA"


# Enrichment of rare pLoF in NMD genes NMDhigh vs low, comparing with random genes?

variants_genes_higher2 <- variants_gnomad_samples %>%
        mutate(samplecheck=if_else(sample_short %in% matrix_analysis$sample_short, 'yes','no')) %>%
        dplyr::filter(samplecheck=='yes') %>%
        dplyr::select(-samplecheck)  %>%
        mutate(samplecheck=if_else(Gene.refGene %in% SelectedGenes$Genes, 'yes','no')) %>%
        dplyr::filter(samplecheck=='yes') %>%
        dplyr::select(-samplecheck)  %>%
        group_by(sample_short)
samples_with_rare_variants <- unique(variants_genes_higher2$sample_short)

NMDeff_method <- "ASE_stopgain_0.2"
NMDeff_method <- "endogenous_NMD_global_2_shared"
NMDeff_quantiles <- quantile(sample_NMD_efficiencies_TCGA[,NMDeff_method],seq(0,1,0.01),na.rm=TRUE)
NMDeff_high <- NMDeff_quantiles["10%"]
NMDeff_low <- NMDeff_quantiles["90%"]
sample_NMD_efficiencies_TCGA$NMDeff_class <- NA
sample_NMD_efficiencies_TCGA[which(sample_NMD_efficiencies_TCGA[,NMDeff_method] <= NMDeff_high),"NMDeff_class"] <- "NMDeff_high"
sample_NMD_efficiencies_TCGA[which(sample_NMD_efficiencies_TCGA[,NMDeff_method] > NMDeff_low),"NMDeff_class"] <- "NMDeff_low"

df <- sample_NMD_efficiencies_TCGA[,c("sample","NMDeff_class",NMDeff_method)]
df$samples_pLOF <- NA
#aggregate(endogenous_NMD_global_2_shared ~ NMDeff_class, data=df, median)
df[which(df$sample %in% samples_with_rare_variants),"samples_pLOF"] <- "yes"
df[which(!df$sample %in% samples_with_rare_variants),"samples_pLOF"] <- "no"
#df[which(is.na(df$ASE_stopgain_0.2)),"samples_pLOF"] <- NA

#table(df$samples_pLOF)
fisher.test(table(df$samples_pLOF,df$NMDeff_class))

# 3.1) TCGA



# 3.2) GTEx