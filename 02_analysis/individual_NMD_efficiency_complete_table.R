NMD_eff_list <- function(NMD_method, NB_res_tissue_path, NB_res_tissue_var, tissues, VAF = NULL, dataset) {
    NB_res_list <- list()
    for (i in seq(1:length(tissues))) {
        tissue <- tissues[i]
        print(tissue)
        nb_error <- FALSE
        tryCatch( {
            file_path <- gsub("\\[X\\]",tissue,paste0(NB_res_tissue_path,paths[paths$folder_or_object==NB_res_tissue_var,"path_or_filename"]))
            if (NMD_method %in% c("ASE","PTCs")) {
                file_path <- gsub("\\[X2\\]",VAF,file_path)
            } 
            #else {
            #    file_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/GTEx/tissues/",tissue,"/NB_results/",NMD_method,"/GTEx-",tissue,"_NB_coeff_outliers.txt")
            #}
            NB_res <- read.table(file = file_path, header = TRUE, sep = "\t", row.names = 1)
            if (dataset == "TCGA") {
                NB_res$cancer_subtype <- paste0(tissue,"-",NB_res$subtype)
            } else if (dataset == "GTEx") {
            }
            colnames(NB_res) <- paste0(NMD_method,"_",colnames(NB_res))
            colnames(NB_res) <- gsub("\\.","\\_",colnames(NB_res))
            NB_res_list[[tissue]] <- NB_res }
            ,error = function(e) {nb_error <<- TRUE})
        if (isTRUE(nb_error)) {
            print(nb_error)
            next
        }
    }
    return(NB_res_list)
}

# Do not activate any conda
library("TCGAbiolinks")
library("readxl")
library("dplyr")

paths <- read.table(file = "/home/gpalou/projects/NMD/scripts/NMD_efficiency/analysis_results/NMD_efficiency_table_PATHS_cluster.txt", header = TRUE, sep = ",", stringsAsFactors = FALSE)
TMB_TCGA_path <- paths[paths$folder_or_object=="TMB_TCGA_path","path_or_filename"]
TCGA_names_path <- paths[paths$folder_or_object=="TCGA_names_path","path_or_filename"]
GTEx_names_path <- paths[paths$folder_or_object=="GTEx_names_path","path_or_filename"]
cancers_stratification_metadata_path <- paths[paths$folder_or_object=="cancers_stratification_metadata_path","path_or_filename"]
TCGA_samples_metadata_path <- paths[paths$folder_or_object=="TCGA_metadata_path","path_or_filename"]
GTEx_lib_size_path <- paths[paths$folder_or_object=="GTEx_lib_size_path","path_or_filename"]
GTEx_samples_metadata_path <- paths[paths$folder_or_object=="GTEx_samples_metadata_path","path_or_filename"]

NMD_genesets <- paste0("endogenous_",c("NMD_Karousis","NMD_Colombo","NMD_Tani","non_NMD_neg_control_with_NMD_features",
                    "NMD_Courtney","NMD_ensembl","NMD_global","NMD_global_2_shared","SMG6","SMG7","non_NMD_neg_control"))
ASE_sets <- c("ASE_synonymous","ASE_stopgain","ASE_stopgain_NMD_evading")
PTC_sets <- c("PTCs_synonymous","PTCs_stopgain_NMD_triggering","PTCs_stopgain_NMD_evading")

# 2.1) TCGA
# Sample size TCGA
TCGA_RNAseq_number_samples <- read.table(file = paste0(TCGA_names_path,paths[paths$folder_or_object=="TCGA_number_samples","path_or_filename"]),
                                header = TRUE, stringsAsFactors = FALSE)
rownames(TCGA_RNAseq_number_samples) <- paste0("TCGA-",TCGA_RNAseq_number_samples$cancer)
cancer_types <- rownames(TCGA_RNAseq_number_samples)
TCGA_tissues <- read.table(file = paste0(TCGA_names_path,paths[paths$folder_or_object=="TCGA_names","path_or_filename"]),sep = "\t", stringsAsFactors=FALSE)$V1

# 2.1.1) TCGA NMD endogenous NMDeff
#args <- commandArgs(trailingOnly=TRUE)

TCGA_NB_res_cancer_path <- paths[paths$folder_or_object=="TCGA_NB_res_cancer_endogenous_path","path_or_filename"]
TCGA_NB_res_cancer_var <- "TCGA_endogenous_NB_res_cancer"
method <- "endogenous"
TCGA_endogenous_NB_res_list <- NMD_eff_list(NMD_method = "endogenous", NB_res_tissue_path = TCGA_NB_res_cancer_path, dataset = "TCGA",
                                            NB_res_tissue_var = TCGA_NB_res_cancer_var, tissues = TCGA_tissues)
# Obtain Data Frame
cols <- c(NMD_genesets,paste0("endogenous_",c("cancer_subtype","purity","LF","num_NMD_targets")))
endogenous_NB_res_list <- lapply(TCGA_endogenous_NB_res_list,function(df){df[,cols]})
endogenous_NB_res_df  <- do.call(rbind, endogenous_NB_res_list)

# 2.1.2) TCGA PTCs NB results
VAF <- 0.01
TCGA_NB_res_cancer_path <- paths[paths$folder_or_object=="TCGA_NB_res_cancer_PTCs_path","path_or_filename"]
TCGA_NB_res_cancer_var <- "TCGA_PTCs_NB_res_cancer"
method <- "PTCs"
TCGA_PTCs_NB_res_list <- NMD_eff_list(NMD_method = "PTCs", NB_res_tissue_path = TCGA_NB_res_cancer_path, dataset = "TCGA", 
                                    NB_res_tissue_var = TCGA_NB_res_cancer_var, tissues = TCGA_tissues, VAF = VAF)
# Obtain Data Frame
cols <- c("PTCs_synonymous","PTCs_stopgain_NMD_triggering","PTCs_stopgain_NMD_evading","PTCs_num_PTCs")
PTCs_NB_res_list <- lapply(TCGA_PTCs_NB_res_list,function(df){df[,cols]})
PTCs_NB_res_df  <- do.call(rbind, PTCs_NB_res_list)
colnames(PTCs_NB_res_df)[4] <- "PTC_num_PTCs"
# NAs
if (length(which(rowSums(is.na(PTCs_NB_res_df)) == 4))!= 0) {
  PTCs_NB_res_df <- PTCs_NB_res_df[-which(rowSums(is.na(PTCs_NB_res_df)) == 4),]
}

# 2.1.3) TCGA ASE NB results
TCGA_NB_res_cancer_path <- paths[paths$folder_or_object=="TCGA_NB_res_cancer_ASE_path","path_or_filename"]
TCGA_NB_res_cancer_var <- "TCGA_ASE_NB_res_cancer"
method <- "ASE"
# 2.1.3.1) VAF 0.01
ASE_VAF <- 0.01
TCGA_ASE_NB_res_list <- NMD_eff_list(NMD_method = "ASE", NB_res_tissue_path = TCGA_NB_res_cancer_path, dataset = "TCGA",
                                    NB_res_tissue_var = TCGA_NB_res_cancer_var, tissues = TCGA_tissues, VAF = ASE_VAF)
# Obtain Data Frame
#cols <- c(ASE_sets,"subtype","purity","LF","num.PTCs")
ASE_NB_res_list <- lapply(TCGA_ASE_NB_res_list,function(df){df[,c(ASE_sets,"ASE_num_PTCs")]})
ASE_NB_res_df_0.01  <- do.call(rbind, ASE_NB_res_list)
colnames(ASE_NB_res_df_0.01)[4] <- "ASE_num_PTCs"
# NAs
ASE_NB_res_df_0.01 <- ASE_NB_res_df_0.01[-which(rowSums(is.na(ASE_NB_res_df_0.01)) == 4),]
colnames(ASE_NB_res_df_0.01) <- paste0(colnames(ASE_NB_res_df_0.01),"_",ASE_VAF)
# 2.1.3.2) VAF 0.2
ASE_VAF <- 0.2
TCGA_ASE_NB_res_list <- NMD_eff_list(NMD_method = "ASE", NB_res_tissue_path = TCGA_NB_res_cancer_path, dataset = "TCGA",
                                    NB_res_tissue_var = TCGA_NB_res_cancer_var, tissues = TCGA_tissues, VAF = ASE_VAF)
# Obtain Data Frame
#cols <- c(ASE_sets,"subtype","purity","LF","num.PTCs")
ASE_NB_res_list <- lapply(TCGA_ASE_NB_res_list,function(df){df[,c(ASE_sets,"ASE_num_PTCs")]})
ASE_NB_res_df_0.2  <- do.call(rbind, ASE_NB_res_list)
colnames(ASE_NB_res_df_0.2)[4] <- "ASE_num_PTCs"
# NAs
ASE_NB_res_df_0.2 <- ASE_NB_res_df_0.2[-which(rowSums(is.na(ASE_NB_res_df_0.2)) == 4),]
colnames(ASE_NB_res_df_0.2) <- paste0(colnames(ASE_NB_res_df_0.2),"_",ASE_VAF)
# 2.1.4) Merge
ASE_NB_res_df <- merge(ASE_NB_res_df_0.2,ASE_NB_res_df_0.01, by = "row.names", all.x = TRUE)

# 2.1.5) Merge all NMDeff from all 3 methods in one data.frame
NMD_efficiencies_TCGA <- merge(endogenous_NB_res_df,ASE_NB_res_df, by.x = "row.names", by.y = "Row.names", all.x = TRUE)
NMD_efficiencies_TCGA <- merge(NMD_efficiencies_TCGA,PTCs_NB_res_df, by.x = "Row.names", by.y = "row.names", all.x = TRUE)
NMD_efficiencies_TCGA$cancer_type <- gsub("(TCGA-\\w{2,4})[_\\.].*","\\1",NMD_efficiencies_TCGA$Row.names)
NMD_efficiencies_TCGA$sample <- gsub("TCGA-\\w{2,4}[_\\.].*(TCGA.*)","\\1",NMD_efficiencies_TCGA$Row.names)
NMD_efficiencies_TCGA$cancer_type_strat <- gsub("[\\.]","",gsub("(TCGA-\\w{2,4}[_\\.].*)TCGA.*","\\1",NMD_efficiencies_TCGA$Row.names))
NMD_efficiencies_TCGA$cancer_type_strat <- gsub(" $","",NMD_efficiencies_TCGA$cancer_type_strat)
NMD_efficiencies_TCGA$cancer_type_strat <- gsub(" ","_",NMD_efficiencies_TCGA$cancer_type_strat)
colnames(NMD_efficiencies_TCGA)[colnames(NMD_efficiencies_TCGA) == "endogenous_cancer_subtype"] <- "cancer_subtype"

# 2.1.6) NMDeff mean by ASE + Endogenous
NMD_efficiencies_TCGA$NMDeff_mean <- rowMeans(scale(NMD_efficiencies_TCGA[,c('ASE_stopgain_0.2', 'endogenous_NMD_global_2_shared')]), na.rm=FALSE)

# 2.1.6) Add TMB/TIB
TMB_TCGA <- read.csv(paste0(TMB_TCGA_path,paths[paths$folder_or_object=="TMB_TCGA","path_or_filename"]), row.names = 1)
# Some samples are duplicated
TMB_TCGA <- TMB_TCGA[-which(duplicated(substr(rownames(TMB_TCGA),1,12))),]
rownames(TMB_TCGA) <- substr(rownames(TMB_TCGA),1,12)
# Add to NMD df
NMD_efficiencies_TCGA <- merge(NMD_efficiencies_TCGA,TMB_TCGA, by.x ="sample", by.y = "row.names", all.x = TRUE)
NMD_efficiencies_TCGA$Row.names <- NULL
# 2.1.7) Add CNVB, sample library size and purity/LF remove thresholds
TCGA_CNVB <- data.frame()
for (cancer in cancer_types) {
  TCGA_cancer_samples_metadata_path <- gsub("\\[X\\]",cancer,paste0(TCGA_samples_metadata_path,paths[paths$folder_or_object=="TCGA_samples_metadata","path_or_filename"]))
  TCGA_samples_metadata <- read.table(TCGA_cancer_samples_metadata_path, header = TRUE)
  TCGA_samples_metadata_filt <- TCGA_samples_metadata[,c("sample","CNV_burden","sample_lib_size","LF_remove","purity_remove")]
  if (nrow(TCGA_CNVB) == 0) {
    TCGA_CNVB <- TCGA_samples_metadata_filt
  } else {
    TCGA_CNVB <- rbind(TCGA_CNVB,TCGA_samples_metadata_filt)
  }
}
NMD_efficiencies_TCGA_filt <- merge(NMD_efficiencies_TCGA,TCGA_CNVB, by = "sample", all.x = TRUE)
NMD_efficiencies_TCGA_filt <- NMD_efficiencies_TCGA_filt[!duplicated(NMD_efficiencies_TCGA_filt),]

# 2.1.8) Add SBS signatures (POLE mutations)
TCGA_genetic_sample_phenotypes_path <- "/g/strcombio/fsupek_cancer1/gpalou/TCGA_genetic_sample_phenotypes"
SBS <- read.csv(file = paste0(TCGA_genetic_sample_phenotypes_path,"/TCGA_WES_sigProfiler_SBS_signatures_in_samples.csv"), 
                header = TRUE)
SBS[,"Sample.Names"] <- gsub("-",".", substr(SBS[,"Sample.Names"],1,12), fixed = TRUE)
SBS_filt <- SBS[SBS$Accuracy > 0.80,c("Sample.Names","SBS10a","SBS10b")]
SBS_filt <- SBS_filt[!duplicated(SBS_filt$Sample.Names),]
NMD_efficiencies_TCGA_filt <- merge(NMD_efficiencies_TCGA_filt,SBS_filt, by.x = "sample", by.y = "Sample.Names", all.x = TRUE)
NMD_efficiencies_TCGA_filt$POLE <- ""
NMD_efficiencies_TCGA_filt[which(NMD_efficiencies_TCGA_filt$SBS10a >= 50 & NMD_efficiencies_TCGA_filt$SBS10b >= 50),"POLE"] <- "SBS10ab"
NMD_efficiencies_TCGA_filt$POLE <- factor(NMD_efficiencies_TCGA_filt$POLE)
#NMD_efficiencies_TCGA_filt$MSI_POLE <- paste0(NMD_efficiencies_TCGA_filt$MSI_status,"_",NMD_efficiencies_TCGA_filt$POLE)

# 2.1.9) Add MSI labels
TCGA_MSI <- read.table(paste0(cancers_stratification_metadata_path,paths[paths$folder_or_object=="MSI_status_MSIsensor","path_or_filename"]),
                        header = TRUE, sep = "\t")
TCGA_MSI$TCGA_cancer <- paste0("TCGA-",TCGA_MSI$TCGA_cancer)
TCGA_MSI$TCGA_cancer <- NULL
NMD_efficiencies_TCGA_filt$sample <- gsub("\\.","-",NMD_efficiencies_TCGA_filt$sample)
NMD_efficiencies_TCGA_filt <- merge(NMD_efficiencies_TCGA_filt,TCGA_MSI, by.x = "sample", by.y = "TCGA_sample", all.x = TRUE)

# MSI only for UCEC, COAD, STAD
row <- NMD_efficiencies_TCGA_filt$cancer_type %in% c("TCGA-COAD","TCGA-UCEC","TCGA-STAD")
NMD_efficiencies_TCGA_filt[row,"cancer_type_MSI"] <- paste0(NMD_efficiencies_TCGA_filt[row,"cancer_type"])
NMD_efficiencies_TCGA_filt[grep("NA",NMD_efficiencies_TCGA_filt[,"cancer_type_MSI"]),"cancer_type_MSI"] <- NA

# 2.1.10) Add full TCGA barcode and split batch effects
TCGA_RNAseq_metadata <- read.table(file = "/g/strcombio/fsupek_home/gpalou/data/TCGA_fastq_RNAseq/TCGA_fastq_RNAseq_metadata.txt", header = TRUE, sep = "\t")
# Portion
TCGA_RNAseq_metadata$batch_portion <- substr(TCGA_RNAseq_metadata$TCGA_full_barcode,18,19)
# Plate
TCGA_RNAseq_metadata$batch_plate <- substr(TCGA_RNAseq_metadata$TCGA_full_barcode,22,25)
# Center
TCGA_RNAseq_metadata$batch_center <- substr(TCGA_RNAseq_metadata$TCGA_full_barcode,27,28)
# Vial
TCGA_RNAseq_metadata$batch_vial <- substr(TCGA_RNAseq_metadata$TCGA_full_barcode,16,16)
# Merge
TCGA_RNAseq_metadata_filt <- TCGA_RNAseq_metadata[,c("batch_portion","batch_plate","batch_center","batch_vial","TCGA_full_barcode","submitter_id")]
NMD_efficiencies_TCGA_filt <- merge(NMD_efficiencies_TCGA_filt,TCGA_RNAseq_metadata_filt, by.x = "sample", by.y = "submitter_id", all.x = TRUE)

# 2.1.12) Cancer stratification
# Create a unique variable: cancer_type_strat
copy <- NMD_efficiencies_TCGA_filt
# All
TCGA_subtypes <- TCGAbiolinks::PanCancerAtlas_subtypes()
DT::datatable(
    data = TCGA_subtypes,
    filter = 'top',
    options = list(scrollX = TRUE, keys = TRUE, pageLength = 5),
    rownames = FALSE
)
head(data.frame(TCGA_subtypes))
TCGA_subtypes$TCGA_sample <- substr(TCGA_subtypes$pan.samplesID,1,12)
TCGA_subtypes_filt <- TCGA_subtypes[,!colnames(TCGA_subtypes) %in% c("cancer.type","pan.samplesID")]
# Some Individuals have >1 sample and a different cancer subtype assignation, let's remove them
dup_samples <- TCGA_subtypes_filt[duplicated(TCGA_subtypes_filt$TCGA_sample),"TCGA_sample"][[1]]
TCGA_subtypes_filt <- TCGA_subtypes_filt[!TCGA_subtypes_filt$TCGA_sample %in% dup_samples,]
# Add
NMD_efficiencies_TCGA_filt <- merge(NMD_efficiencies_TCGA_filt, TCGA_subtypes_filt, by.x = c("sample"), by.y = "TCGA_sample", all.x = TRUE)

### BRCA ###
tags <- gsub("BRCA.","BRCA_",NMD_efficiencies_TCGA_filt[which(NMD_efficiencies_TCGA_filt$cancer_type %in% "TCGA-BRCA"),"Subtype_Selected"])
NMD_efficiencies_TCGA_filt[which(NMD_efficiencies_TCGA_filt$cancer_type %in% "TCGA-BRCA"),"cancer_type_strat"] <- tags
table(NMD_efficiencies_TCGA_filt$cancer_type_strat)

### BLCA ###
TCGAquery_subtype_BLCA <- TCGAquery_subtype(tumor = "blca")
TCGAquery_subtype_BLCA <- TCGAquery_subtype_BLCA[,c("patient","mRNA cluster")]
NMD_efficiencies_TCGA_filt <- merge(NMD_efficiencies_TCGA_filt, TCGAquery_subtype_BLCA, by.x = "sample", by.y= "patient", all.x = TRUE)
tags <- NMD_efficiencies_TCGA_filt[which(NMD_efficiencies_TCGA_filt$cancer_type %in% "TCGA-BLCA"),"mRNA cluster"]
tags <- paste0("BLCA_",tags)
tags <- gsub("Basal_squamous","Basal_scc",tags)
tags <- gsub("Luminal_papillary","Lum_pap",tags)
tags <- gsub("Luminal_infiltrated","Lum_inf",tags)
tags <- gsub("Luminal$","Lum",tags)
NMD_efficiencies_TCGA_filt[which(NMD_efficiencies_TCGA_filt$cancer_type %in% "TCGA-BLCA"),"cancer_type_strat"] <- tags
table(NMD_efficiencies_TCGA_filt$cancer_type_strat)
table(NMD_efficiencies_TCGA_filt[-which(is.na(NMD_efficiencies_TCGA_filt$ASE_stopgain_0.2) | NMD_efficiencies_TCGA_filt$ASE_num_PTCs_0.2 < 3),"cancer_type_strat"])
# Remove low sample size subtypes
NMD_efficiencies_TCGA_filt[NMD_efficiencies_TCGA_filt$cancer_type_strat %in% c("BLCA_Lum","BLCA_Neuronal"),"cancer_type_strat"] <- NA
NMD_efficiencies_TCGA_filt[,"mRNA cluster"] <- NULL

### SARC ###
TCGAquery_subtype_SARC <- TCGAquery_subtype(tumor = "sarc")
TCGAquery_subtype_SARC <- TCGAquery_subtype_SARC[,c("patient","short histo")]
NMD_efficiencies_TCGA_filt <- merge(NMD_efficiencies_TCGA_filt, TCGAquery_subtype_SARC, by.x = "sample", by.y= "patient", all.x = TRUE)
tags <- NMD_efficiencies_TCGA_filt[which(NMD_efficiencies_TCGA_filt$cancer_type %in% "TCGA-SARC"),"short histo"]
tags <- paste0("SARC_",tags)
tags <- gsub("STLMS","Muscle",tags)
tags <- gsub("ULMS","Muscle",tags)
tags <- gsub("DDLPS","Fat",tags)
NMD_efficiencies_TCGA_filt[which(NMD_efficiencies_TCGA_filt$cancer_type %in% "TCGA-SARC"),"cancer_type_strat"] <- tags
table(NMD_efficiencies_TCGA_filt$cancer_type_strat)
table(NMD_efficiencies_TCGA_filt[-which(is.na(NMD_efficiencies_TCGA_filt$ASE_stopgain_0.2) | NMD_efficiencies_TCGA_filt$ASE_num_PTCs_0.2 < 3),"cancer_type_strat"])
# Remove low sample size subtypes
NMD_efficiencies_TCGA_filt[NMD_efficiencies_TCGA_filt$cancer_type_strat %in% c("SARC_MFS","SARC_MPNST","SARC_NA","SARC_SS","SARC_UPS"),"cancer_type_strat"] <- NA
NMD_efficiencies_TCGA_filt[,"short histo"] <- NULL

### ESCA ###
TCGAquery_subtype_ESCA <- TCGAquery_subtype(tumor = "esca")
TCGAquery_subtype_ESCA <- TCGAquery_subtype_ESCA[,c("patient","Histological Type")]
NMD_efficiencies_TCGA_filt <- merge(NMD_efficiencies_TCGA_filt, TCGAquery_subtype_ESCA, by.x = "sample", by.y= "patient", all.x = TRUE)
tags <- NMD_efficiencies_TCGA_filt[which(NMD_efficiencies_TCGA_filt$cancer_type %in% "TCGA-ESCA"),"Histological Type"]
tags <- paste0("ESCA_",tags)
tags <- gsub("AC","ac",tags)
tags <- gsub("ESCC","scc",tags)
NMD_efficiencies_TCGA_filt[which(NMD_efficiencies_TCGA_filt$cancer_type %in% "TCGA-ESCA"),"cancer_type_strat"] <- tags
table(NMD_efficiencies_TCGA_filt$cancer_type_strat)
table(NMD_efficiencies_TCGA_filt[-which(is.na(NMD_efficiencies_TCGA_filt$ASE_stopgain_0.2) | NMD_efficiencies_TCGA_filt$ASE_num_PTCs_0.2 < 3),"cancer_type_strat"])
# Remove low sample size subtypes
NMD_efficiencies_TCGA_filt[NMD_efficiencies_TCGA_filt$cancer_type_strat %in% c("ESCA_NA"),"cancer_type_strat"] <- NA
NMD_efficiencies_TCGA_filt[,"Histological Type"] <- NULL

### HNSC HPV ###
# TCGAquery_subtype_HNSC <- TCGAquery_subtype(tumor = "hnsc")
TCGA_subtype_HNSC <- read_excel("/g/strcombio/fsupek_cancer1/gpalou/TCGA_cancer_subtype/TCGA_HNSC_smoke.xlsx",
                   sheet = "HPV_Summary", skip = 0)
TCGA_subtype_HNSC <- TCGA_subtype_HNSC[,c("Barcode","Final_HPV_Status")]
NMD_efficiencies_TCGA_filt <- merge(NMD_efficiencies_TCGA_filt, TCGA_subtype_HNSC, by.x = "sample", by.y= "Barcode", all.x = TRUE)
tags <- NMD_efficiencies_TCGA_filt[which(NMD_efficiencies_TCGA_filt$cancer_type %in% "TCGA-HNSC"),"Final_HPV_Status"]
tags <- paste0("HNSC_",tags)
tags <- gsub("Positive","HPV_pos",tags)
tags <- gsub("Negative","HPV_neg",tags)
tags <- gsub("NA","HPV_NA",tags)
NMD_efficiencies_TCGA_filt[which(NMD_efficiencies_TCGA_filt$cancer_type %in% "TCGA-HNSC"),"cancer_type_strat"] <- tags
table(NMD_efficiencies_TCGA_filt$cancer_type_strat)
table(NMD_efficiencies_TCGA_filt[-which(is.na(NMD_efficiencies_TCGA_filt$ASE_stopgain_0.2) | NMD_efficiencies_TCGA_filt$ASE_num_PTCs_0.2 < 3),"cancer_type_strat"])
# Remove low sample size subtypes
NMD_efficiencies_TCGA_filt[,"Final_HPV_Status"] <- NULL

### MSI, MSS and POLE ###
tags <- NMD_efficiencies_TCGA_filt[which(NMD_efficiencies_TCGA_filt$cancer_type %in% c("TCGA-UCEC","TCGA-COAD","TCGA-STAD")),"Subtype_Selected"]
tags <- gsub(".*\\.","",tags)
add <- gsub("TCGA-","",NMD_efficiencies_TCGA_filt[which(NMD_efficiencies_TCGA_filt$cancer_type %in% c("TCGA-UCEC","TCGA-COAD","TCGA-STAD")),"cancer_type"])
tags <- paste0(add,"_",tags)
tags <- gsub("_CIN","",tags)
tags <- gsub("_GS","",tags)
tags <- gsub("_HM-SNV","",tags)
tags <- gsub("_NA","",tags)
tags <- gsub("_EBV","",tags)
tags <- gsub("_CN_LOW","",tags)
tags <- gsub("_CN_HIGH","",tags)

NMD_efficiencies_TCGA_filt[which(NMD_efficiencies_TCGA_filt$cancer_type %in% c("TCGA-UCEC","TCGA-COAD","TCGA-STAD")),"cancer_type_strat"] <- tags
table(NMD_efficiencies_TCGA_filt$cancer_type_strat)
table(NMD_efficiencies_TCGA_filt[-which(is.na(NMD_efficiencies_TCGA_filt$ASE_stopgain_0.2) | NMD_efficiencies_TCGA_filt$ASE_num_PTCs_0.2 < 3),"cancer_type_strat"])

### GBM treatment and MGMT/IDH status ###
GBM_samples_metadata <- read_excel(paste0(cancers_stratification_metadata_path,paths[paths$folder_or_object=="TCGA_GBM_stratification","path_or_filename"])
                            , sheet = "Clinical Data", skip = 2)
GBM_samples_metadata <- GBM_samples_metadata[,c("Case ID","MGMT Status","IDH1 status","Expression\nSubclass","Therapy Class")]
colnames(GBM_samples_metadata) <- c("sample","GBM_MGMT_status","GBM_IDH1_status","GBM_expression_subclass","GBM_therapy_class")
NMD_efficiencies_TCGA_filt <- merge(NMD_efficiencies_TCGA_filt,GBM_samples_metadata, by = "sample", all.x = TRUE)

### LGG ###
LGG_samples_metadata <- read_excel(paste0(cancers_stratification_metadata_path,paths[paths$folder_or_object=="TCGA_LGG_stratification","path_or_filename"])
                            ,sheet = "Clinical.Output")
LGG_samples_metadata <- LGG_samples_metadata[,c("Tumor","primary_pharma_therapy","death","new_tumor","ProgFreeSurvEvent","histological_type","IDH/1p19q Subtype")]
colnames(LGG_samples_metadata) <- c("sample",paste0("LGG_",c("primary_pharma_therapy","death","new_tumor","ProgFreeSurvEvent","histological_type","IDH_1p19q_Subtype")))
NMD_efficiencies_TCGA_filt <- merge(NMD_efficiencies_TCGA_filt,LGG_samples_metadata, by.x = "sample", by.y = "sample", all.x = TRUE)

# 2.1.13) Add Samples metadata
TCGA_samples_metadata_2 <- read.csv(file = "/g/strcombio/fsupek_cancer1/gpalou/TCGA_metadata/TCGA_patients_info.csv", header = TRUE)
colnames(TCGA_samples_metadata_2) <- c("days_to_death","age","vital_status","disease_type","sex","race","ethnicity","sample","cancer_type","primary_site")
TCGA_samples_metadata_2 <- TCGA_samples_metadata_2[,c("days_to_death","age","vital_status","sex","race","ethnicity","sample")]
# Merge
NMD_efficiencies_TCGA_filt <- merge(NMD_efficiencies_TCGA_filt,TCGA_samples_metadata_2, by = "sample", all.x = TRUE)

# 2.1.14) Randomization of the sample-NMDeff phenotypes for each cancer type
NMD_phenotypes <- c(NMD_genesets,
                    c(paste0(ASE_sets,"_0.2"),paste0(ASE_sets,"_0.01")),
                    PTC_sets)

## Randomized each column separately
# Function to randomize column values
# randomize_column <- function(column) {
#   sample(column)
# }
# NMD_efficiencies_TCGA_final <- NMD_efficiencies_TCGA_filt %>%
#                 group_by(cancer_type) %>%
#                 mutate(across(NMD_phenotypes, randomize_column, .names = "{.col}_randomized"))

### MANUAL
# NMDeff_all_cancers_randomized <- c()
# for (cancer in unique(NMD_efficiencies_TCGA_filt$cancer_type)) { 
#   print(cancer)
#   NMDeff_cancer <- NMD_efficiencies_TCGA_filt[NMD_efficiencies_TCGA_filt$cancer_type %in% cancer,]
#   NMDeff_cancer_randomized <- NMDeff_cancer[sample(1:nrow(NMDeff_cancer)),]
#   if (length(NMDeff_all_cancers_randomized) == 0) {
#     NMDeff_all_cancers_randomized <- NMDeff_cancer_randomized
#   } else {
#     NMDeff_all_cancers_randomized <- rbind(NMDeff_all_cancers_randomized,NMDeff_cancer_randomized)
#   }
# }
# NMDeff_all_cancers_randomized <- NMDeff_all_cancers_randomized[,NMD_phenotypes]
# colnames(NMDeff_all_cancers_randomized) <- paste0(colnames(NMDeff_all_cancers_randomized),"_randomized")
# NMD_efficiencies_TCGA_final <- cbind(NMD_efficiencies_TCGA_filt,NMDeff_all_cancers_randomized)

set.seed(333)
NMD_efficiencies_TCGA_filt <- NMD_efficiencies_TCGA_filt[order(NMD_efficiencies_TCGA_filt$cancer_type),]
NMD_efficiencies_TCGA_randomized <- NMD_efficiencies_TCGA_filt %>%
                mutate(ASE_stopgain_0.2 = ifelse(ASE_num_PTCs_0.2 < 3, NA, ASE_stopgain_0.2)) %>%
                mutate(ASE_stopgain_0.01 = ifelse(ASE_num_PTCs_0.01 < 3, NA, ASE_stopgain_0.01)) %>%
                group_by(cancer_type) %>%
                select(NMD_phenotypes) %>%
                slice(sample(n())) %>%
                ungroup(cancer_type)
NMD_efficiencies_TCGA_randomized <- NMD_efficiencies_TCGA_randomized[,NMD_phenotypes]
colnames(NMD_efficiencies_TCGA_randomized) <- paste0(colnames(NMD_efficiencies_TCGA_randomized),"_randomized")
NMD_efficiencies_TCGA_final <- cbind(NMD_efficiencies_TCGA_filt,NMD_efficiencies_TCGA_randomized)

NMD_efficiencies_TCGA_final$cancer_type_strat <- gsub("TCGA-","",NMD_efficiencies_TCGA_final$cancer_type_strat)
# Save
output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/cancers/pancancer/NMD_efficiencies_TCGA.txt")
write.table(NMD_efficiencies_TCGA_final, file = output_path, sep = "\t", quote = TRUE, col.names = TRUE, row.names = FALSE)

# 2.2.1) GTEx endogenous NB results
GTEx_tissues <- read.table(file = paste0(GTEx_names_path,paths[paths$folder_or_object=="GTEx_names","path_or_filename"]),sep = "\t", stringsAsFactors=FALSE)$V1
GTEx_NB_res_tissue_path <- paths[paths$folder_or_object=="GTEx_NB_res_tissue_endogenous_path","path_or_filename"]
GTEx_NB_res_tissue_var <- "GTEx_endogenous_NB_res_tissue"
method <- "endogenous"
GTEx_endogenous_NB_res_list <- NMD_eff_list(NMD_method = "endogenous", NB_res_tissue_path = GTEx_NB_res_tissue_path, dataset = "GTEx",
                                            NB_res_tissue_var = GTEx_NB_res_tissue_var, tissues = GTEx_tissues)
# Obtain Data Frame
cols <- c(NMD_genesets,"endogenous_num_NMD_targets")
endogenous_NB_res_list <- lapply(GTEx_endogenous_NB_res_list,function(df){df[,cols]})
endogenous_NB_res_df  <- do.call(rbind, endogenous_NB_res_list)

# 2.2.2) GTEx ASE NB results
GTEx_NB_res_tissue_path <- paths[paths$folder_or_object=="GTEx_NB_res_tissue_ASE_path","path_or_filename"]
GTEx_NB_res_tissue_var <- "GTEx_ASE_NB_res_tissue"
method <- "ASE"
# 2.2.2.1) VAF 0.01
ASE_VAF <- 0.01
GTEx_ASE_NB_res_list <- NMD_eff_list(NMD_method = method, NB_res_tissue_path = GTEx_NB_res_tissue_path, dataset = "GTEx",
                                            NB_res_tissue_var = GTEx_NB_res_tissue_var, tissues = GTEx_tissues, VAF = ASE_VAF)
# Obtain Data Frame
cols <- c(ASE_sets,"ASE_num_PTCs")
ASE_NB_res_list <- lapply(GTEx_ASE_NB_res_list,function(df){df[,cols]})
ASE_NB_res_df_0.01  <- do.call(rbind, ASE_NB_res_list)
# NAs
filter <- as.numeric(which(rowSums(is.na(ASE_NB_res_df_0.01[,ASE_sets])) == length(ASE_sets)))
ASE_NB_res_df_0.01 <- ASE_NB_res_df_0.01[-filter,]
colnames(ASE_NB_res_df_0.01) <- paste0(colnames(ASE_NB_res_df_0.01),"_",ASE_VAF)

# 2.2.2.2) VAF 0.2
ASE_VAF <- 0.2
GTEx_ASE_NB_res_list <- NMD_eff_list(NMD_method = method, NB_res_tissue_path = GTEx_NB_res_tissue_path, dataset = "GTEx",
                                            NB_res_tissue_var = GTEx_NB_res_tissue_var, tissues = GTEx_tissues, VAF = ASE_VAF)
# Obtain Data Frame
cols <- c(ASE_sets,"ASE_num_PTCs")
ASE_NB_res_list <- lapply(GTEx_ASE_NB_res_list,function(df){df[,cols]})
ASE_NB_res_df_0.2  <- do.call(rbind, ASE_NB_res_list)
# NAs
filter <- as.numeric(which(rowSums(is.na(ASE_NB_res_df_0.2[,ASE_sets])) == length(ASE_sets)))
ASE_NB_res_df_0.2 <- ASE_NB_res_df_0.2[-filter,]
colnames(ASE_NB_res_df_0.2) <- paste0(colnames(ASE_NB_res_df_0.2),"_",ASE_VAF)

# 2.2.2.3) Merge
ASE_NB_res_df <- merge(ASE_NB_res_df_0.2,ASE_NB_res_df_0.01, by = "row.names", all.x = TRUE)

# # 2.2.3) GTEx PTCs NB results
PTCs_VAF <- 0.01
GTEx_NB_res_tissue_PTCs_path <- paths[paths$folder_or_object=="GTEx_NB_res_tissue_PTCs_path","path_or_filename"]
PTC_sets <- c("PTCs_synonymous","PTCs_stopgain","PTCs_stopgain_NMD_evading")
GTEx_NB_res_cancer_var <- "GTEx_PTCs_NB_res_tissue"
method <- "PTCs"
GTEx_PTCs_NB_res_list <- NMD_eff_list(NMD_method = "PTCs", NB_res_tissue_path = GTEx_NB_res_tissue_PTCs_path, dataset = "GTEx", 
                                    NB_res_tissue_var = GTEx_NB_res_cancer_var, tissues = GTEx_tissues, VAF = PTCs_VAF)
# Obtain Data Frame
cols <- c("PTCs_synonymous","PTCs_stopgain_NMD_triggering","PTCs_stopgain_NMD_evading","PTCs_num_PTCs")
PTCs_NB_res_list <- lapply(GTEx_PTCs_NB_res_list,function(df){df[,cols]})
PTCs_NB_res_df  <- do.call(rbind, PTCs_NB_res_list)
colnames(PTCs_NB_res_df)[4] <- "PTC_num_PTCs"
# NAs
if (length(which(rowSums(is.na(PTCs_NB_res_df)) == 4))!= 0) {
  PTCs_NB_res_df <- PTCs_NB_res_df[-which(rowSums(is.na(PTCs_NB_res_df)) == 4),]
}

# Merge
NMD_efficiencies_GTEx <- merge(endogenous_NB_res_df,ASE_NB_res_df, by.x = "row.names", by.y = "Row.names", all.x = TRUE)
NMD_efficiencies_GTEx <- merge(NMD_efficiencies_GTEx,PTCs_NB_res_df, by.x = "Row.names", by.y = "row.names", all.x = TRUE)
# Tissue
rownames(NMD_efficiencies_GTEx) <- NMD_efficiencies_GTEx$Row.names
NMD_efficiencies_GTEx$Row.names <- NULL
NMD_efficiencies_GTEx$sample <- gsub(".*(GTEX-.*)","\\1",rownames(NMD_efficiencies_GTEx))
NMD_efficiencies_GTEx$tissue <- gsub("(.*)\\.GTEX-.*","\\1",rownames(NMD_efficiencies_GTEx))

# 2.2.4) Acronyms
acronyms <- c("ADPSBQ","ADPVSC","ADRNLG","ARTAORT","ARTCRN","ARTTBL","BLDDER","BRNAMY","BRNACC","BRNCDT","BRNCHB","BRNCHA","BRNCTXA","BRNCTXB",
              "BRNHPP","BRNHPT","BRNNCC","BRNPTM","BRNSPC","BRNSNG","BREAST","FIBRBLS","LCL","CML","CVXECT","CVSEND","CLNSGM","CLNTRN","ESPGEJ","ESPMCS","ESPMSL",
              "FLLPNT","HRTAA","HRTLV","KDNCTX","KDNMDL","LIVER","LUNG","SLVRYG","MSCLSK","NERVET","OVARY","PNCREAS","PTTARY","PRSTTE","SKINNS","SKINS",
              "SNTTRM","SPLEEN","STMACH","TESTIS","THYROID","UTERUS","VAGINA","WHLBLD")
tissues_acronyms <- data.frame(tissues = GTEx_tissues, acronyms = acronyms )
NMD_efficiencies_GTEx <- merge(NMD_efficiencies_GTEx,tissues_acronyms, by.x = "tissue", by.y = "tissues", all.x = TRUE)

# 2.2.5) NMDeff mean by ASE + Endogenous
NMD_efficiencies_GTEx$NMDeff_mean <- rowMeans(scale(NMD_efficiencies_GTEx[,c('ASE_stopgain_0.2', 'endogenous_NMD_global_2_shared')]), na.rm=FALSE)

# 2.2.6) Sample library size
GTEx_sample_lib_size_df <- data.frame()
for (GTEx_tissue in GTEx_tissues) {
  if ( GTEx_tissue == "Cells_Leukemia_cell_line_CML") {next}
  GTEx_lib_size_tissue_path <- gsub("\\[X1\\]",GTEx_tissue,paste0(GTEx_lib_size_path,paths[paths$folder_or_object=="GTEx_lib_size","path_or_filename"]))
  GTEx_lib_size_tissue <- read.table(GTEx_lib_size_tissue_path, header = TRUE)
  GTEx_lib_size_tissue$sample <- gsub("\\.","-",GTEx_lib_size_tissue$sample)
  GTEx_lib_size_tissue$sample_filt <- gsub("(GTEX\\-\\w{4,5})\\-.*","\\1",GTEx_lib_size_tissue$sample)
  GTEx_lib_size_tissue$tissue <- GTEx_tissue
  colnames(GTEx_lib_size_tissue)[1] <- "sample_full_barcode"
  if (nrow(GTEx_sample_lib_size_df) == 0) {
    GTEx_sample_lib_size_df <- GTEx_lib_size_tissue
  } else {
    GTEx_sample_lib_size_df <- rbind(GTEx_sample_lib_size_df,GTEx_lib_size_tissue)
  }

}
NMD_efficiencies_GTEx <- merge(NMD_efficiencies_GTEx,GTEx_sample_lib_size_df, by.x = c("sample","tissue"), by.y = c("sample_filt","tissue"), all.x = TRUE)

# 2.2.7) Metadata: Age, sex
GTEx_samples_metadata_all <- paste0(GTEx_samples_metadata_path,paths[paths$folder_or_object=="GTEx_samples_metadata","path_or_filename"])
GTEx_samples_metadata <- read.table(file = GTEx_samples_metadata_all, header = TRUE, sep = "\t")
colnames(GTEx_samples_metadata) <- c("sample","sex","age","death_group")
NMD_efficiencies_GTEx <- merge(NMD_efficiencies_GTEx,GTEx_samples_metadata, by = c("sample"), all.x = TRUE)

# 2.2.8) Randomization of the sample-NMDeff phenotypes
set.seed(333)
NMD_efficiencies_GTEx <- NMD_efficiencies_GTEx[order(NMD_efficiencies_GTEx$tissue),]
NMD_efficiencies_GTEx_randomized <- NMD_efficiencies_GTEx %>%
                mutate(ASE_stopgain_0.2 = ifelse(ASE_num_PTCs_0.2 < 3, NA, ASE_stopgain_0.2)) %>%
                mutate(ASE_stopgain_0.01 = ifelse(ASE_num_PTCs_0.01 < 3, NA, ASE_stopgain_0.01)) %>%
                group_by(tissue) %>%
                select(NMD_phenotypes) %>%
                slice(sample(n())) %>%
                ungroup()
NMD_efficiencies_GTEx_randomized <- NMD_efficiencies_GTEx_randomized[,NMD_phenotypes]
colnames(NMD_efficiencies_GTEx_randomized) <- paste0(colnames(NMD_efficiencies_GTEx_randomized),"_randomized")
NMD_efficiencies_GTEx_final <- cbind(NMD_efficiencies_GTEx,NMD_efficiencies_GTEx_randomized)

# Save
output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/GTEx/pantissue/NMD_efficiencies_GTEx.txt")
write.table(NMD_efficiencies_GTEx_final, file = output_path, sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

# 2.2.5) GTEx - TCGA tissue comparisons
tissues <- c("Kidney","Liver","Thyroid","Skin","Lung","Brain","Uterus","Ovary","Prostate","Pancreas","Esophagus","Adrenal_gland","Breast","Bladder","Stomach","Colon","Muscle","Blood")
TCGA_cancers <- c("KIRC|KIRP","LIHC","THCA","SKCM","LUSC|LUAD","GBM|LGG","UCEC_MSS","OV","PRAD","PAAD","ESCA","ACC","BRCA","BLCA","STAD_MSS","COAD_MSS","SARC_muscle","LAML")
GTEx_tissues <- c("KDN","LIVER","THYROID","SKIN","LUNG","BRN","UTERUS","OVARY","PRSTTE","PNCREAS","ESP","ADR","BREAST","BLDDER","STMACH","CLN","MSCLSK","WHLBLD")
TCGA_GTEx_tissues <- data.frame(tissues,GTEx_tissues,TCGA_cancers)
output_path <- paste0("/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/TCGA_GTEx_match.txt")
write.table(TCGA_GTEx_tissues, file = output_path, sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)