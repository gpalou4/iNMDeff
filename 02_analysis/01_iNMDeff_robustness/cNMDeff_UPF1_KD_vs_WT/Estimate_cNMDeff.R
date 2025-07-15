# conda activate /home/dnaro/.conda/envs/NMD_regression_3
# library("edgeR")
.libPaths( rev( .libPaths() ) )
# Bayesian for NB regression
library("V8")
# library("shinystan")
library("rstan")
library("rstanarm")
library("readxl")
library("dplyr")
library("ggpubr")
library("ggsignif")

################################ SCRIPT ################################
################################ SCRIPT ################################
################################ SCRIPT ################################

glm_nb_regression <- function(nb_df) {

    # NMD target variable converted to numeric   
    nb_df$NMD_target <- ifelse(nb_df$NMD_target == "NMD_target", 1, 0 )
    print(head(nb_df))
    print(dim(nb_df))
    # Bayesian NB regression model
    glm_nb_model <- "stan_glm(raw_gene_exp ~ as.factor(NMD_target) + as.factor(ensembl_gene_id) + length"
    # Raw model
    glm_nb_model = paste(glm_nb_model, ", data = nb_df, family = neg_binomial_2, verbose = FALSE)")

    # Perform the regression
    tryCatch( { 
        nb_res <- NA
        nb_res <- eval(parse(text = glm_nb_model))
    }, error = function(e) {
        print("################################################################################################")
        print("ERROR in NB regression. Model:")
        print(e)
        nb_error_3 <<- TRUE
        print(glm_nb_model)
        print("################################################################################################")
    })
    return(nb_res)
}

args <- commandArgs(trailingOnly=TRUE)
NMD_gene <- args[1]
print(paste0("Processing NMD gene --> ",NMD_gene))

# 1) Data

# 1.1) Cell Line UPF1 KD
CL_UPF1KD_path <- "/g/strcombio/fsupek_cancer1/gpalou/cell_lines_UPF1_KD/ENSEMBL/rsem.merged.transcript_counts.tsv"
CL_metadata_path <- "/g/strcombio/fsupek_cancer1/gpalou/cell_lines_UPF1_KD/CL_metadata_subset_and_controls.xlsx"
CL_UPF1KD <- read.table(file = CL_UPF1KD_path, header = TRUE, sep = "\t", row.names = 1)
CL_UPF1KD <- CL_UPF1KD[,-1]
CL_UPF1KD <- CL_UPF1KD[grep(".*PAR.*",rownames(CL_UPF1KD), invert = TRUE),]
rownames(CL_UPF1KD) <- sub("\\..*","", rownames(CL_UPF1KD))
CL_metadata <- data.frame(read_excel(CL_metadata_path))
# Same order as CL data
CL_metadata$ID <- factor(paste0(CL_metadata$GSE_project,"_R",CL_metadata$replicate))
CL_metadata$Type <- as.factor(CL_metadata$Type)
CL_metadata <- CL_metadata[order(CL_metadata$ID),]
CL_UPF1KD <- CL_UPF1KD[,order(colnames(CL_UPF1KD))]

# Change names
# CL 1: HeLa # DISCARD!! #
# cols <- grep("GSE152435_C",colnames(CL_UPF1KD))
# colnames(CL_UPF1KD)[cols] <- gsub("_C_","_WT1_",colnames(CL_UPF1KD)[cols])
# cols <- grep("GSE152435_R.*",colnames(CL_UPF1KD))
# colnames(CL_UPF1KD)[cols] <- gsub("_R","_KD1_R",colnames(CL_UPF1KD)[cols])
# # CL 2: HeLa (Colombo)
# cols <- grep("GSE86148_C.*",colnames(CL_UPF1KD))
# colnames(CL_UPF1KD)[cols] <- gsub("_C_","_WT2_",colnames(CL_UPF1KD)[cols])
# cols <- grep("GSE86148_R.*",colnames(CL_UPF1KD))
# colnames(CL_UPF1KD)[cols] <- gsub("_R","_KD2_R",colnames(CL_UPF1KD)[cols])
# # CL 3: HepG2 (ENCODE)
# cols <- grep("GSE88148",colnames(CL_UPF1KD))
# colnames(CL_UPF1KD)[cols] <- gsub("_C_","_WT3_",colnames(CL_UPF1KD)[cols])
# cols <- grep("GSE88466",colnames(CL_UPF1KD))
# colnames(CL_UPF1KD)[cols] <- gsub("_R","_KD3_R",colnames(CL_UPF1KD)[cols])
# # CL 4: K562 (ENCODE)
# cols <- grep("GSE88266",colnames(CL_UPF1KD))
# colnames(CL_UPF1KD)[cols] <- gsub("_C_","_WT4",colnames(CL_UPF1KD)[cols])
# cols <- grep("GSE88140",colnames(CL_UPF1KD))
# colnames(CL_UPF1KD)[cols] <- gsub("_R","_KD4_R",colnames(CL_UPF1KD)[cols])

head(CL_UPF1KD)

# 1.2) Transcripts length
ensembl_v88_transcripts_length <- read.table(file = "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/ensembl_v88_transcripts_length.txt", 
                                header = TRUE, sep = "\t")
# 1.3) Protein coding transcripts
ensembl_v88_coding_transcripts <- read.table(file = "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/ensembl_v88_coding_transcripts.txt",
                                  header = FALSE, sep = "\t",colClasses = "vector")   

# 1.4) NMD targets Consensus
path <- "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/NMD_targets/ENSEMBL_transcripts/NMD_global_2_shared_ensembl_final_old.txt"
NMD_targets <- read.table(file = path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# 2) Filters

# 2.1) Filter low-expressed genes
CL_UPF1KD_filt <- CL_UPF1KD[rowSums(log2(CL_UPF1KD) >= 1) >= round(length(colnames(CL_UPF1KD)) * 0.50),]

# 2.2) Filter out non-coding transcripts
CL_UPF1KD_filt <- CL_UPF1KD_filt[ rownames(CL_UPF1KD_filt) %in% 
                  unique(ensembl_v88_coding_transcripts$V1),]                                   
dim(CL_UPF1KD_filt)

# 2.3) NMD targets
CL_UPF1KD_filt <- CL_UPF1KD_filt[rownames(CL_UPF1KD_filt) %in% NMD_targets$ensembl_transcript_id,]

# 3) Negative binomial regression to estimate cell NMD efficiency (cNMDeff)

# Remove one gene at a time to compare with original NMDtargets-control dataframe

# NMD_genes <- c("all",unique(NMD_targets$ensembl_gene_id))

all_nb_coeff_res <- c()

# for (NMD_gene in NMD_genes) {
print(NMD_gene)
# print((which(NMD_genes %in% NMD_gene)/length(NMD_genes) ) * 100)
if (NMD_gene == "all") {
  NMD_targets_filt <- NMD_targets
} else {
  NMD_targets_filt <- NMD_targets[!NMD_targets$ensembl_gene_id %in% NMD_gene,]
}

for (cell_line in colnames(CL_UPF1KD_filt)) {

    print(cell_line)
    nb_coeff_res <- data.frame(sample = NA, num_NMD_targets = NA, error = NA, nb_coeff = NA,
            p_value = NA, SD = NA, CI_2.5 = NA, CI_97.5 = NA)

    # Create the data frame for the regression
    nb_df <- CL_UPF1KD_filt[,colnames(CL_UPF1KD_filt) %in% cell_line, drop = FALSE]
    # Add ensembl gene ID
    nb_df <- merge(nb_df,NMD_targets[,c("ensembl_transcript_id","ensembl_gene_id","final_consensus")], 
            by.x = "row.names", by.y = "ensembl_transcript_id", all.x = TRUE)
    # Add transcript length
    nb_df <- merge(nb_df,ensembl_v88_transcripts_length, 
            by.x = "Row.names", by.y = "transcript_id", all.x = TRUE)

    # Colnames
    colnames(nb_df) <- c("ensembl_transcript_id","raw_gene_exp","ensembl_gene_id","NMD_target","length")
    head(nb_df)
    nb_df$raw_gene_exp <- round(nb_df$raw_gene_exp)

    # NB regression --> cNMDeff
    nb_res <- glm_nb_regression(nb_df = nb_df)
    # Save output
    if (is.na(nb_res)) {
        nb_res_pval <- NA
        nb_res_coeff <- NA
        nb_coeff_res[1,"error"] <- "NB regression error" 
        nb_coeff_res[1,"nb_coeff"] <- nb_res_coeff
        nb_coeff_res[1,"SD"] <- NA
        nb_coeff_res[1,"CI_2.5"] <- NA
        nb_coeff_res[1,"CI_97.5"] <- NA
    } else if (!is.na(nb_res)) {
        nb_coeff <- as.numeric(coef(nb_res)[2])
        # Store NB coefficient
        nb_coeff_res[1,"nb_coeff"] <- nb_coeff
        nb_coeff_res[1,"SD"] <- as.numeric(nb_res$stan_summary[2,"sd"])
        #nb_coeff_res[1,paste0(NMD.gene.set.name,".pval")] <- nb.res.coeff
        nb_coeff_res[1,"CI_2.5"] <- as.numeric(nb_res$stan_summary[2,"2.5%"])
        nb_coeff_res[1,"CI_97.5"] <- as.numeric(nb_res$stan_summary[2,"97.5%"])
    }
    nb_coeff_res$cell_line <- cell_line
    nb_coeff_res$NMD_gene_excluded <- NMD_gene
    all_nb_coeff_res <- rbind(all_nb_coeff_res,nb_coeff_res)
}

# }

# Save
path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/Estimate_cNMDeff/excluding_genes/cNMDeff_",NMD_gene,".txt")
write.table(all_nb_coeff_res, file = path,
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
