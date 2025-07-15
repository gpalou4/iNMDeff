
reg_glm_model <- function(glm_char,TCGA_cancer, NMD_eff_merge, CNV_PCs, CNV_PCs_char) {

    if ( CNV_PCs != 0) {
        CNV_PCs <- CNV_PCs_char
    } 
    if (TCGA_cancer == "pancancer") {
        glm_model <- paste0("glm(",glm_char," + endogenous_purity + TMB + as.factor(cancer_subtype) + CNV_burden + as.factor(sex) + age +  sample_lib_size + ",CNV_PCs," , data = NMD_eff_merge, family = \"gaussian\", na.action = na.exclude)")
    } else {
        subtype_levels <- levels(as.factor(as.character((NMD_eff_merge$cancer_subtype))))
        if ( length(subtype_levels) == 0  || (TCGA_cancer == "TCGA-BRCA_basal") ) {
            glm_model <- paste0("glm(",glm_char," + endogenous_purity + ",CNV_PCs,", data = NMD_eff_merge , family = \"gaussian\", na.action = na.exclude)")
        } else {
            if (TCGA_cancer %in% c("TCGA-THYM","TCGA-DLBC")) {
                glm_model <- paste0("glm(",glm_char," + as.factor(cancer_subtype) + ",CNV_PCs," ,data = NMD_eff_merge, family = \"gaussian\", na.action = na.exclude)")
            } else {
                glm_model <- paste0("glm(",glm_char," + endogenous_purity + as.factor(cancer_subtype) + ",CNV_PCs," ,data = NMD_eff_merge, family = \"gaussian\", na.action = na.exclude)")
            }
        }
    }
    if ( CNV_PCs == 0) {
        glm_model <- gsub(" \\+ 0","",glm_model)
    }
    return(glm_model)
}

CNV_PCs_glm_char_function <- function(CNV_PCs) {
    PCs <- c()
    for (i in 1:as.numeric(CNV_PCs)) {
        PC <- paste0("Dim.",i)
        PCs <- c(PCs,PC)
    }
    CNV_PCs_glm_char <- paste(PCs,collapse=" + ")
    print(CNV_PCs_glm_char)
    return(CNV_PCs_glm_char)
}

# Open files

args <- commandArgs(trailingOnly=TRUE)
print(paste0("ARGUMENTS USED --> ",args))
TCGA_cancer <- args[1]
NMD_method <- args[2]
test <- args[3]
CNV_PCs <- as.numeric(args[4])
VAF <- args[5]
alpha <- as.numeric(args[6])
robust_SPCA <- args[7]
PCA_PCs <- as.numeric(args[8])
scale <- TRUE
center <- TRUE

# TCGA_cancer <- "TCGA-GBM" 
# TCGA_cancer <- "pancancer" 
# NMD_method <- "endogenous"
# test <- "CGC"
# CNV_PCs <- 1
# VAF <- 0.2
# alpha <- "3e-04"
# robust_SPCA <- "no"
# PCA_PCs <- 100

# 1) Data
# 1.1) CGC
CGC <- read.csv(file = "/g/strcombio/fsupek_cancer1/gpalou/COSMIC/cancer_gene_census_updated.tsv", 
                    header = TRUE, stringsAsFactors = FALSE)# 1.2) TCGA tissues
TCGA_tissues <- read.table("/g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/TCGA_projects_names.txt", 
                sep = "\t", header = FALSE, stringsAsFactors = FALSE)$V1

# 1.2) Number of optimal PCs
if (NMD_method == "ASE") {
    NMD_method_VAF <- paste0(NMD_method,"_",VAF)
} else {
    NMD_method_VAF <- NMD_method
}

best_PC <- read.table(file = paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/somatic_associations/optimizing_PCs/",NMD_method_VAF,"_best_PC_all_TCGA_cancers.txt"),
            header = TRUE, stringsAsFactors = TRUE, sep = "\t")

# 1.3) NMDeff and CGC somatic mut data
if (NMD_method == "ASE") {
    NMD_geneset <- paste0(NMD_method,"_stopgain_",VAF)
} else if (NMD_method == "endogenous") {
    NMD_geneset <- paste0(NMD_method,"_NMD_global_2_shared")
} else if (NMD_method == "PTC") {
    NMD_geneset <- paste0(NMD_method,"_stopgain_NMD_triggering")
}
# TCGA NMDeff
sample_NMD_efficiencies_TCGA_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/cancers/pancancer/NMD_efficiencies_TCGA.txt"
sample_NMD_efficiencies_TCGA <- read.table(file = sample_NMD_efficiencies_TCGA_path, header = TRUE, sep = "\t")
# Convert some columns to factors
factor_cols <- c("cancer_type","cancer_type_strat","cancer_subtype","LF_remove","purity_remove", "MSI_status",
                "batch_portion","batch_plate","batch_center","batch_vial","TCGA_full_barcode")
sample_NMD_efficiencies_TCGA[factor_cols] <- lapply(sample_NMD_efficiencies_TCGA[factor_cols], factor) 
cancers <- unique(sample_NMD_efficiencies_TCGA$cancer_type_strat)

if (TCGA_cancer == "pancancer") {
    print(TCGA_cancer)
    # CGC from cancer
    CGC_all_cancers_path <- paste0("/g/strcombio/fsupek_cancer1/gpalou/TCGA_somatic_variants/TCGA_CGC_somatic_mut_CNV.txt")
    CGC_all_cancers <- read.csv(file = CGC_all_cancers_path, sep ="\t", header = TRUE)
    # NMD efficiency from all cancers (merge)
    NMD_eff <- sample_NMD_efficiencies_TCGA
    NMD_eff$sample <- gsub("\\-","\\.",NMD_eff$sample)
    CGC_cancer <- CGC_all_cancers
    # Remove samples with <3 PTCs
    if (NMD_method %in% c("ASE")) {
        NMD_eff[which(NMD_eff[,paste0(NMD_method,"_num_PTCs_",VAF)] < 3),c(paste0(NMD_method,"_stopgain_",VAF))] <- NA
    }
    # PCA(somatic_CNV)
    tryCatch({
        TCGA_CNV_PCA <- read.table(file = paste0("/g/strcombio/fsupek_cancer1/gpalou/TCGA_CNV_PCA/",TCGA_cancer,"_sparse_PCA_ind_",alpha,"_robust_",robust_SPCA,"_num_PCs_",PCA_PCs,".txt"),
                                            header = TRUE, stringsAsFactors = FALSE, sep = "\t")
        # Remove PCs with 0
        cols <- colnames(TCGA_CNV_PCA)[which( colSums(TCGA_CNV_PCA) != 0 )]
        TCGA_CNV_PCA <- TCGA_CNV_PCA[,cols]
        # Scale PCs
        TCGA_CNV_PCA <- data.frame(scale(TCGA_CNV_PCA, scale = scale, center = center))
    }, error = function(e){
        print(e)
        }
    )
} else {
    # CGC from cancer
    TCGA_cancer_original <- strsplit(TCGA_cancer, "_")[[1]][1]
    print(TCGA_cancer_original)
    CGC_cancer_path <- paste0("/g/strcombio/fsupek_cancer1/gpalou/TCGA_somatic_variants/",gsub("TCGA-","",TCGA_cancer_original),"/",TCGA_cancer_original,"_CGC_somatic_mut_CNV.txt")
    CGC_cancer <- read.csv(file = CGC_cancer_path, sep ="\t", header = TRUE, stringsAsFactors = FALSE)
    # NMD efficiency from cancer
    NMD_eff <- sample_NMD_efficiencies_TCGA[which(sample_NMD_efficiencies_TCGA$cancer_type %in% TCGA_cancer),]
    NMD_eff$sample <- gsub("\\-","\\.",NMD_eff$sample)
    # Remove samples with <3 PTCs
    if (NMD_method %in% c("ASE")) {
        NMD_eff[which(NMD_eff[,paste0(NMD_method,"_num_PTCs_",VAF)] < 3),c(paste0(NMD_method,"_stopgain_",VAF))] <- NA
    }
    tryCatch({
        #TCGA_CNV_PCA <- read.table(file = paste0("/g/strcombio/fsupek_cancer1/gpalou/TCGA_CNV_PCA/",TCGA_cancer_original,"_sparse_PCA_ind_",alpha,"_robust_",robust_SPCA,".txt"),
        #                                    header = TRUE, stringsAsFactors = FALSE, sep = "\t")
        TCGA_CNV_PCA <- read.table(file = paste0("/g/strcombio/fsupek_cancer1/gpalou/TCGA_CNV_PCA/pancancer_sparse_PCA_ind_",alpha,"_robust_",robust_SPCA,"_num_PCs_",PCA_PCs,".txt"),
                                            header = TRUE, stringsAsFactors = FALSE, sep = "\t")
        # Remove samples that are not found in our cancer type
        TCGA_CNV_PCA <- TCGA_CNV_PCA[rownames(TCGA_CNV_PCA) %in% NMD_eff$sample,]
        # Remove PCs with 0
        cols <- colnames(TCGA_CNV_PCA)[which( colSums(TCGA_CNV_PCA) != 0 )]
        TCGA_CNV_PCA <- TCGA_CNV_PCA[,cols]
        # Scale PCs
        TCGA_CNV_PCA <- data.frame(scale(TCGA_CNV_PCA, scale = scale, center = center))
    }, error = function(e){
        print(e)
        }
    )
}

sample_size <- nrow(TCGA_CNV_PCA)
nPCs <- ncol(TCGA_CNV_PCA)

# Number of PCs for each mutation
row <- best_PC$TCGA_cancer %in% TCGA_cancer
if (sample_size <= 80) {
    CNV_PCs_ob <- best_PC[row,"best_PC_som_CNV"]
    if (!is.na(CNV_PCs_ob)) {
        CNV_PCs_glm_char_som_CNV <- CNV_PCs_glm_char_function(CNV_PCs_ob)
    } else {
        CNV_PCs_glm_char_som_CNV <- NULL
    }
} else {
    CNV_PCs_ob <- nPCs
    CNV_PCs_glm_char_som_CNV <- paste0(colnames(TCGA_CNV_PCA), collapse =" + ")
    print(CNV_PCs_glm_char_som_CNV)
}

# if (sample_size <= 80) {
#     CNV_PCs_ob <- best_PC[row,"best_PC_som_mut_miss_tr"]
#     if (!is.na(CNV_PCs_ob)) {
#         CNV_PCs_glm_char_som_mut_miss_tr <- CNV_PCs_glm_char_function(1)
#     } else {
#         CNV_PCs_glm_char_som_mut_miss_tr <- NULL
#     } 
# } else {
#     CNV_PCs_ob <- nPCs
# }
CNV_PCs_glm_char_som_mut_miss_tr <- ""

# if (sample_size <= 80) {
#     CNV_PCs_ob <- best_PC[row,"best_PC_som_mut_synonymous"]
#     if (!is.na(CNV_PCs_ob)) {
#         CNV_PCs_glm_char_som_mut_synonymous <- CNV_PCs_glm_char_function(1)
#     } else {
#         CNV_PCs_glm_char_som_mut_synonymous <- NULL
#     }
# } else {
#     CNV_PCs_ob <- nPCs
# }
CNV_PCs_glm_char_som_mut_synonymous <- ""

# 2) Regressions

# 2.1) per gene

if (test == "CGC") {

    df_ass <- CGC_cancer[,c("Gene.Symbol","Synonyms")]
    colnames(df_ass) <- c("Gene_symbol","ENSEMBL_gene")
    cols <- c("som_mut_truncating_coeff","som_mut_truncating_pval","som_mut_missense_coeff","som_mut_missense_pval", "som_mut_synonymous_coeff",
                "som_mut_synonymous_pval","som_CNV_amp_coeff","som_CNV_amp_pval","som_CNV_del_coeff","som_CNV_del_pval","error")
    df_ass[,cols] <- NA
    number_comparisons_to_reduce <- 0

    for (ENSEMBL_gene in as.character(unique(df_ass$ENSEMBL_gene))) {

        if (is.na(ENSEMBL_gene)) {
            print("No ENSEMBL gene ID ... skipping")
            next
        }

        #ENSEMBL_gene <- df_ass[df_ass$Gene_symbol%in%"SMG5","ENSEMBL_gene"]

        row_num <- which(df_ass[,"ENSEMBL_gene"] %in% ENSEMBL_gene)
        print(ENSEMBL_gene)
        CGC_cancer_gene <- CGC_cancer[CGC_cancer$Synonyms %in% ENSEMBL_gene,]
        # 2.1) For somatic truncating mutations
        som_mut_index <- grep("somatic_mut_truncating",colnames(CGC_cancer_gene))
        CGC_cancer_gene_som_mut <- CGC_cancer_gene[,som_mut_index, drop = FALSE]
        # Create data frame
        reg_df_som_mut_truncating <- data.frame(t(CGC_cancer_gene_som_mut))
        colnames(reg_df_som_mut_truncating) <- "somatic_mut_truncating"
        reg_df_som_mut_truncating$TCGA_sample <- gsub("(TCGA.*)_somatic_mut_truncating","\\1",rownames(reg_df_som_mut_truncating))
        # 2.2) For somatic missense mutations
        som_mut_index <- grep("somatic_mut_missense",colnames(CGC_cancer_gene))
        CGC_cancer_gene_som_mut <- CGC_cancer_gene[,som_mut_index, drop = FALSE]
        # Create data frame
        reg_df_som_mut_missense <- data.frame(t(CGC_cancer_gene_som_mut))
        colnames(reg_df_som_mut_missense) <- "somatic_mut_missense"
        reg_df_som_mut_missense$TCGA_sample <- gsub("(TCGA.*)_somatic_mut_missense","\\1",rownames(reg_df_som_mut_missense))
        # 2.3) For somatic synonymous mutations
        som_mut_index <- grep("somatic_mut_synonymous",colnames(CGC_cancer_gene))
        CGC_cancer_gene_som_mut <- CGC_cancer_gene[,som_mut_index, drop = FALSE]
        # Create data frame
        reg_df_som_mut_synonymous <- data.frame(t(CGC_cancer_gene_som_mut))
        colnames(reg_df_som_mut_synonymous) <- "somatic_mut_synonymous"
        reg_df_som_mut_synonymous$TCGA_sample <- gsub("(TCGA.*)_somatic_mut_synonymous","\\1",rownames(reg_df_som_mut_synonymous))
        # 2.4) For CNV
        som_CNV_index <- grep("somatic_CNV",colnames(CGC_cancer_gene))
        CGC_cancer_gene_som_CNV <- CGC_cancer_gene[,som_CNV_index, drop = FALSE]
        # Create data frame
        reg_df_som_CNV <- data.frame(t(CGC_cancer_gene_som_CNV))
        colnames(reg_df_som_CNV) <- "somatic_CNV"
        reg_df_som_CNV$TCGA_sample <- gsub("(TCGA.*)_somatic_CNV","\\1",rownames(reg_df_som_CNV))
        # 2.5) Merge and add NMDeff
        NMD_eff_merge <- merge(NMD_eff,reg_df_som_CNV, by.x = "sample", by.y = "TCGA_sample", all.x = TRUE)
        NMD_eff_merge <- merge(NMD_eff_merge,reg_df_som_mut_truncating, by.x = "sample", by.y = "TCGA_sample", all.x = TRUE)
        NMD_eff_merge <- merge(NMD_eff_merge,reg_df_som_mut_missense, by.x = "sample", by.y = "TCGA_sample", all.x = TRUE)
        NMD_eff_merge <- merge(NMD_eff_merge,reg_df_som_mut_synonymous, by.x = "sample", by.y = "TCGA_sample", all.x = TRUE)
        # 2.6) Add PCA CNV
        NMD_eff_merge <- merge(NMD_eff_merge,TCGA_CNV_PCA, by.x = "sample", by.y = "row.names", all.x = TRUE)

        # 3) Regressions for somatic mutations
        # Remove overlapping somatic CNV
        NMD_eff_merge_filt <- NMD_eff_merge[which(NMD_eff_merge$somatic_CNV == "no"),]
        # Normalize TMB
        NMD_eff_merge_filt$TMB <- sqrt(NMD_eff_merge_filt$TMB)
        # Scale all variables
        numeric_cols <- unlist(lapply(NMD_eff_merge_filt, is.numeric), use.names = FALSE)
        numeric_cols[1:30] <- FALSE
        NMD_eff_merge_filt[numeric_cols] <- lapply(NMD_eff_merge_filt[numeric_cols], function(x) c(scale(x)))

        # Check sample size
        som_mut_truncating_levels <- unique(na.omit(NMD_eff_merge_filt$somatic_mut_truncating))
        som_mut_missense_levels <- unique(na.omit(NMD_eff_merge_filt$somatic_mut_missense))
        som_mut_synonymous_levels <- unique(na.omit(NMD_eff_merge_filt$somatic_mut_synonymous))
        som_mut_truncating_size <- sum(na.omit((NMD_eff_merge_filt$somatic_mut_truncating == "yes")))
        som_mut_missense_size <- sum(na.omit((NMD_eff_merge_filt$somatic_mut_missense == "yes")))
        som_mut_synonymous_size <- sum(na.omit((NMD_eff_merge_filt$somatic_mut_synonymous == "yes")))

        n <- 2
        # 3.1) Truncating and Missense mut
        somMutReg <- TRUE
        if ( som_mut_truncating_size >= n & som_mut_missense_size == 0 ) {
            glm_char <- paste0(NMD_geneset," ~ somatic_mut_truncating")
        } else if ( som_mut_truncating_size == 0 & som_mut_missense_size >= n ) { # It has only missense mut
            glm_char <-  paste0(NMD_geneset," ~ somatic_mut_missense")
        } else if ( som_mut_truncating_size >=n & som_mut_missense_size >= n ) { # It has both
            glm_char <-  paste0(NMD_geneset," ~ somatic_mut_truncating + somatic_mut_missense")
        } else if ( som_mut_truncating_size == 0 & som_mut_missense_size == 0 ) { # It has no som mut
            print("No somatic mutations in this gene, skipping regression")
            number_comparisons_to_reduce <- number_comparisons_to_reduce + 1
            somMutReg <- FALSE
        } else {
            print("No enough sample size in this gene, skipping regression")
            number_comparisons_to_reduce <- number_comparisons_to_reduce + 1
            somMutReg <- FALSE
        }

        # Regression
        if (isTRUE(somMutReg) && !is.null(CNV_PCs_glm_char_som_mut_miss_tr)) {
            glm_model <- reg_glm_model(glm_char = glm_char, TCGA_cancer = TCGA_cancer, NMD_eff_merge = NMD_eff_merge_filt, CNV_PCs = CNV_PCs, CNV_PCs_char = CNV_PCs_glm_char_som_mut_miss_tr)
            glm_model <- gsub("\\+  ,",",",glm_model)
            tryCatch({
                glm_res <- eval(parse(text=glm_model))
                glm_res <- summary(glm_res)
                if ( (TCGA_cancer != "pancancer") & (CNV_PCs == 2) ) {
                    # Update number of CNV PCs used for the regression
                    CNV_PCs_updated <- as.numeric(gsub("Dim.","",rownames(glm_res$coefficients))[length(rownames(glm_res$coefficients))])
                    PCs <- c()
                    for (i in 1:as.numeric(CNV_PCs_updated)) {
                        PC <- paste0("Dim.",i)
                        PCs <- c(PCs,PC)
                    }
                    CNV_PCs_glm_char_tmp <- paste(PCs,collapse=" + ")
                    # Re-do the regression
                    glm_model <- reg_glm_model(glm_char = glm_char, TCGA_cancer = TCGA_cancer, NMD_eff_merge = NMD_eff_merge_filt, CNV_PCs = CNV_PCs, CNV_PCs_char = CNV_PCs_glm_char_tmp)
                    glm_res <- eval(parse(text=glm_model))
                    glm_res <- summary(glm_res)
                    # If still coefficients have no p-value, reduce PCs one by one
                    if (!is.na(glm_res$coefficients[2,2])) {
                        PCop <- TRUE
                    } else {
                        PCop <- FALSE
                        PC <- CNV_PCs_updated
                        while ( !isTRUE(PCop) ) {
                            PC <- PC - 1 
                            CNV_PCs_updated_tmp <- PC
                            PCs <- c()
                            for (i in 1:as.numeric(CNV_PCs_updated_tmp)) {
                                PC <- paste0("Dim.",i)
                                PCs <- c(PCs,PC)
                            }
                            CNV_PCs_glm_char_tmp <- paste(PCs,collapse=" + ")
                            glm_model <- reg_glm_model(glm_char = glm_char, TCGA_cancer = TCGA_cancer, NMD_eff_merge = NMD_eff_merge_filt, CNV_PCs = CNV_PCs, CNV_PCs_char = CNV_PCs_glm_char_tmp)
                            glm_res <- eval(parse(text=glm_model))
                            glm_res <- summary(glm_res)
                            if (!is.na(glm_res$coefficients[1,2])) {
                                PCop <<- TRUE
                            }
                        }
                    }
                }
                hasMutMissReg <- grep("missense",rownames(glm_res$coefficients))
                hasMutTruncReg <- grep("truncating",rownames(glm_res$coefficients))
                if (length(hasMutMissReg) != 0) {
                    df_ass[row_num,"som_mut_missense_coeff"] <- glm_res$coefficients[hasMutMissReg,"Estimate"]
                    df_ass[row_num,"som_mut_missense_pval"] <- glm_res$coefficients[hasMutMissReg,"Pr(>|t|)"]
                } 
                if (length(hasMutTruncReg) != 0) {
                    df_ass[row_num,"som_mut_truncating_coeff"] <- glm_res$coefficients[hasMutTruncReg,"Estimate"]
                    df_ass[row_num,"som_mut_truncating_pval"] <- glm_res$coefficients[hasMutTruncReg,"Pr(>|t|)"]
                } 
            }, error = function(e){
                df_ass[row_num,"error"] <<- "regression error, probably contrasts"
                print(e)
                }
            )
        }
        # 3.2) Synonymous
        # It has no Syn mut
        if ( som_mut_synonymous_size <= 1 ) { 
            print("No somatic synonymous mut in this gene, skipping regression")
            number_comparisons_to_reduce <- number_comparisons_to_reduce + 1
            somMutSynReg <- FALSE
        } else {
            glm_char <-  paste0(NMD_geneset," ~ somatic_mut_synonymous")
            somMutSynReg <- TRUE
        }
        # Regression
        if (isTRUE(somMutSynReg) && !is.null(CNV_PCs_glm_char_som_mut_synonymous)) {
            glm_model <- reg_glm_model(glm_char = glm_char, TCGA_cancer = TCGA_cancer, NMD_eff_merge = NMD_eff_merge_filt, CNV_PCs = CNV_PCs, CNV_PCs_char = CNV_PCs_glm_char_som_mut_synonymous)
            glm_model <- gsub("\\+  ,",",",glm_model)
            tryCatch({
                glm_res <- eval(parse(text=glm_model))
                glm_res <- summary(glm_res)
                if ( (TCGA_cancer != "pancancer") & (CNV_PCs == 2) ) {
                    # Update number of CNV PCs used for the regression
                    CNV_PCs_updated <- as.numeric(gsub("Dim.","",rownames(glm_res$coefficients))[length(rownames(glm_res$coefficients))])
                    PCs <- c()
                    for (i in 1:as.numeric(CNV_PCs_updated)) {
                        PC <- paste0("Dim.",i)
                        PCs <- c(PCs,PC)
                    }
                    CNV_PCs_glm_char_tmp <- paste(PCs,collapse=" + ")
                    # Re-do the regression
                    glm_model <- reg_glm_model(glm_char = glm_char, TCGA_cancer = TCGA_cancer, NMD_eff_merge = NMD_eff_merge_filt, CNV_PCs = CNV_PCs, CNV_PCs_char = CNV_PCs_glm_char_tmp)
                    glm_res <- eval(parse(text=glm_model))
                    glm_res <- summary(glm_res)
                    # If still coefficients have no p-value, reduce PCs one by one
                    if (!is.na(glm_res$coefficients[2,2])) {
                        PCop <- TRUE
                    } else {
                        PCop <- FALSE
                        PC <- CNV_PCs_updated
                        while ( !isTRUE(PCop) ) {
                            PC <- PC - 1 
                            CNV_PCs_updated_tmp <- PC
                            PCs <- c()
                            for (i in 1:as.numeric(CNV_PCs_updated_tmp)) {
                                PC <- paste0("Dim.",i)
                                PCs <- c(PCs,PC)
                            }
                            CNV_PCs_glm_char_tmp <- paste(PCs,collapse=" + ")
                            glm_model <- reg_glm_model(glm_char = glm_char, TCGA_cancer = TCGA_cancer, NMD_eff_merge = NMD_eff_merge_filt, CNV_PCs = CNV_PCs, CNV_PCs_char = CNV_PCs_glm_char_tmp)
                            glm_res <- eval(parse(text=glm_model))
                            glm_res <- summary(glm_res)
                            if (!is.na(glm_res$coefficients[1,2])) {
                                PCop <<- TRUE
                            }
                        }
                    }
                }
                hasMutSynReg <- grep("synonymous",rownames(glm_res$coefficients))
                if (length(hasMutSynReg) != 0) {
                    df_ass[row_num,"som_mut_synonymous_coeff"] <- glm_res$coefficients[hasMutSynReg,"Estimate"]
                    df_ass[row_num,"som_mut_synonymous_pval"] <- glm_res$coefficients[hasMutSynReg,"Pr(>|t|)"]
                } 
            }, error = function(e){
                df_ass[row_num,"error"] <<- "regression error, probably contrasts"
                print(e)
                }
            )
        }
        # 4) Regressions for somatic CNV
        # Remove overlapping somatic truncating/missense mut
        NMD_eff_merge_filt <- NMD_eff_merge[which(NMD_eff_merge$somatic_mut_missense == "no" & NMD_eff_merge$somatic_mut_truncating == "no"),]
        # Normalize TMB
        NMD_eff_merge_filt$TMB <- sqrt(NMD_eff_merge_filt$TMB)
        # Scale all variables
        numeric_cols <- unlist(lapply(NMD_eff_merge_filt, is.numeric), use.names = FALSE)
        numeric_cols[1:30] <- FALSE
        NMD_eff_merge_filt[numeric_cols] <- lapply(NMD_eff_merge_filt[numeric_cols], function(x) c(scale(x)))
        # aggregate(NMD.global.2.shared ~ somatic_CNV , data = NMD_eff_merge, median)
        som_mut_CNV_levels <- unique(na.omit(NMD_eff_merge_filt$somatic_CNV))
        # It has no CNV mut
        if ( length(som_mut_CNV_levels) == 1 ) { 
            print("No somatic CNV in this gene, skipping regression")
            number_comparisons_to_reduce <- number_comparisons_to_reduce + 1
            somCNVReg <- FALSE
        } else {
            glm_char <- paste0(NMD_geneset," ~ relevel(as.factor(somatic_CNV), ref = \"no\")")
            somCNVReg <- TRUE
        }
        # Regression
        if (isTRUE(somCNVReg) && !is.null(CNV_PCs_glm_char_som_CNV)) {
            glm_model <- reg_glm_model(glm_char = glm_char, TCGA_cancer = TCGA_cancer, NMD_eff_merge = NMD_eff_merge_filt, CNV_PCs = CNV_PCs, CNV_PCs_char = CNV_PCs_glm_char_som_CNV)
            tryCatch({
                glm_res <- eval(parse(text=glm_model))
                glm_res <- summary(glm_res)
                if ( (TCGA_cancer != "pancancer") & (CNV_PCs == 1) ) {
                    # Update number of CNV PCs used for the regression
                    #CNV_PCs_updated <- as.numeric(gsub("Dim.","",rownames(glm_res$coefficients))[length(rownames(glm_res$coefficients))])
                    # PCs <- c()
                    # for (i in 1:as.numeric(CNV_PCs_updated)) {
                    #     PC <- paste0("Dim.",i)
                    #     PCs <- c(PCs,PC)
                    # }
                    PCs <- rownames(glm_res$coefficients)[grep("Dim.",rownames(glm_res$coefficients))]
                    CNV_PCs_glm_char_tmp <- paste(PCs,collapse=" + ")
                    # Re-do the regression
                    glm_model <- reg_glm_model(glm_char = glm_char, TCGA_cancer = TCGA_cancer, NMD_eff_merge = NMD_eff_merge_filt, CNV_PCs = CNV_PCs, CNV_PCs_char = CNV_PCs_glm_char_tmp)
                    glm_res <- eval(parse(text=glm_model))
                    glm_res <- summary(glm_res)
                    # If still coefficients have no p-value, reduce PCs one by one
                    if (!is.na(glm_res$coefficients[2,2])) {
                        PCop <- TRUE
                    } else {
                        PCop <- FALSE
                        PC <- length(PCs)
                        while ( !isTRUE(PCop) ) {
                            PC <- PC - 1 
                            CNV_PCs_glm_char_tmp <- paste(PCs[1:PC],collapse=" + ")
                            # CNV_PCs_updated_tmp <- PC
                            # PCs <- c()
                            # for (i in 1:as.numeric(CNV_PCs_updated_tmp)) {
                            #     PC <- paste0("Dim.",i)
                            #     PCs <- c(PCs,PC)
                            # }
                            # CNV_PCs_glm_char_tmp <- paste(PCs,collapse=" + ")
                            glm_model <- reg_glm_model(glm_char = glm_char, TCGA_cancer = TCGA_cancer, NMD_eff_merge = NMD_eff_merge_filt, CNV_PCs = CNV_PCs, CNV_PCs_char = CNV_PCs_glm_char_tmp)
                            glm_res <- eval(parse(text=glm_model))
                            glm_res <- summary(glm_res)
                            if (!is.na(glm_res$coefficients[1,2])) {
                                PCop <<- TRUE
                            }
                        }
                    }
                }
                hasCNVampReg <- grep("\\)amp",rownames(glm_res$coefficients))
                hasCNVdelReg <- grep("\\)del",rownames(glm_res$coefficients))
                if (length(hasCNVampReg) != 0) {
                    df_ass[row_num,"som_CNV_amp_coeff"] <- glm_res$coefficients[hasCNVampReg,"Estimate"]
                    df_ass[row_num,"som_CNV_amp_pval"] <- glm_res$coefficients[hasCNVampReg,"Pr(>|t|)"]
                } 
                if (length(hasCNVdelReg) != 0) {
                    df_ass[row_num,"som_CNV_del_coeff"] <- glm_res$coefficients[hasCNVdelReg,"Estimate"]
                    df_ass[row_num,"som_CNV_del_pval"] <- glm_res$coefficients[hasCNVdelReg,"Pr(>|t|)"]
                } 
            }, error = function(e){
                df_ass[row_num,"error"] <<- "regression error, probably contrasts"
                print(e)
                }
            )
        }
    }
    # 5) Write results
    output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/somatic_associations/final_associations/",TCGA_cancer,"_",NMD_method_VAF,"_CGC_somatic_mut_CNV_PCs_",CNV_PCs,".txt")
    write.table(df_ass, file = output_path, sep = "\t", quote = FALSE,
                        col.names = TRUE, row.names = FALSE)
} else if (test == "subtypes") { # 4.2) per RNAseq NMF cluster subtype

    subtype_num <- max(na.omit(NMD_eff$subtype))
    df_ass <- data.frame(subtype = 1:subtype_num)
    cols <- c("subtype_coeff","subtype_pval","error")
    df_ass[,cols] <- NA
    number_comparisons_to_reduce <- 0
    # Check if cancer has subtypes
    if ( sum(!is.na(NMD_eff$subtype)) == 0 ) {
        print("This cancer has no subtypes, skipping")
        subtypeReg <- FALSE
    } else{
        subtypeReg <- TRUE
    }
    if (isTRUE(subtypeReg)) {
            # Regression
            glm_model <- paste0("glm(",NMD_geneset," ~ as.factor(subtype) + endogenous_purity,data = NMD_eff, family = \"gaussian\", na.action = na.exclude)")
            tryCatch({
                glm_res <- eval(parse(text=glm_model))
                glm_res <- summary(glm_res)
                for (i in 2:subtype_num) {
                    subtype <- which(df_ass[,"subtype"] %in% i)
                    print(paste0("RNAseq NMF cluster subtype --> ",subtype))
                    hasSubtypeReg <- grep(subtype,rownames(glm_res$coefficients))
                    if (length(hasSubtypeReg) != 0) {
                        df_ass[hasSubtypeReg,"subtype_coeff"] <- glm_res$coefficients[hasSubtypeReg,"Estimate"]
                        df_ass[hasSubtypeReg,"subtype_pval"] <- glm_res$coefficients[hasSubtypeReg,"Pr(>|t|)"]
                    } 
                }}, error = function(e){
                        df_ass[hasSubtypeReg,"error"] <<- "regression error, probably contrasts"
                        print(e)
                }
            )
            }
        # Write results
        output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/subtypes/",TCGA_cancer,"_",NMD_method_VAF,"_","subtypes.txt")
        write.table(df_ass, file = output_path, sep = "\t", quote = FALSE,
                            col.names = TRUE, row.names = FALSE)
}
