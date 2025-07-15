library("ggplot2")
library("caret")
library("glmnet")
library("preprocessCore")
library("sva")
library("viridis")
library("ggpubr")

# Compute R^2 from true and predicted values
eval_results <- function(true, predicted, df) {
    SSE <- sum((predicted - true)^2)
    SST <- sum((true - mean(true, na.rm = TRUE))^2)
    R_square <- 1 - SSE / SST
    RMSE <- sqrt(SSE / nrow(df))
    # Model performance metrics
    data.frame(
        RMSE = RMSE,
        Rsquare = R_square
    )
}

# FUNCION TO QN + COMBAT
qn <- function(dfnorm, dataset) {
    tar <- normalize.quantiles.determine.target(x = t(as.matrix(dfnorm[dfnorm$type == dataset, !colnames(dfnorm) %in% c("type")]))) # only tumor samples for target
    dfnorm[, !colnames(dfnorm) %in% c("type")] <- t(normalize.quantiles.use.target(x = t(as.matrix(dfnorm[, !colnames(dfnorm) %in% c("type")])), target = tar, copy = FALSE)) # apply QN to tumors+cell lines
    # Combat
    cleandat <- ComBat(
        dat = t(dfnorm[, !colnames(dfnorm) %in% c("type")]), batch = dfnorm$type,
        mod = NULL, # ref.batch = dataset,
        par.prior = TRUE
    )
    cleandat[1:5, 1:5]
    ma_fi <- data.frame(t(cleandat))
    ma_fi$type <- dfnorm$type
    dim(ma_fi)
    return(ma_fi)
}

# 1) Data
# 1.1) Obtain the Validation gene expression Cell Line datasets
# Girish et al 2022 1q Loss
input_path <- "/g/strcombio/fsupek_cancer1/gpalou/papers/Girish_1qloss/"
cell_line_files <- list.files(input_path)
cell_lines <- c("A2058", "HCT116", "AGS", "A2780")
CL_list <- list()
for (cell_line in cell_lines) {
    cell_line_files_filt <- cell_line_files[grep(cell_line, cell_line_files)]
    all_CL_df <- c()
    for (CL_file_char in cell_line_files_filt) {
        print(CL_file_char)
        CL_file <- read.table(file = paste0(input_path, CL_file_char), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
        CL_file <- CL_file[, c(1:2)]
        colnames(CL_file)[1] <- "gene_id"
        cols <- grep("1q\\.loss|WT", colnames(CL_file))
        if (length(cols) == 0) {
            next
        }
        print(colnames(CL_file))
        if (length(all_CL_df) == 0) {
            all_CL_df <- CL_file
        } else {
            all_CL_df <- merge(all_CL_df, CL_file, by = "gene_id", all = TRUE)
        }
    }
    colnames(all_CL_df)[-1] <- paste0(cell_line, "_", colnames(all_CL_df)[-1])
    CL_list[[cell_line]] <- all_CL_df
}

# DepMap gene expression from ~1k CLs
output_path <- "/g/strcombio/fsupek_cancer1/gpalou/DepMap/OmicsExpressionProteinCodingGenesTPMLogp1.csv"
RNAseq_DepMap_TPM_CL <- read.csv(file = output_path, header = TRUE)
colnames(RNAseq_DepMap_TPM_CL)[1] <- "DepMap_ID"

# 1.2) Classification of CL based on arm-level CNV by --> Cohen Sharir Nat 2021
output_path <- "/g/strcombio/fsupek_cancer1/gpalou/DepMap/Cohen_Sharir_Nat_2021_DepMap_CNV_CL.csv"
DepMap_CNV_arm_level_type <- read.csv(file = output_path, header = TRUE, sep = ";")

# 1.3) Gene-level TPM RNAseq data
# TCGA
output_path <- "/g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/primary_tumor/TCGA_RNAseq_matrix_TPM_gene.txt"
# write.table(RNAseq_TCGA_TPM_all, file = output_path, sep = "\t", quote = FALSE,
#             col.names = TRUE, row.names = TRUE)
RNAseq_TCGA_TPM_all <- read.table(file = output_path, header = TRUE, sep = "\t", row.names = 1)
print("Dimensions -->")
print(dim(RNAseq_TCGA_TPM_all))
# Check Samples with NAs
na_counts <- colSums(is.na(RNAseq_TCGA_TPM_all))
cols_na <- names(na_counts[na_counts > 0])
RNAseq_TCGA_TPM_all <- RNAseq_TCGA_TPM_all[, !colnames(RNAseq_TCGA_TPM_all) %in% cols_na]
# Check Genes with NAs
na_counts <- rowSums(is.na(RNAseq_TCGA_TPM_all))
rows_na <- names(na_counts[na_counts > 0])
print("Dimensions -->")
print(dim(RNAseq_TCGA_TPM_all))

# GTEx
GTEx_path <- "/g/strcombio/fsupek_cancer1/gpalou/GTEx/RNAseq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct"
RNAseq_GTEx_TPM_all <- read.table(file = GTEx_path, header = TRUE, sep = "\t", row.names = 1, skip = 2)
print("Dimensions -->")
print(dim(RNAseq_GTEx_TPM_all))
RNAseq_GTEx_TPM_all <- RNAseq_GTEx_TPM_all[-grep("PAR", rownames(RNAseq_GTEx_TPM_all)), ]
rownames(RNAseq_GTEx_TPM_all) <- gsub("(.*)\\..*", "\\1", rownames(RNAseq_GTEx_TPM_all))
RNAseq_GTEx_TPM_all$Description <- NULL

# 1.4) samples NMDeff pantissue
# TCGA
NMD_genesets <- c("endogenous_NMD_global_2_shared", "ASE_stopgain_0.01", "ASE_stopgain_0.2", "PTCs_stopgain_NMD_triggering")
# PTC // ASE // Endogenous
sample_NMD_efficiencies_TCGA_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/cancers/pancancer/NMD_efficiencies_TCGA.txt")
sample_NMD_efficiencies_TCGA <- read.table(file = sample_NMD_efficiencies_TCGA_path, header = TRUE, sep = "\t")
# Convert some columns to factors
factor_cols <- c(
    "cancer_type", "cancer_type_strat", "cancer_subtype", "LF_remove", "purity_remove", "MSI_status",
    "batch_portion", "batch_plate", "batch_center", "batch_vial", "TCGA_full_barcode"
)
sample_NMD_efficiencies_TCGA[factor_cols] <- lapply(sample_NMD_efficiencies_TCGA[factor_cols], factor)
cancers <- unique(sample_NMD_efficiencies_TCGA$cancer_type_strat)
# Filters for TCGA
sample_NMD_efficiencies_TCGA[which(sample_NMD_efficiencies_TCGA$ASE_num_PTCs_0.2 < 3), c("ASE_stopgain_0.2")] <- NA
sample_NMD_efficiencies_TCGA[which(sample_NMD_efficiencies_TCGA$ASE_num_PTCs_0.01 < 3), c("ASE_stopgain_0.01")] <- NA

# GTEx
sample_NMD_efficiencies_GTEx_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/GTEx/pantissue/NMD_efficiencies_GTEx.txt"
sample_NMD_efficiencies_GTEx <- read.table(file = sample_NMD_efficiencies_GTEx_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# Convert some columns to factors
factor_cols <- c("tissue", "sample")
sample_NMD_efficiencies_GTEx[factor_cols] <- lapply(sample_NMD_efficiencies_GTEx[factor_cols], factor)
GTEx_tissues <- unique(sample_NMD_efficiencies_GTEx$acronyms)
# Filters for GTEx
sample_NMD_efficiencies_GTEx[which(sample_NMD_efficiencies_GTEx$ASE_num_PTCs_0.2 < 3), c("ASE_stopgain_0.2")] <- NA
sample_NMD_efficiencies_GTEx[which(sample_NMD_efficiencies_GTEx$ASE_num_PTCs_0.01 < 3), c("ASE_stopgain_0.01")] <- NA

# 1.5) NMD target genes
NMD_genes <- read.table(
    file = "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/NMD_targets/ENSEMBL_transcripts/NMD_global_2_shared_ensembl.txt",
    header = TRUE, sep = "\t"
)
NMD_ensembl_genes <- as.character(unique(NMD_genes$ensembl_gene_id))

# 1.6) CNV-PCs pancancer
scale <- TRUE
center <- TRUE
alpha <- "3e-04"
num_PCs <- "100"
tryCatch(
    {
        TCGA_CNV_PCA_ind <- read.table(
            file = paste0("/g/strcombio/fsupek_cancer1/gpalou/TCGA_CNV_PCA/pancancer_sparse_PCA_ind_", alpha, "_robust_no_num_PCs_", num_PCs, ".txt"),
            header = TRUE, stringsAsFactors = FALSE, sep = "\t"
        )
        rownames(TCGA_CNV_PCA_ind) <- gsub("\\.", "-", rownames(TCGA_CNV_PCA_ind))
        # Remove PCs with 0
        cols <- colnames(TCGA_CNV_PCA_ind)[which(colSums(TCGA_CNV_PCA_ind) != 0)]
        TCGA_CNV_PCA_ind <- TCGA_CNV_PCA_ind[, cols]
        print("PCA dimensions --> ")
        print(dim(TCGA_CNV_PCA_ind))
        nPCs <- ncol(TCGA_CNV_PCA_ind)
        # Scale PCs
        TCGA_CNV_PCA_ind <- data.frame(scale(TCGA_CNV_PCA_ind, scale = scale, center = center))
    },
    error = function(e) {
        print(e)
    }
)
# Merge
sample_NMD_efficiencies_TCGA <- merge(sample_NMD_efficiencies_TCGA, TCGA_CNV_PCA_ind, by.x = "sample", by.y = "row.names", all.x = TRUE)

# 1.7) ENSEMBL gene ID - Gene Symbol conversor table
ensembl_gene_symbol <- read.table(file = "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/ensembl_v88_gene_transcript_genesymbol.txt",
                                    header = TRUE, sep = "\t", stringsAsFactors = FALSE)
ensembl_gene_symbol <- unique(ensembl_gene_symbol[,c(1:2)])
ensembl_gene_symbol <- ensembl_gene_symbol[!duplicated(ensembl_gene_symbol$gene_name),]

# 2) Shared predictors (genes) across the datasets (discovery and validation datasets)
#CL_genes <- intersect(intersect(CL_list[[1]]$gene_id, CL_list[[3]]$gene_id), CL_list[[4]]$gene_id)
# Convert gene symbols to ensembl gene id in CL DepMap
colnames(RNAseq_DepMap_TPM_CL) <- gsub("(.*)\\.\\..*","\\1",colnames(RNAseq_DepMap_TPM_CL))
gene_symbols <- data.frame(gene_symbol = colnames(RNAseq_DepMap_TPM_CL))
df <- merge(gene_symbols, ensembl_gene_symbol, by.x = "gene_symbol", by.y = "gene_name", all.x = TRUE, sort = FALSE)
df <- na.omit(df)
colnames(RNAseq_DepMap_TPM_CL)[match(df$gene_symbol,colnames(RNAseq_DepMap_TPM_CL))] <- as.character(df$gene_id)
rownames(RNAseq_DepMap_TPM_CL) <- RNAseq_DepMap_TPM_CL$DepMap_ID
RNAseq_DepMap_TPM_CL$DepMap_ID <- NULL
RNAseq_DepMap_TPM_CL_t <- t(RNAseq_DepMap_TPM_CL)
# Shared genes
CL_genes <- colnames(RNAseq_DepMap_TPM_CL)
TCGA_genes <- rownames(RNAseq_TCGA_TPM_all)
GTEx_genes <- rownames(RNAseq_GTEx_TPM_all)
all_shared_genes <- intersect(intersect(CL_genes, TCGA_genes), GTEx_genes)
datasets <- c("TCGA", "GTEx")
models <- c("LASSO", "ridge")
NMDeff_methods <- c("ASE_stopgain_0.2", "endogenous_NMD_global_2_shared")

# 3) Create the sample-NMDeff ridge model using gene-level expression data

NMDeff_gene_model_validation <- function(dataset, model, target_var, genes_type) {

    # 3.1) Preparing data
    if (dataset == "TCGA") {
        RNAseq_TPM <- RNAseq_TCGA_TPM_all
    } else if (dataset == "GTEx") {
        RNAseq_TPM <- RNAseq_GTEx_TPM_all
    }
    # Log2(TPM+1) as in DepMap data
    RNAseq_TPM <- log2(RNAseq_TPM+1)

    # Filter genes
    # NMD genes vs all
    if (genes_type == "all") {
        RNAseq_TPM_filt <- RNAseq_TPM
    } else if (genes_type == "NMD") {
        RNAseq_TPM_filt <- RNAseq_TPM[rownames(RNAseq_TPM) %in% NMD_ensembl_genes, ]
    }
    # Low Expressed
    RNAseq_TPM_filt2 <- RNAseq_TPM_filt[rowSums(log2(RNAseq_TPM_filt) >= 0.5) >= round(length(colnames(RNAseq_TPM_filt)) * 0.50), ]
    print(dim(RNAseq_TPM_filt2))

    # Low Variance
    RNAseq_TPM_all_t <- data.frame(t(RNAseq_TPM_filt2))
    RNAseq_TPM_all_t <- na.omit(RNAseq_TPM_all_t)
    low_var_cols <- nearZeroVar(RNAseq_TPM_all_t)
    if (length(low_var_cols) != 0) {
        RNAseq_TPM_all_t_filt <- RNAseq_TPM_all_t[, -low_var_cols]
    } else {
        RNAseq_TPM_all_t_filt <- RNAseq_TPM_all_t
    }
    print(dim(RNAseq_TPM_all_t_filt))

    # 3.2) Merge TCGA/GTEx with Cell Line data
    dat_merged <- merge(t(RNAseq_TPM_all_t_filt), RNAseq_DepMap_TPM_CL_t, by.x = "row.names", by.y = "row.names")
    # dat_merged <- merge(t(RNAseq_TPM_all_t_filt), CL_list[[1]], by.x = "Row.names", by.y = "gene_id")
    # dat_merged <- merge(dat_merged, CL_list[[3]], by.x = "Row.names", by.y = "gene_id")
    # dat_merged <- merge(dat_merged, CL_list[[4]], by.x = "Row.names", by.y = "gene_id")
    dat_merged <- data.frame(dat_merged)
    rownames(dat_merged) <- dat_merged$Row.names
    dat_merged$Row.names <- NULL
    RNAseq_TPM_cancer_and_CL <- data.frame(t(dat_merged))
    RNAseq_TPM_cancer_and_CL$type <- "CL"
    RNAseq_TPM_cancer_and_CL[grep(substr(dataset, 1, 3), rownames(RNAseq_TPM_cancer_and_CL)), "type"] <- dataset
    table(RNAseq_TPM_cancer_and_CL$type)
    # 3.3) Data alignment between tumors and cell lines
    # (From Marina)
    # For the alignment of TCGA and CL data, we first applied quantile normalization (R package preprocessCore 1.46.0) and next applied ComBat (R package sva 3.32.1),
    # a batch effect correction method. We used ComBat as if our dataset was the TCGA and CL data combined, and the batch effects labels were whether a sample
    # belongs to TCGA or CL (for MET) or a sample belongs to TCGA, GDSC, or CLLE (for GE). We applied this method for GE, MET, CNA, MS96, and RMD. For validation,
    # we calculated a PCA, subsampling TCGA data to match the number of CL samples (stratified by cancer types). In addition, we calculated Elastic Net classifiers
    # to predict (in the processed dataset) TCGA versus CL and calculated the AUC and AUPRC to check whether the process of alignment is being successful or not.
    # In addition to the chosen adjustment method, we tested other approaches based on canonical correlation analysis, partial least squares, and PCA, which did
    # not exceed accuracy of ComBat and therefore were not examined further.

    # apply the function to our matrix
    RNAseq_TPM_cancer_and_CL[1:5, 1:5]
    RNAseq_TPM_cancer_and_CL_norm <- qn(RNAseq_TPM_cancer_and_CL, dataset)
    RNAseq_TPM_cancer_and_CL_norm[1:5, 1:5]
    #RNAseq_TPM_cancer_and_CL_norm <- RNAseq_TPM_cancer_and_CL

    # Merge Gene expression and NMDeff (our phenotype to predict)
    if (dataset == "TCGA") {
        rownames(RNAseq_TPM_cancer_and_CL_norm)[grep("TCGA", rownames(RNAseq_TPM_cancer_and_CL_norm))] <- gsub("\\.", "-", substr(rownames(RNAseq_TPM_cancer_and_CL_norm)[grep("TCGA", rownames(RNAseq_TPM_cancer_and_CL_norm))], 1, 12))
        dat <- merge(RNAseq_TPM_cancer_and_CL_norm, sample_NMD_efficiencies_TCGA, by.x = "row.names", by.y = "sample", all.x = TRUE)
        rownames(dat) <- dat$Row.names
        dat$Row.names <- NULL
        if (target_var == "ASE_stopgain_0.2") {
            outlier_samples <- c("TCGA-JY-A6FG")
        } else if (target_var == "endogenous_NMD_global_2_shared") {
            outlier_samples <- c("TCGA-OR-A5J1", "TCGA-OR-A5J3", "TCGA-JY-A6FD", "TCGA-OR-A5J6", "TCGA-OR-A5J5", "TCGA-OR-A5J2", "TCGA-OR-A5J7", "TCGA-JY-A6FG")
        }
        dat <- dat[!rownames(dat) %in% outlier_samples, ]
        cols <- grep(paste0("ENSG|", target_var, "|type|sample$|Dim.*|cancer_subtype"), colnames(dat))
        dat <- dat[, cols]
    } else if (dataset == "GTEx") {
        rownames(RNAseq_TPM_cancer_and_CL_norm)[grep("GTE", rownames(RNAseq_TPM_cancer_and_CL_norm))] <- gsub("\\.", "-", rownames(RNAseq_TPM_cancer_and_CL_norm)[grep("GTE", rownames(RNAseq_TPM_cancer_and_CL_norm))])
        dat <- merge(RNAseq_TPM_cancer_and_CL_norm, sample_NMD_efficiencies_GTEx, by.x = "row.names", by.y = "sample_full_barcode", all.x = TRUE)
        rownames(dat) <- dat$Row.names
        dat$Row.names <- NULL
        cols <- grep(paste0("ENSG|", target_var, "|type|sample_full_barcode$|tissue"), colnames(dat))
        dat <- dat[, cols]
    }
    # Remove NAs (except in CLs)
    dat <- dat[(!(is.na(dat[, target_var])) & (dat$type != "CL") | (dat$type == "CL")), ]

    # Check and remove genes that correlate with covariates of interest

    # covariate <- "CNV_burden"
    # gene_cols <- grep("ENSG",colnames(dat))
    # cov_col <- grep("covariate",colnames(dat))
    # dat[,c(gene_cols,)]

    # 2.1.2) Data partitioning
    # Remove CL from the data
    dat_to_model <- dat[dat$type != "CL", ]
    df_results <- data.frame(
        dataset = NA, NMDeff_method = NA, model = NA, genes_type = NA, Rsquare_train = NA, Rsquare_test = NA,
        RMSE_train = NA, RMSE_test = NA, train_ss = NA, test_ss = NA, genes_size = NA
    )
    # Train vs test
    set.seed(100)
    index <- sample(1:nrow(dat_to_model), 0.7 * nrow(dat_to_model))
    train <- dat_to_model[index, ] # Create the training data
    test <- dat_to_model[-index, ] # Create the test data
    dim(train)
    dim(test)
    df_results$train_ss <- nrow(train)
    df_results$test_ss <- nrow(test)
    df_results$genes_size <- ncol(train)
    # Scale
    numeric_cols <- grep(paste0("ENSG|Dim.|sample_lib_size|type"), colnames(dat_to_model))
    # pre_proc_val <- preProcess(train[,numeric_cols], method = c("center", "scale")) --> Not needed, already QN
    # train[,numeric_cols] <- predict(pre_proc_val, train[,numeric_cols])
    # test[,numeric_cols] <- predict(pre_proc_val, test[,numeric_cols])
    # summary(train)

    # 2.1.3) Create dummy variables
    # The glmnet function does not work with dataframes, so we need to create a numeric matrix
    # for the training features and a vector of target values.
    options(expressions = 5e5)
    cols_reg <- grep(paste0("ENSG|", target_var, "$"), colnames(dat_to_model))
    # cols_reg <- grep(paste0("ENSG|cancer_type$|tissue|", target_var, "$"), colnames(dat_to_model))
    # cols_reg <- grep(paste0("ENSG|cancer_subtype|tissue|sample_lib_size|Dim\\..*|",target_var),colnames(dat_to_model))
    dummyVarsModel <- paste0("dummyVars(", target_var, " ~ ., data = dat_to_model[,cols_reg])")
    dummies <- eval(parse(text = dummyVarsModel))
    train_dummies <- predict(dummies, newdata = train[, cols_reg])
    test_dummies <- predict(dummies, newdata = test[, cols_reg])
    # cols_reg <- grep(paste0("ENSG|cancer_type.TCGA-SKCM|cancer_type.TCGA-OV|cancer_type.TCGA-STAD|", target_var, "$"), colnames(train_dummies))
    # train_dummies <- train_dummies[,cols_reg]
    # test_dummies <- test_dummies[,cols_reg]
    print(dim(train_dummies))
    print(dim(test_dummies))

    x <- as.matrix(train_dummies)
    y_train <- train[, target_var]
    x_test <- as.matrix(test_dummies)
    y_test <- test[, target_var]
    lambdas <- 10^seq(2, -3, by = -.1)

    if (model == "LASSO") {
        # 2.2) ##################### LASSO (coeff to 0) #######################
        # Setting alpha = 1 implements lasso regression
        reg_model_cv <- cv.glmnet(x, y_train, alpha = 1, lambda = lambdas, standardize = TRUE, nfolds = 5)
        lambda_best <- reg_model_cv$lambda.min
        # Final model
        final_reg_model <- glmnet(x, y_train, alpha = 1, lambda = lambda_best, standardize = TRUE)
    } else if (model == "ridge") {
        # 2.3) ####################### Ridge (some coeff to 0) #######################
        reg_model_cv <- cv.glmnet(x, y_train, alpha = 0, lambda = lambdas)
        lambda_best <- reg_model_cv$lambda.min
        # Final Model
        final_reg_model <- glmnet(x, y_train, alpha = 0, family = "gaussian", lambda = lambda_best)
    } else if (model == "elastic") {
        # 2.4) ####################### Elastic Net (between) #######################
        # Set training control
        train_cont <- trainControl(
            method = "repeatedcv",
            number = 10,
            repeats = 5,
            search = "random",
            verboseIter = TRUE
        )
        # Train the model
        train_model <- paste0("train(", target_var, " ~ .,
                                    data = train,
                                    method = 'glmnet',
                                    preProcess = c('center', 'scale'),
                                    tuneLength = 10,
                                    trControl = train_cont)")
        final_reg_model <- eval(parse(text = train_model))
        # Best tuning parameter
        final_reg_model$bestTune
        lambda_best <- NA
    }
    print(summary(final_reg_model))
    print(head(coef(final_reg_model)))
    # Prediction on train and test data
    if (is.na(lambda_best)) {
        predictions_train <- predict(final_reg_model, newx = x)
        predictions_test <- predict(final_reg_model, newx = x_test)
    } else {
        predictions_train <- predict(final_reg_model, s = lambda_best, newx = x)
        predictions_test <- predict(final_reg_model, s = lambda_best, newx = x_test)
    }
    # Evaluation on train and test data
    train_res <- eval_results(y_train, predictions_train, train)
    test_res <- eval_results(y_test, predictions_test, test)
    df_results$Rsquare_train <- train_res$Rsquare
    df_results$RMSE_train <- train_res$RMSE
    df_results$Rsquare_test <- test_res$Rsquare
    df_results$RMSE_test <- test_res$RMSE
    df_results
    # 4) Apply the model to each validation dataset (Cell Line)
    cols_reg <- grep(paste0("ENSG"), colnames(dat))
    x_CL <- dat[dat$type == "CL", cols_reg]
    # Swap NAs by 0
    x_CL[is.na(x_CL)] <- 0

    # Tissue
    # x_CL$CL <- gsub("_.*", "", rownames(x_CL))
    # x_CL[,"cancer_type.TCGA-OV"] <- 0
    # x_CL[,"cancer_type.TCGA-SKCM"] <- 0
    # x_CL[,"cancer_type.TCGA-STAD"] <- 0
    # x_CL[1:8,"cancer_type.TCGA-SKCM"] <- 1
    # x_CL[9:16,"cancer_type.TCGA-OV"] <- 1
    # x_CL[17:22,"cancer_type.TCGA-STAD"] <- 1
    CL_predictions <- data.frame(predict(final_reg_model, newx = as.matrix(x_CL)))

    # 5) Plot 
    rownames(CL_predictions) <- gsub("\\.","-",rownames(CL_predictions))
    CL_predictions <- merge(CL_predictions,DepMap_CNV_arm_level_type[,c("DepMap_ID","X1q")], by.x = "row.names", by.y = "DepMap_ID", all.x = TRUE)
    colnames(CL_predictions) <- c("CL","NMDeff","Chr_1q_CNV")   
    # CL_predictions$CL <- gsub("_.*", "", rownames(CL_predictions))
    # CL_predictions$condition <- unlist(lapply(strsplit(rownames(CL_predictions), "_"), function(x) {
    #     x[2]
    # }))
    # CL_predictions$condition <- ifelse(CL_predictions$condition == "WT","1q_trisomy","1q_disomy")
    CL_predictions <- na.omit(CL_predictions)
    CL_predictions$Chr_1q_CNV <- factor(CL_predictions$Chr_1q_CNV, levels = c(1,0,-1))
    if(target_var == "endogenous_NMD_global_2_shared") {
        NMD_method_char <- "Endogenous"
    } else {
        NMD_method_char <- "ASE"
    }
    CL_predictions$NMDeff <- scale(CL_predictions$NMDeff)
    combinations <- combn(names(table(CL_predictions$Chr_1q_CNV)), 2, simplify = FALSE)
    
    p <- ggplot(CL_predictions, aes(x = Chr_1q_CNV, y = -(NMDeff), fill = Chr_1q_CNV)) +
            geom_violin() + scale_fill_brewer(palette = "Dark2", direction = -1) + 
            geom_boxplot(width=0.3, color="black", alpha=0.2) +
            scale_x_discrete(labels = c("Amplification","Neutral","Deletion")) +
            xlab("CNA state - Chr 1q") + ylab("NMD efficiency") +
            #ggtitle(paste0("Cell-NMDeff based on ",dataset," - ",NMD_method_char)) +
            theme_bw(base_size = 25) +
            # geom_jitter(aes(fill = Chr_1q_CNV), 
            #             position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0, dodge.width = 0.9),
            #             alpha = 0.25, size = 1) +
            theme(axis.text.x = element_text(size = 35),
                axis.title.x = element_text(size = 40),
                axis.text.y = element_text(size = 35),
                axis.title.y = element_text(size = 40),
                plot.title = element_text(size = 45, hjust = 0.5),
                legend.position = "none") +
            stat_compare_means(comparisons = combinations, size = 12,
                              label.y = c(1,1.5,2),
                              label = "p.format", method = "wilcox.test", hide.ns = TRUE) +
        annotate("text",
                x = 1:length(table(CL_predictions$Chr_1q_CNV)),
                y = aggregate( -(NMDeff) ~ Chr_1q_CNV, CL_predictions, median)[ , 2],
                label = table(CL_predictions$Chr_1q_CNV),
                col = "black",
                vjust = - 1,
                size = 12)
    output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/CNV_associations/chr1q_validation/",dataset,"_",target_var,"_",model,"_genes_",genes_type,"_CL_NMDeff_predictions_Girish_chr1qLoss.png")
    png(output_path, width = 4500, height = 3500, res = 300)
    # p <- cowplot::plot_grid(plotlist=list(p1,p2), labels = "AUTO", align = "v", ncol = 2, nrow = 1)
    print(p)
    dev.off()
    return(list(df_results = df_results, CL_predictions = CL_predictions))
}

for (dataset in datasets) {
    for (model in models) {
        for (target_var in NMDeff_methods) {
            NMDeff_gene_model_validation(dataset = dataset, model = model, 
                                        target_var = target_var, genes_type = "all")
        }
    }
}

res_df <- NMDeff_gene_model_validation(dataset = "TCGA", model = "ridge", 
                            target_var = "endogenous_NMD_global_2_shared", genes_type = "all")
write.table(res_df[["CL_predictions"]], file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig3/Fig3E_2.txt", 
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
saveRDS(res_df[["CL_predictions"]], "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig3/Fig3E_2.RData")

# Supplementary
res_df <- NMDeff_gene_model_validation(dataset = "GTEx", model = "ridge", 
                            target_var = "endogenous_NMD_global_2_shared", genes_type = "all")
write.table(res_df[["CL_predictions"]], file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig12/SuppFig12B.txt", 
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
saveRDS(res_df[["CL_predictions"]], "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig12/SuppFig12B.RData")
