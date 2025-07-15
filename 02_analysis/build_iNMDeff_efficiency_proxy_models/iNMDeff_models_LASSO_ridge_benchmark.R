library("ggplot2")
library("caret")
library("glmnet")

# Compute R^2 from true and predicted values
eval_results <- function(true, predicted, df) {
  SSE <- sum((predicted - true)^2)
  SST <- sum((true - mean(true,na.rm =TRUE))^2)
  R_square <- 1 - SSE / SST
  RMSE = sqrt(SSE/nrow(df))
  # Model performance metrics
  data.frame(
    RMSE = RMSE,
    Rsquare = R_square
  )
}

# Start R with --> R --max-ppsize=5000000
# 1) Data

# 1.1) samples NMDeff pantissue
# TCGA
NMD_genesets <- c("endogenous_NMD_global_2_shared","ASE_stopgain_0.01","ASE_stopgain_0.2","PTCs_stopgain_NMD_triggering")
# 1.1) sample NMD efficiencies TCGA
# PTC // ASE // Endogenous
sample_NMD_efficiencies_TCGA_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/cancers/pancancer/NMD_efficiencies_TCGA.txt")
sample_NMD_efficiencies_TCGA <- read.table(file = sample_NMD_efficiencies_TCGA_path, header = TRUE, sep = "\t")
# Convert some columns to factors
factor_cols <- c("cancer_type","cancer_type_strat","cancer_subtype","LF_remove","purity_remove", "MSI_status",
                "batch_portion","batch_plate","batch_center","batch_vial","TCGA_full_barcode")
sample_NMD_efficiencies_TCGA[factor_cols] <- lapply(sample_NMD_efficiencies_TCGA[factor_cols], factor) 
cancers <- unique(sample_NMD_efficiencies_TCGA$cancer_type_strat)
# Filters for TCGA
sample_NMD_efficiencies_TCGA[which(sample_NMD_efficiencies_TCGA$ASE_num_PTCs_0.2 < 3),c("ASE_stopgain_0.2")] <- NA
sample_NMD_efficiencies_TCGA[which(sample_NMD_efficiencies_TCGA$ASE_num_PTCs_0.01 < 3),c("ASE_stopgain_0.01")] <- NA
#sample_NMD_efficiencies_TCGA_filt <- sample_NMD_efficiencies_TCGA[sample_NMD_efficiencies_TCGA$PTC_num_PTCs >= 8,]

# GTEx
sample_NMD_efficiencies_GTEx_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/GTEx/pantissue/NMD_efficiencies_GTEx.txt"
sample_NMD_efficiencies_GTEx <- read.table(file = sample_NMD_efficiencies_GTEx_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# Convert some columns to factors
factor_cols <- c("tissue","sample")
sample_NMD_efficiencies_GTEx[factor_cols] <- lapply(sample_NMD_efficiencies_GTEx[factor_cols], factor) 
GTEx_tissues <- unique(sample_NMD_efficiencies_GTEx$acronyms)
# Filters for GTEx
sample_NMD_efficiencies_GTEx[which(sample_NMD_efficiencies_GTEx$ASE_num_PTCs_0.2 < 3),c("ASE_stopgain_0.2")] <- NA
sample_NMD_efficiencies_GTEx[which(sample_NMD_efficiencies_GTEx$ASE_num_PTCs_0.01 < 3),c("ASE_stopgain_0.01")] <- NA

# 1.2) Pantissue TPM gene expression
# TCGA
# TCGA_tissues <- as.character(read.table(file = "/g/strcombio/fsupek_home/gpalou/data/TCGA_RNAseq_quantification/TCGA_projects_names.txt",sep = "\t")$V1)
# print("1) RNAseq TPM")
# RNAseq_sample_names <- list()
# RNAseq_TCGA_TPM_all <- c()
# for (i in seq(1:length(TCGA_tissues))) {
#   TCGA_cancer <- as.character(TCGA_tissues[i])
#   print(paste0(i," --> ",TCGA_cancer))
#   RNAseq_error <- FALSE
#   tryCatch( {   RNAseq_TCGA_TPM <- read.table(file = gsub("\\[X\\]",TCGA_cancer, "/g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/primary_tumor/[X]/[X]_RNAseq_matrix_TPM_gene.txt"),
#                                               header = TRUE, sep = "\t", row.names = 1)
#   colnames(RNAseq_TCGA_TPM) <- gsub("\\.","-",substr(colnames(RNAseq_TCGA_TPM),1,12))
#   # keep track of sample names for that cancer
#   RNAseq_sample_names[[TCGA_cancer]] <- colnames(RNAseq_TCGA_TPM)
#   # Join matrices in one
#   if (length(RNAseq_TCGA_TPM_all)==0) {RNAseq_TCGA_TPM_all <- RNAseq_TCGA_TPM}
#   else{RNAseq_TCGA_TPM_all <- cbind(RNAseq_TCGA_TPM_all,RNAseq_TCGA_TPM)}
#   # Check rownames are the same
#   print(paste0("ENSEMBL transcripts IDs are the same? ",table(rownames(RNAseq_TCGA_TPM)==rownames(RNAseq_TCGA_TPM_all))))
#   }
#   ,error = function(e) {
#     RNAseq_error <<- TRUE
#     print(e)
#     })
#   if (isTRUE(RNAseq_error)) {next}
# }

output_path <- "/g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/primary_tumor/TCGA_RNAseq_matrix_TPM_gene.txt"
# write.table(RNAseq_TCGA_TPM_all, file = output_path, sep = "\t", quote = FALSE,
#             col.names = TRUE, row.names = TRUE)
RNAseq_TCGA_TPM_all <- read.table(file = output_path, header = TRUE, sep = "\t", row.names = 1)
print("Dimensions -->")
print(dim(RNAseq_TCGA_TPM_all))
# Check Samples with NAs
na_counts <- colSums(is.na(RNAseq_TCGA_TPM_all))
cols_na <- names(na_counts[na_counts > 0])
RNAseq_TCGA_TPM_all <- RNAseq_TCGA_TPM_all[,!colnames(RNAseq_TCGA_TPM_all) %in% cols_na]
# Check Genes with NAs
na_counts <- rowSums(is.na(RNAseq_TCGA_TPM_filt))
rows_na <- names(na_counts[na_counts > 0])
print("Dimensions -->")
print(dim(RNAseq_TCGA_TPM_all))

# GTEx
GTEx_path <- "/g/strcombio/fsupek_cancer1/gpalou/GTEx/RNAseq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct"
RNAseq_GTEx_TPM_all <- read.table(file = GTEx_path, header = TRUE, sep = "\t", row.names = 1, skip = 2)
print("Dimensions -->")
print(dim(RNAseq_GTEx_TPM_all))
RNAseq_GTEx_TPM_all <- RNAseq_GTEx_TPM_all[-grep("PAR",rownames(RNAseq_GTEx_TPM_all)),]
rownames(RNAseq_GTEx_TPM_all) <- gsub("(.*)\\..*","\\1",rownames(RNAseq_GTEx_TPM_all))
RNAseq_GTEx_TPM_all$Description <- NULL

# Check Samples with NAs
# na_counts <- colSums(is.na(RNAseq_GTEx_TPM_all))
# cols_na <- names(na_counts[na_counts > 0])
# RNAseq_GTEx_TPM_filt <- RNAseq_GTEx_TPM_all[,!colnames(RNAseq_GTEx_TPM_all) %in% cols_na]
# # Check Genes with NAs
# na_counts <- rowSums(is.na(RNAseq_GTEx_TPM_all))
# rows_na <- names(na_counts[na_counts > 0])

# 1.3) NMD target genes
NMD_genes <- read.table(file = "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/NMD_targets/ENSEMBL_transcripts/NMD_global_2_shared_ensembl.txt",
                        header = TRUE, sep = "\t")
NMD_ensembl_genes <- as.character(unique(NMD_genes$ensembl_gene_id))

# 1.4) CNV-PCs pancancer
scale <- TRUE
center <- TRUE
alpha <- "3e-04"
num_PCs <- "100"
tryCatch({
    TCGA_CNV_PCA_ind <- read.table(file = paste0("/g/strcombio/fsupek_cancer1/gpalou/TCGA_CNV_PCA/pancancer_sparse_PCA_ind_",alpha,"_robust_no_num_PCs_",num_PCs,".txt"),
                                        header = TRUE, stringsAsFactors = FALSE, sep = "\t")
    rownames(TCGA_CNV_PCA_ind) <- gsub("\\.","-",rownames(TCGA_CNV_PCA_ind))
    # Remove PCs with 0
    cols <- colnames(TCGA_CNV_PCA_ind)[which( colSums(TCGA_CNV_PCA_ind) != 0 )]
    TCGA_CNV_PCA_ind <- TCGA_CNV_PCA_ind[,cols]
    print("PCA dimensions --> ")
    print(dim(TCGA_CNV_PCA_ind))
    nPCs <- ncol(TCGA_CNV_PCA_ind)
    # Scale PCs
    TCGA_CNV_PCA_ind <- data.frame(scale(TCGA_CNV_PCA_ind, scale = scale, center = center))
}, error = function(e){
    print(e)
    }
)
# Merge
sample_NMD_efficiencies_TCGA <- merge(sample_NMD_efficiencies_TCGA,TCGA_CNV_PCA_ind, by.x = "sample", by.y = "row.names", all.x = TRUE)

# 2) NMDeff gene level model

NMDeff_gene_level_model <- function(dataset, target_var, model, genes_type = "all") {

  benchmark_df <- data.frame(dataset = NA, NMDeff_method = NA, model = NA, genes_type = NA, Rsquare_train = NA, Rsquare_test = NA,
                            RMSE_train = NA, RMSE_test = NA, train_ss = NA, test_ss = NA, genes_size = NA)
  benchmark_df[,1:4] <- c(dataset,target_var,model,genes_type)
  
  # 2.1) Preparing data
  if (dataset == "TCGA") { 
    RNAseq_TPM <- RNAseq_TCGA_TPM_all
  } else if (dataset == "GTEx") {
    RNAseq_TPM <- RNAseq_GTEx_TPM_all
  }
  # 2.1.1) Filter genes
  # NMD genes?
  if (genes_type == "all") {
    RNAseq_TPM_filt <- RNAseq_TPM
  } else if (genes_type == "NMD") {
    RNAseq_TPM_filt <- RNAseq_TPM[rownames(RNAseq_TPM) %in% NMD_ensembl_genes,]
  }
  # Low Expressed
  RNAseq_TPM_filt2 <- RNAseq_TPM_filt[rowSums(log2(RNAseq_TPM_filt) >= 0.5) >= round(length(colnames(RNAseq_TPM_filt)) * 0.50),]
  print(dim(RNAseq_TPM_filt2))

  #################
  # RNAseq_TPM_filt2 <- RNAseq_TPM_filt2[sample(1:nrow(RNAseq_TPM_filt2),500),]
  #################
  
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

  # 2.1.1) Merge Gene expression and NMDeff

  if (dataset == "TCGA") {
    rownames(RNAseq_TPM_all_t_filt) <- gsub("\\.","-",substr(rownames(RNAseq_TPM_all_t_filt),1,12))
    dat <- merge(sample_NMD_efficiencies_TCGA,RNAseq_TPM_all_t_filt, by.x = "sample", by.y = "row.names", all.x = TRUE)
    if (target_var == "ASE_stopgain_0.2") {
      outlier_samples <- c("TCGA-JY-A6FG")
    } else if (target_var == "endogenous_NMD_global_2_shared") {
      outlier_samples <- c("TCGA-OR-A5J1","TCGA-OR-A5J3","TCGA-JY-A6FD","TCGA-OR-A5J6","TCGA-OR-A5J5","TCGA-OR-A5J2","TCGA-OR-A5J7","TCGA-JY-A6FG")
    }
    dat <- dat[!dat$sample %in% outlier_samples,]
    rownames(dat) <- dat$sample
    cols = grep(paste0("ENSG|",target_var,"|sample$|Dim.*|cancer_subtype"),colnames(dat))
    dat <- dat[,cols]
    dat$sample <- NULL
  } else if (dataset == "GTEx") {
    rownames(RNAseq_TPM_all_t_filt) <- gsub("\\.","-",rownames(RNAseq_TPM_all_t_filt))
    dat <- merge(sample_NMD_efficiencies_GTEx,RNAseq_TPM_all_t_filt, by.x = "sample_full_barcode", by.y = "row.names", all.x = TRUE)
    rownames(dat) <- dat$sample_full_barcode
    cols <- grep(paste0("ENSG|",target_var,"|sample_full_barcode$|tissue"),colnames(dat))
    dat <- dat[,cols]
    dat$sample_full_barcode <- NULL
  }

  #################
  # dat <- dat[sample(1:nrow(dat),500),]
  #################

  # Remove NAs
  dat <- dat[!is.na(dat[,target_var]),]
  # 2.1.2) Data partitioning and scale numeric variables (gene expression)
  # Train vs test
  set.seed(100) 
  index <- sample(1:nrow(dat), 0.7*nrow(dat)) 
  train <- dat[index,] # Create the training data 
  test <- dat[-index,] # Create the test data
  dim(train)
  dim(test)
  benchmark_df$train_ss <- nrow(train)
  benchmark_df$test_ss <- nrow(test)
  benchmark_df$genes_size <- ncol(train)
  # Scale
  numeric_cols <- grep(paste0("ENSG|Dim.|sample_lib_size"),colnames(dat))
  pre_proc_val <- preProcess(train[,numeric_cols], method = c("center", "scale"))
  train[,numeric_cols] <- predict(pre_proc_val, train[,numeric_cols])
  test[,numeric_cols] <- predict(pre_proc_val, test[,numeric_cols])
  #summary(train)

  # 2.1.3) Create dummy variables
  # The glmnet function does not work with dataframes, so we need to create a numeric matrix 
  # for the training features and a vector of target values.
  options(expressions = 5e5)
  cols_reg <- grep(paste0("ENSG|cancer_subtype|tissue|sample_lib_size|Dim\\..*|",target_var),colnames(dat))
  dummyVarsModel <- paste0("dummyVars(",target_var," ~ ., data = dat[,cols_reg])")
  dummies <- eval(parse(text=dummyVarsModel))
  train_dummies <- predict(dummies, newdata = train[,cols_reg])
  test_dummies <- predict(dummies, newdata = test[,cols_reg])
  print(dim(train_dummies))
  print(dim(test_dummies))

  x = as.matrix(train_dummies)
  y_train = train[,target_var]
  x_test = as.matrix(test_dummies)
  y_test = test[,target_var]
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
    final_reg_model <- glmnet(x, y_train, alpha = 0, family = 'gaussian', lambda = lambda_best)
  } else if (model == "elastic") {
    # 2.4) ####################### Elastic Net (between) #######################
    # Set training control
    train_cont <- trainControl(method = "repeatedcv",
                                  number = 10,
                                  repeats = 5,
                                  search = "random",
                                  verboseIter = TRUE)
    # Train the model
    train_model <- paste0("train(",target_var," ~ .,
                              data = train,
                              method = 'glmnet',
                              preProcess = c('center', 'scale'),
                              tuneLength = 10,
                              trControl = train_cont)")
    final_reg_model <- eval(parse(text=train_model))
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
  benchmark_df$Rsquare_train <- train_res$Rsquare
  benchmark_df$RMSE_train <- train_res$RMSE
  benchmark_df$Rsquare_test <- test_res$Rsquare
  benchmark_df$RMSE_test <- test_res$RMSE
  return(benchmark_df)
}

# y_predicted <- predict(best_model, s = best_lambda, newx = x)

# #find SST and SSE
# sst <- sum((y - mean(y))^2)
# sse <- sum((y_predicted - y)^2)

# #find R-Squared
# rsq <- 1 - sse/sst
# rsq

datasets <- c("TCGA","GTEx")
#models <- c("LASSO","ridge","elastic")
models <- c("LASSO","ridge")
NMDeff_methods <- c("ASE_stopgain_0.2","endogenous_NMD_global_2_shared")

benchmark_df_final <- data.frame(dataset = NA, NMDeff_method = NA, model = NA, genes_type = NA, Rsquare_train = NA, Rsquare_test = NA,
                            RMSE_train = NA, RMSE_test = NA, train_ss = NA, test_ss = NA, genes_size = NA)
for (dataset in datasets) {
  for (NMDeff_method in NMDeff_methods) {
    for (model in models) {
      benchmark_df <- NMDeff_gene_level_model(dataset = dataset, target_var = NMDeff_method, model = model, genes_type = "all")
      print(benchmark_df)
      benchmark_df_final <- rbind(benchmark_df_final,benchmark_df)
    }
  }
}

# Modify dataset to plot
benchmark_df_final <- benchmark_df_final[-1,]
yvars <- c("Rsquare_train","Rsquare_test")
benchmark_df_final_stacked <- stack(benchmark_df_final[,yvars])
benchmark_df_final_stacked$ind <- gsub("Rsquare_","",benchmark_df_final_stacked$ind)
colnames(benchmark_df_final_stacked)[2] <- "type"
final_df <- cbind(benchmark_df_final_stacked,
            rbind(benchmark_df_final[,!colnames(benchmark_df_final) %in% yvars],benchmark_df_final[,!colnames(benchmark_df_final) %in% yvars]))
output_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/NMDeff_gene_model/NMDeff_gene_model_results.txt"
# write.table(final_df, file = output_path, sep = "\t", quote = FALSE,
#             col.names = TRUE, row.names = TRUE)
final_df <- read.table(file = output_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# Save
write.table(final_df, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig12/SuppFig12A.txt", 
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
saveRDS(final_df, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig12/SuppFig12A.RData")

# Remove ASE
library(dplyr)
library(RColorBrewer)

final_df_filt <- final_df %>%
                filter(NMDeff_method == "endogenous_NMD_global_2_shared")
# 3) Plot
png(paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/NMDeff_gene_model/NMDeff_gene_model_barplot.png"), width = 4500, height = 3500, res = 300)
p <- ggplot(data = final_df_filt, aes(x = model, y = values, fill = type)) +
  xlab("Model") + ylab("Rsquare") + theme_classic() +
  facet_wrap(~ dataset) + scale_fill_brewer(palette = "Paired", direction = -1) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  theme(plot.title = element_text(hjust = 0.5, size = 40),
        axis.title.x = element_text(color="black", size=35),
        axis.text.x = element_text(color="black", size=30),
        axis.title.y = element_text(color="black", size=35),
        axis.text.y = element_text(color="black", size=30),
        strip.text.x = element_text(size = 35),
        legend.title=element_blank(),
        legend.text=element_text(size=30),
        legend.position = "top",
        legend.key.size = unit(1, "cm"))
print(p)
dev.off()

