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

NMDeff_gene_model_validation <- function(dataset, model, target_var, dataset_to_predict) {

    # 1) Preparing dataset data
    if (dataset == "TCGA") {
        RNAseq_TPM <- RNAseq_TCGA_TPM_all
    } else if (dataset == "GTEx") {
        RNAseq_TPM <- RNAseq_GTEx_TPM_all
    }
    # Normalize by Log2(TPM+1) as in the datasets we want to predict
    if (dataset_to_predict == "Motzer_RCC") {
      RNAseq_TPM[] <- lapply(RNAseq_TPM, function(x) ifelse(x < 0.01, 0.01, x))
      RNAseq_TPM <- log2(RNAseq_TPM)
    } else {
      RNAseq_TPM <- log2(RNAseq_TPM+1)
    }
    # Filter genes
    RNAseq_TPM_filt <- RNAseq_TPM
    # Low Expressed
    RNAseq_TPM_filt2 <- RNAseq_TPM_filt[rowSums(RNAseq_TPM_filt >= 0.5) >= round(length(colnames(RNAseq_TPM_filt)) * 0.50), ]
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

    # 2) Merge TCGA/GTEx with the other dataset
    if (dataset_to_predict == "Liu_SKCM") {
      Immuno_dataset <- Liu_SKCM
      Immuno_dataset$gene_id <- NULL
    } else if (dataset_to_predict == "Braun_RCC") {
      Immuno_dataset <- Braun_RCC
      Immuno_dataset$gene_id <- NULL
    } else if (dataset_to_predict == "Carrol_EAC") {
      Immuno_dataset <- Carrol_EAC
    } else if (dataset_to_predict == "Motzer_RCC") {
      Immuno_dataset <- Motzer_RCC
      Immuno_dataset$gene_id <- NULL
    } else if (dataset_to_predict == "Riaz_SKCM") {
      Immuno_dataset <- Riaz_SKCM
    } else if (dataset_to_predict == "Hartwig") {
      Immuno_dataset <- RNAseq_Hartwig_TPM_all
    }
    
    dat_merged <- merge(t(RNAseq_TPM_all_t_filt), Immuno_dataset, by.x = "row.names", by.y = "row.names")
    dat_merged <- data.frame(dat_merged)
    rownames(dat_merged) <- dat_merged$Row.names
    dat_merged$Row.names <- NULL
    RNAseq_TPM_cancer_and_immuno_dat <- data.frame(t(dat_merged))
    RNAseq_TPM_cancer_and_immuno_dat$type <- "Immuno_dat"
    RNAseq_TPM_cancer_and_immuno_dat[grep(substr(dataset, 1, 3), rownames(RNAseq_TPM_cancer_and_immuno_dat)), "type"] <- dataset
    table(RNAseq_TPM_cancer_and_immuno_dat$type)

    # Save tmp matrix
    out_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/phenotypes/immunotherapy/RNAseq_",dataset,"_and_",dataset_to_predict,".txt")
    write.table(RNAseq_TPM_cancer_and_immuno_dat, file = out_path, sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
    # RNAseq_TPM_cancer_and_immuno_dat <- read.table(file = out_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

    # 3) Data alignment between tumors and the other dataset
    # apply the function to our matrix
    rows <- which(RNAseq_TPM_cancer_and_immuno_dat$type == "Immuno_dat")[1:5]
    RNAseq_TPM_cancer_and_immuno_dat[rows-2, 1:5]
    RNAseq_TPM_cancer_and_immuno_dat <- qn(RNAseq_TPM_cancer_and_immuno_dat, dataset)
    RNAseq_TPM_cancer_and_immuno_dat[rows-2, 1:5]

    # Merge Gene expression and NMDeff (our phenotype to predict)
    if (dataset == "TCGA") {
        rownames(RNAseq_TPM_cancer_and_immuno_dat)[grep("TCGA", rownames(RNAseq_TPM_cancer_and_immuno_dat))] <- gsub("\\.", "-", substr(rownames(RNAseq_TPM_cancer_and_immuno_dat)[grep("TCGA", rownames(RNAseq_TPM_cancer_and_immuno_dat))], 1, 12))
        dat <- merge(RNAseq_TPM_cancer_and_immuno_dat, sample_NMD_efficiencies_TCGA, by.x = "row.names", by.y = "sample", all.x = TRUE)
        rownames(dat) <- dat$Row.names
        dat$Row.names <- NULL
        if (target_var == "ASE_PTC_NMD_triggering_0.2") {
            outlier_samples <- c("TCGA-JY-A6FG")
        } else if (target_var == "endogenous_NMD_Consensus") {
            outlier_samples <- c("TCGA-OR-A5J1", "TCGA-OR-A5J3", "TCGA-JY-A6FD", "TCGA-OR-A5J6", "TCGA-OR-A5J5", "TCGA-OR-A5J2", "TCGA-OR-A5J7", "TCGA-JY-A6FG")
        }
        dat <- dat[!rownames(dat) %in% outlier_samples, ]
        cols <- grep(paste0("ENSG|", target_var, "$|^type$|sample$|Dim.*"), colnames(dat))
        # cols <- grep(paste0("ENSG|", target_var, "$|^type$|sample$|Dim.*|cancer_subtype"), colnames(dat))
        dat <- dat[, cols]
    } else if (dataset == "GTEx") {
        rownames(RNAseq_TPM_cancer_and_immuno_dat)[grep("GTE", rownames(RNAseq_TPM_cancer_and_immuno_dat))] <- gsub("\\.", "-", rownames(RNAseq_TPM_cancer_and_immuno_dat)[grep("GTE", rownames(RNAseq_TPM_cancer_and_immuno_dat))])
        dat <- merge(RNAseq_TPM_cancer_and_immuno_dat, sample_NMD_efficiencies_GTEx, by.x = "row.names", by.y = "sample_full_barcode", all.x = TRUE)
        rownames(dat) <- dat$Row.names
        dat$Row.names <- NULL
        cols <- grep(paste0("ENSG|", target_var, "$|^type$|sample_full_barcode$"), colnames(dat))
        # cols <- grep(paste0("ENSG|", target_var, "$|^type$|sample_full_barcode$|tissue"), colnames(dat))
        dat <- dat[, cols]
    }
    # Remove NAs (except in Immuno_dat)
    dat <- dat[(!(is.na(dat[, target_var])) & (dat$type != "Immuno_dat") | (dat$type == "Immuno_dat")), ]

    ######
    # Check and remove genes that correlate with covariates of interest
    ######

    # 4) Data partitioning
    # Remove Immuno_dat from the data
    dat_to_model <- dat[dat$type != "Immuno_dat", ]
    # Columns are factors (for some reason), so, transform into numeric
    # dat_to_model <- dat_to_model %>%
    #   mutate(across(where(~ starts_with("ENSG") && is.factor(.) && !any(is.na(as.numeric(as.character(.))))), 
    #                 ~ as.numeric(as.character(.))))

    df_results <- data.frame(
        dataset = NA, NMDeff_method = NA, model = NA, Rsquare_train = NA, Rsquare_test = NA,
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
    if (dataset_to_predict == "Carrol_EAC") {
      set.seed(100)
      cols_keep <- sample(grep(paste0("ENSG"), colnames(dat_to_model)))[1:16000]
      cols_reg <- c(cols_keep,grep(paste0(target_var), colnames(dat_to_model)))
    } else {
      cols_reg <- grep(paste0("ENSG|", target_var, "$"), colnames(dat_to_model))
    }
    dummyVarsModel <- paste0("dummyVars(", target_var, " ~ ., data = dat_to_model[,cols_reg])")
    dummies <- eval(parse(text = dummyVarsModel))
    train_dummies <- predict(dummies, newdata = train[,cols_reg])
    test_dummies <- predict(dummies, newdata = test[, cols_reg])
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
    if (dataset_to_predict == "Carrol_EAC") {
      set.seed(100)
      cols_reg <- sample(grep(paste0("ENSG"), colnames(dat_to_model)))[1:16000]
    } else {
      cols_reg <- grep(paste0("ENSG"), colnames(dat))
    }
    
    x_Immuno_dat <- dat[dat$type == "Immuno_dat", cols_reg]
    x_Immuno_dat <- x_Immuno_dat[-1,]
    # Swap NAs by 0
    x_Immuno_dat[is.na(x_Immuno_dat)] <- 0
    # Convert columns to numeric
    # x_Immuno_dat <- x_Immuno_dat %>%
    #   mutate(across(where(~ starts_with("ENSG") && is.factor(.) && !any(is.na(as.numeric(as.character(.))))), 
    #                 ~ as.numeric(as.character(.))))

    # Merge
    Immuno_dat_predictions <- data.frame(predict(final_reg_model, newx = as.matrix(x_Immuno_dat)))   
    if (dataset_to_predict == "Liu_SKCM") {
      Immuno_dat_predictions <- merge(Immuno_dat_predictions, Liu_SKCM_metadata, by.x = "row.names", by.y = "sampleID", all.x = TRUE)
      colnames(Immuno_dat_predictions)[colnames(Immuno_dat_predictions) %in% c("s0","BR")] <- c("NMDeff","treatment_response")
      # Immuno_dat_predictions <- Immuno_dat_predictions[c("s0","BR","TimeToBR")]
      # Immuno_dat_predictions <- Immuno_dat_predictions[Immuno_dat_predictions$BR != "MR",]
      combinations <- combn(names(table(Immuno_dat_predictions$treatment_response)), 2, simplify = FALSE)
    } else if (dataset_to_predict == "Braun_RCC") {
      rownames(Immuno_dat_predictions) <- gsub("\\.","-",rownames(Immuno_dat_predictions))
      Immuno_dat_predictions <- merge(Immuno_dat_predictions, Braun_RCC_metadata, by.x = "row.names", by.y = "RNA_ID", all.x = TRUE)
      colnames(Immuno_dat_predictions)[colnames(Immuno_dat_predictions) %in% c("s0","ORR")] <- c("NMDeff","treatment_response")
      combinations <- combn(names(table(Immuno_dat_predictions$treatment_response)), 2, simplify = FALSE)
    } else if (dataset_to_predict == "Carrol_EAC") {
      rownames(Immuno_dat_predictions) <- gsub("\\.","-",rownames(Immuno_dat_predictions))
      Immuno_dat_predictions$sample_type <- substr(rownames(Immuno_dat_predictions),10,50)
      Immuno_dat_predictions$sampleID <- substr(rownames(Immuno_dat_predictions),1,8)
      Immuno_dat_predictions <- merge(Immuno_dat_predictions, Carrol_EAC_metadata, by.x = "sampleID", by.y = "sampleID", all.x = TRUE)
      colnames(Immuno_dat_predictions)[colnames(Immuno_dat_predictions) %in% c("s0","Clinical_benefit")] <- c("NMDeff","treatment_response")
      combinations <- combn(names(table(Immuno_dat_predictions$treatment_response)), 2, simplify = FALSE)
    } else if (dataset_to_predict == "Motzer_RCC") {
      Immuno_dat_predictions <- merge(Immuno_dat_predictions, Motzer_RCC_metadata, by.x = "row.names", by.y = "ID", all.x = TRUE)
      colnames(Immuno_dat_predictions)[colnames(Immuno_dat_predictions) %in% c("s0","PDL1FL")] <- c("NMDeff","treatment_response")
      combinations <- combn(names(table(Immuno_dat_predictions$treatment_response)), 2, simplify = FALSE)
    } else if (dataset_to_predict == "Riaz_SKCM") {
      Immuno_dat_predictions$patientID <- gsub("(Pt[0-9]{1,3})\\_.*","\\1",rownames(Immuno_dat_predictions))
      Immuno_dat_predictions$patient_full_ID <- rownames(Immuno_dat_predictions)
      Immuno_dat_predictions <- merge(Immuno_dat_predictions, Riaz_SKCM_metadata, by.x = "patientID", by.y = "Patient", all.x = TRUE)
      colnames(Immuno_dat_predictions)[colnames(Immuno_dat_predictions) %in% c("s0","Response")] <- c("NMDeff","treatment_response")
      combinations <- combn(names(table(Immuno_dat_predictions$treatment_response)), 2, simplify = FALSE)
    } else if (dataset_to_predict == "Hartwig") {
      Immuno_dat_predictions <- merge(Immuno_dat_predictions, Hartwig_metadata, by.x = "row.names", by.y = "sampleId", all.x = TRUE)
      colnames(Immuno_dat_predictions)[colnames(Immuno_dat_predictions) %in% c("s0")] <- c("NMDeff")
      Immuno_dat_predictions$treatment_response <- Immuno_dat_predictions$firstResponse
      Immuno_dat_predictions$treatment_response <- ifelse(Immuno_dat_predictions$treatment_response == "null",NA,Immuno_dat_predictions$treatment_response)
      combinations <- combn(names(table(Immuno_dat_predictions$treatment_response)), 2, simplify = FALSE)
    }
    
    # 5) Plot 
    p <- Immuno_dat_predictions %>% #filter(ImmunoPhenotype == "No_IF") %>%
        ggplot(aes(x = factor(treatment_response), y = NMDeff, fill = factor(treatment_response))) +
        geom_violin() + scale_fill_brewer(palette = "Dark2") + xlab("") +
        geom_boxplot(aes(fill = factor(treatment_response)), width=0.2, color="black", alpha=0.2) +
        geom_jitter(aes(fill = factor(treatment_response)), 
                    position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0, dodge.width = 0.9),
                    alpha = 0.25, size = 1) +        theme_bw(base_size = 25) +
        theme(axis.text.x = element_text(size = 18, angle = 45, hjust = 1)) +
      annotate("text",
              x = 1:length(table(Immuno_dat_predictions$treatment_response)),
              y = aggregate( NMDeff ~ treatment_response, Immuno_dat_predictions, median)[ , 2],
              label = table(Immuno_dat_predictions$treatment_response),
              col = "black",
              vjust = - 1,
              size = 5) #+
      # stat_compare_means(comparisons = combinations, size = 12,
      #                   # label.y = c(2.5,2.25,2.5,2.75,3),
      #                   label = "p.format", method = "wilcox.test", hide.ns = TRUE)# +
    output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/phenotypes/immunotherapy/validation_datasets/",dataset_to_predict,"/",dataset,"_",target_var,"_",model,"_Immuno_dat_NMDeff_predictions.png")
    png(output_path, width = 4500, height = 3000, res = 300)
    # p <- cowplot::plot_grid(plotlist=list(p1,p2), labels = "AUTO", align = "v", ncol = 2, nrow = 1)
    print(p)
    dev.off()
    return(Immuno_dat_predictions)
}

survival_curve <- function(Immuno_NMDeff, dataset_to_predict, percentile = NULL, KM_output = "yes", 
                survival_var, iNMDeff_model_type, binary_probability = NULL) {

  list_plots <- list()
  list_KM_objects <- list()
    for (NMD_method in unique(Immuno_NMDeff$NMD_method)) {
      print(NMD_method)
      if (dataset_to_predict == "Motzer_RCC") {
        Immuno_NMDeff_df <- Immuno_NMDeff[Immuno_NMDeff$NMD_method == NMD_method & Immuno_NMDeff$TRT01P != "Sunitinib",]
        # Immuno_NMDeff_df <- Immuno_NMDeff[Immuno_NMDeff$NMD_method == NMD_method,]
        if (survival_var == "OS") {
          survival_var_char <- ""
          status_variable <- ""
        } else if (survival_var == "PFS") {
          survival_var_char <- "PFS_P"
          status_variable <- "PFS_P_CNSR"
        } 
      } else if (dataset_to_predict == "Liu_SKCM") {
        Immuno_NMDeff_df <- Immuno_NMDeff[Immuno_NMDeff$NMD_method == NMD_method,]
        if (survival_var == "OS") {
          survival_var_char <- "OS"
          status_variable <- "dead"
        } else if (survival_var == "PFS") {
          survival_var_char <- "PFS"
          status_variable <- "progressed"
        }
      } else if (dataset_to_predict == "Carrol_EAC") {
        Immuno_NMDeff_df <- Immuno_NMDeff[Immuno_NMDeff$NMD_method == NMD_method,]
        if (survival_var == "OS") {
          survival_var_char <- "OS"
          status_variable <- "Status"
        } else if (survival_var == "PFS") {
          survival_var_char <- "PFS"
          status_variable <- "Progressed"
        }
      } else if (dataset_to_predict == "TCGA") {
        Immuno_NMDeff_df <- Immuno_NMDeff[Immuno_NMDeff$NMD_method == NMD_method,]
        if (survival_var == "OS") {
          survival_var_char <- "days_to_last_follow_up"
          status_variable <- "status"
        } else if (survival_var == "PFS") {
          survival_var_char <- "PFS"
          status_variable <- "PFS_status"
        }
      } else if (dataset_to_predict == "Riaz_SKCM") {
        Immuno_NMDeff_df <- Immuno_NMDeff[Immuno_NMDeff$NMD_method == NMD_method,]
        if (survival_var == "OS") {
          survival_var_char <- "time_to_death_weeks"
          status_variable <- "status"
        } else if (survival_var == "PFS") {
          survival_var_char <- ""
          status_variable <- ""
        }
      }
      Immuno_NMDeff_df <- data.frame(Immuno_NMDeff_df)

      if (iNMDeff_model_type == "continous") {
        NMD_classification <- "NMDeff_mean"
        Immuno_NMDeff_df$NMD_classification <- NA
        percentiles <- quantile(Immuno_NMDeff_df[,NMD_classification],probs = seq(0, 1, 0.01), na.rm = TRUE)
        thres_low <- as.numeric(percentiles[paste0(percentile,"%")])
        thres_high <- as.numeric(percentiles[paste0(100-percentile,"%")])
        Immuno_NMDeff_df[which(Immuno_NMDeff_df[,NMD_classification] <= thres_low),"NMD_classification"] <- "Low"
        Immuno_NMDeff_df[which(Immuno_NMDeff_df[,NMD_classification] >= thres_high),"NMD_classification"] <- "High"
        # thres_high <- as.numeric(percentiles[paste0(100-percentile,"%")])
        # Immuno_NMDeff_df[which(Immuno_NMDeff_df[,NMD_classification] < thres_high),"NMD_classification"] <- "Low"
        # Immuno_NMDeff_df[which(Immuno_NMDeff_df[,NMD_classification] >= thres_high),"NMD_classification"] <- "High"
       } else if (iNMDeff_model_type == "binary") {
          Immuno_NMDeff_df$NMD_classification <- case_when(
            Immuno_NMDeff_df$NMD_High_prob >= binary_probability ~ "NMD-High",
            Immuno_NMDeff_df$NMD_High_prob <= (1-binary_probability) ~ "NMD-Low",
            TRUE ~ "NMD-Mid"
          ) 
          Immuno_NMDeff_df <- Immuno_NMDeff_df %>%
            filter(NMD_classification != "NMD-Mid") %>%
            mutate(NMD_classification = factor(NMD_classification))
        table(Immuno_NMDeff_df$NMD_classification)
       }
      classification_table <- table(Immuno_NMDeff_df$NMD_classification)
      if (length(classification_table) > 1) {
        surv_obj <<- NULL
        surv_obj <<- with(Immuno_NMDeff_df, Surv( as.numeric(eval(parse(text = survival_var_char))), as.numeric(eval(parse(text = status_variable))) ))
        km_fit <<- survfit(surv_obj ~ NMD_classification, data = Immuno_NMDeff_df)
        print(km_fit)
        list_KM_objects[[NMD_method]] <- list(km_fit = km_fit, df = Immuno_NMDeff_df, surv_obj = surv_obj)

        p <- ggsurvplot(km_fit, data = Immuno_NMDeff_df, risk.table = TRUE,
                  pval = TRUE, conf.int = FALSE,
                  ylim = c(0, 1), #xlim = c(0, 22), 
                  xlab = "Time (days)", ylab = "Survival probability",
                  ggtheme = theme_bw(), surv.median.line = "v", legend.labs = c("High", "Low"), 
                  legend.title = paste0("iNMDeff")) + ggtitle(paste0(NMD_method))
        list_plots[[paste0("TCGA_",NMD_method)]] <- p$plot
      }
  }
  if (!is.null(KM_output)) {
    return(list_KM_objects)
  } else {
    return(list_plots)
  }
  
}

# Function to calculate TPM
# calculateTPM <- function(counts, lengths) {
#   # Convert counts to reads per kilobase (RPK)
#   rpk <- sweep(counts, 1, lengths, FUN="/")
#   # Calculate per million scaling factor
#   per_million_scale_factor <- colSums(rpk) / 1e6
#   # Divide the RPK values by the scale factor
#   tpm <- sweep(rpk, 2, per_million_scale_factor, FUN="/")
#   return(tpm)
# }

# conda activate iNMDeff_models
library(dplyr)
library(ggplot2)
library(ggpubr)
library(survival)
library(survminer)
library(caret)
library(preprocessCore)
library(glmnet)
library(readxl)
library(sva)
library(ggrepel)
library(ggpmisc)
library(cowplot)
# Normalize matrix to TPM
library(DGEobj.utils)
library(DESeq2)
# UCSC to Ensembl
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
# logistic regression
library(tidyr)
library(pROC)
# Fisher test to combine p-values
# library(survcomp)
library(MetaKTSP)
library(multiclassPairs)
library(switchBox)

# 1) Data
endogenous_NMD_genesets <- c("endogenous_NMD_Colombo","endogenous_NMD_Karousis","endogenous_NMD_Tani","endogenous_NMD_Courtney","endogenous_NMD_ensembl",
                      "endogenous_NMD_all","endogenous_NMD_Consensus","endogenous_SMG6","endogenous_SMG7",
                      "endogenous_non_NMD_neg_control","endogenous_non_NMD_neg_control_with_NMD_features")
ASE_NMD_genesets <- c("ASE_PTC_NMD_triggering_0.01","ASE_PTC_NMD_evading_0.01","ASE_synonymous_0.01",
                      "ASE_PTC_NMD_triggering_0.2","ASE_PTC_NMD_evading_0.2","ASE_synonymous_0.2")
# 1.1) samples NMDeff
# 1.1.1) TCGA
sample_NMD_efficiencies_TCGA_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/cancers/pancancer/NMD_efficiencies_TCGA.txt"
sample_NMD_efficiencies_TCGA <- read.table(file = sample_NMD_efficiencies_TCGA_path, header = TRUE, sep = "\t")
sample_NMD_efficiencies_TCGA <- modify_NMDeff_dataframe(sample_NMD_efficiencies_TCGA, dataset = "TCGA", scale = FALSE)
# 1.1.2) GTEx
sample_NMD_efficiencies_GTEx_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/GTEx/pantissue/NMD_efficiencies_GTEx.txt"
sample_NMD_efficiencies_GTEx <- read.table(file = sample_NMD_efficiencies_GTEx_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
sample_NMD_efficiencies_GTEx <- modify_NMDeff_dataframe(sample_NMD_efficiencies_GTEx, dataset = "GTEx", scale = FALSE)

# 1.2) TCGA metadata
TCGA_metadata <- read.table("/g/strcombio/fsupek_cancer1/gpalou/TCGA_metadata/TCGA_clinical_all.tsv", 
                            header = TRUE, sep = "\t", stringsAsFactors = FALSE)
TCGA_metadata <- TCGA_metadata[,colnames(TCGA_metadata) %in% c("submitter_id", "gender", "race", "age_at_diagnosis", "vital_status", "days_to_death", "days_to_last_follow_up")]
colnames(TCGA_metadata) <- c("sample", "gender", "race", "age_at_diagnosis", "vital_status", "days_to_death", "days_to_last_follow_up")

# 1.3) TCGA-SKCM immunotherapy data
# df <- GDCquery_clinic(project = "TCGA-SKCM", type = "clinical", save.csv = FALSE)
# df$treatments_radiation_treatment_type
# grep("Radiation Th",df)
TCGA_SKCM_metadata <- read.table("/g/strcombio/fsupek_cancer1/gpalou/TCGA_clinical/TCGA_SKCM_immunotherapy/nationwidechildrens.org_clinical_drug_skcm.txt", 
                            header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# Merge and add survival variables
sample_NMD_efficiencies_TCGA[,c("vital_status","days_to_death","days_to_last_follow_up")] <- NULL
TCGA_NMDeff <- merge(sample_NMD_efficiencies_TCGA, TCGA_metadata, all.x = TRUE)
TCGA_NMDeff[TCGA_NMDeff$vital_status == "dead",]$days_to_last_follow_up <- TCGA_NMDeff[TCGA_NMDeff$vital_status == "dead",]$days_to_death
TCGA_NMDeff$status <- NA
TCGA_NMDeff[TCGA_NMDeff$vital_status == "dead",]$status <- 1
TCGA_NMDeff[TCGA_NMDeff$vital_status == "alive",]$status <- 0
TCGA_NMDeff$days_to_last_follow_up <- as.numeric(TCGA_NMDeff$days_to_last_follow_up)/365.25
TCGA_NMDeff$age_at_diagnosis <- as.numeric(TCGA_NMDeff$age_at_diagnosis)/365.25

# 1.4) RNA-seq from Immunotherapy clinical trials
# 1.4.1) # TCGA
# Gene-level TPM RNAseq data
output_path <- "/g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/primary_tumor/TCGA_RNAseq_matrix_TPM_gene.txt"
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

# 1.4.2) GTEx 
# Gene-level TPM RNAseq data. Not Immunotherapy but normal tissue
GTEx_path <- "/g/strcombio/fsupek_cancer1/gpalou/GTEx/RNAseq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct"
RNAseq_GTEx_TPM_all <- read.table(file = GTEx_path, header = TRUE, sep = "\t", row.names = 1, skip = 2)
print("Dimensions -->")
print(dim(RNAseq_GTEx_TPM_all))
RNAseq_GTEx_TPM_all <- RNAseq_GTEx_TPM_all[-grep("PAR", rownames(RNAseq_GTEx_TPM_all)), ]
rownames(RNAseq_GTEx_TPM_all) <- gsub("(.*)\\..*", "\\1", rownames(RNAseq_GTEx_TPM_all))
RNAseq_GTEx_TPM_all$Description <- NULL

# 1.4.3) Liu 2019 --> Melanoma. 
# RNA-seq gene-level TPM. Immunotherapy clinical trial
Liu_SKCM <- read.table(file = "/g/strcombio/fsupek_cancer1/gpalou/RNAseq_immunotherapy_trials/Liu_2019_melanoma/41591_2019_654_MOESM3_ESM.txt", 
              header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)
Liu_SKCM <- log2(Liu_SKCM+1)
Liu_SKCM <- data.frame(t(Liu_SKCM))
# Metadata
Liu_SKCM_metadata <- read_excel("/g/strcombio/fsupek_cancer1/gpalou/RNAseq_immunotherapy_trials/Liu_2019_melanoma/41591_2019_654_MOESM4_ESM.xlsx",
            skip = 2, sheet = "Supplemental Table 1")
colnames(Liu_SKCM_metadata)[1] <- "sampleID"

# 1.4.4) Braun 2020 --> RCC.
# RNA-seq gene-level. It is quantile normalized, log2(TPM+1), batch corrected with ComBat, etc...
# DO NOT USE THIS!!
# Braun_RCC <- read.table(file = "/g/strcombio/fsupek_cancer1/gpalou/RNAseq_immunotherapy_trials/Braun_2020_RCC/gene_expression.tsv", 
#               header = TRUE, sep = ",", stringsAsFactors = FALSE, row.names = 1)
# # Metadata
# Braun_RCC_metadata <- read_excel("/g/strcombio/fsupek_cancer1/gpalou/RNAseq_immunotherapy_trials/Braun_2020_RCC/41591_2020_839_MOESM2_ESM.xlsx",
#             skip = 1, sheet = "S1_Clinical_and_Immune_Data")
# colnames(Braun_RCC_metadata)[1] <- "sampleID"

# 1.4.5) Carrol 2023 --> EAC
# RNA-seq gene-level raw STAR counts
Carrol_EAC_raw <- read.table(file = "/g/strcombio/fsupek_cancer1/gpalou/RNAseq_immunotherapy_trials/Carrol_2023_EAC/licroxford-carroll_etal_2023-5fdcb6072294/processed_data/bulk_RNAseq_counts/LUD2015-005_RNAseq_featureCounts.tsv", 
              header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)
Carrol_EAC <- Carrol_EAC_raw[,!colnames(Carrol_EAC_raw) %in% c("Chr","Start","End","Strand","Length")]
# Normalize counts to Log2(TPM+1)
Carrol_EAC <- convertCounts(countsMatrix = as.matrix(Carrol_EAC), geneLength = Carrol_EAC_raw$Length, 
                              unit = "TPM", prior.count = 1, log = TRUE )
# Metadata
Carrol_EAC_metadata <- read_excel("/g/strcombio/fsupek_cancer1/gpalou/RNAseq_immunotherapy_trials/Carrol_2023_EAC/mmc7.xlsx",
            skip = 0, sheet = "Sheet1")
colnames(Carrol_EAC_metadata)[1] <- "sampleID"

# 1.4.6) Motzer 2020 --> RCC (JAVELIN Renal 101)
# RNA-seq gene-level (log2 TPM, values < 0.01 set to 0.01 )
Motzer_RCC <- read_excel("/g/strcombio/fsupek_cancer1/gpalou/RNAseq_immunotherapy_trials/Motzer_2020_RCC/41591_2020_1044_MOESM3_ESM.xlsx",
                                    skip = 1, sheet = "S13_Gene_expression_TPM")
Motzer_RCC <- data.frame(Motzer_RCC[!duplicated(Motzer_RCC$HUGO),])
rownames(Motzer_RCC) <- Motzer_RCC$HUGO
Motzer_RCC$HUGO <- NULL
# Motzer_RCC <- data.frame(t(Motzer_RCC))
# Metadata
Motzer_RCC_metadata <- read_excel("/g/strcombio/fsupek_cancer1/gpalou/RNAseq_immunotherapy_trials/Motzer_2020_RCC/41591_2020_1044_MOESM3_ESM.xlsx",
                                    skip = 1, sheet = "S11_Clinical_data")

# 1.4.7) Motzer 2020 --> RCC (IMMOTION151)

# 1.4.8) Riaz 2017 --> Melanoma
# RNA-seq gene-level raw counts (UCSC hg19)
Riaz_SKCM <- read.csv("/g/strcombio/fsupek_cancer1/gpalou/RNAseq_immunotherapy_trials/Riaz_2017_melanoma/GSE91061_BMS038109Sample.hg19KnownGene.raw.csv")
# Metadata
Riaz_SKCM_metadata <- read_excel("/g/strcombio/fsupek_cancer1/gpalou/RNAseq_immunotherapy_trials/Riaz_2017_melanoma/NIHMS907788-supplement-9.xlsx",
                                    skip = 2, sheet = "Table S2")

# 1.4.9) Hartwig
files = list.files("/g/strcombio/fsupek_decider/msalvadores/data_isofox", pattern = "gene_data", full.names = T, recursive = T)
sample_id = gsub(".*data_isofox\\/[A-Z0-9].*\\/([A-Z0-9].*)\\.isf\\.gene\\_data\\.csv","\\1",files)
RNAseq_Hartwig_TPM_all <- c()
c <- 0
for (sample in sample_id) {
  c <- c+1
  print( round((c/length(sample_id))* 100,2))
  RNAseq_sample_path <- paste0("/g/strcombio/fsupek_decider/msalvadores/data_isofox/",sample,"/",sample,".isf.gene_data.csv")
  RNAseq_sample <- read.csv(file = RNAseq_sample_path, header = TRUE, sep = ",", stringsAsFactors = FALSE)
  RNAseq_Hartwig_TPM_sample <- RNAseq_sample[,c("GeneId","RawTPM")]
  rownames(RNAseq_Hartwig_TPM_sample) <- RNAseq_Hartwig_TPM_sample$GeneId
  RNAseq_Hartwig_TPM_sample$GeneId <- NULL
  colnames(RNAseq_Hartwig_TPM_sample)[1] <- sample
  if (length(RNAseq_Hartwig_TPM_all) == 0) {
    RNAseq_Hartwig_TPM_all <- RNAseq_Hartwig_TPM_sample
  } else {
    RNAseq_Hartwig_TPM_all <- merge(RNAseq_Hartwig_TPM_all,RNAseq_Hartwig_TPM_sample, by.x = "row.names", by.y = "row.names")
    rownames(RNAseq_Hartwig_TPM_all) <- RNAseq_Hartwig_TPM_all$Row.names
    RNAseq_Hartwig_TPM_all$Row.names <- NULL
  }
}
print(dim(RNAseq_Hartwig_TPM_all))
RNAseq_Hartwig_TPM_all <- log2(RNAseq_Hartwig_TPM_all+1)
RNAseq_Hartwig_TPM_all <- data.frame(RNAseq_Hartwig_TPM_all)
out_path <- "/g/strcombio/fsupek_cancer1/gpalou/Hartwig/RNAseq_Hartwig_TPM_all.txt"
# write.table(RNAseq_Hartwig_TPM_all, file = out_path, sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
RNAseq_Hartwig_TPM_all <- read.table(file = out_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# Metadata
Hartwig_metadata <- read.table(file = "/g/strcombio/fsupek_cancer1/HARTWIG_vcfs/somatic_update_Nov2021/metadata.tsv",
                                    header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# Metadata 2
Hartwig_metadata_2 <- read_excel("/g/strcombio/fsupek_cancer1/gpalou/Hartwig/Hartwig_metadata_2.xlsx")
Hartwig_metadata_2 <- data.frame(Hartwig_metadata_2)
# 1.5) ENSEMBL gene ID - Gene Symbol conversor table
ensembl_gene_symbol <- read.table(file = "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/ensembl_v88_gene_transcript_genesymbol.txt",
                                    header = TRUE, sep = "\t", stringsAsFactors = FALSE)
ensembl_gene_symbol <- unique(ensembl_gene_symbol[,c(1:2)])
ensembl_gene_symbol <- ensembl_gene_symbol[!duplicated(ensembl_gene_symbol$gene_name),]

# 1.6) TCGA drug therapy

TCGA_cancer_metadata_drug_therapy_all <- c()
cancer_types_original <- as.character(unique(sample_NMD_efficiencies_TCGA$cancer_type))
for (cancer in cancer_types_original) {
  print(cancer)
  cancer <- tolower(cancer)
  # Drug therapy metadata
  if (cancer != "laml") {
    TCGA_cancer_metadata_drug_therapy <- read.csv(paste0("/g/strcombio/fsupek_cancer1/gpalou/TCGA_metadata/TCGA_gdc_clinical/nationwidechildrens.org_clinical_drug_",cancer,".txt"), 
                                header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    cols_keep <- c("bcr_patient_barcode","pharmaceutical_therapy_drug_name","pharmaceutical_therapy_type")
    if ("treatment_best_response" %in% colnames(TCGA_cancer_metadata_drug_therapy)) {
      cols_keep <- c(cols_keep,"treatment_best_response")
      TCGA_cancer_metadata_drug_therapy <- TCGA_cancer_metadata_drug_therapy[-1,cols_keep]
    } else {
      TCGA_cancer_metadata_drug_therapy <- TCGA_cancer_metadata_drug_therapy[-1,cols_keep]
      TCGA_cancer_metadata_drug_therapy$treatment_best_response <- NA
    }
    colnames(TCGA_cancer_metadata_drug_therapy)[colnames(TCGA_cancer_metadata_drug_therapy) %in% "bcr_patient_barcode"] <- "sample"
    if (length(TCGA_cancer_metadata_drug_therapy_all) == 0) {
      TCGA_cancer_metadata_drug_therapy_all <- TCGA_cancer_metadata_drug_therapy
    } else {
      TCGA_cancer_metadata_drug_therapy_all <- rbind(TCGA_cancer_metadata_drug_therapy_all,TCGA_cancer_metadata_drug_therapy)
    }
  }
  if (cancer == "coad") {next} else if (cancer == "read") { cancer <- "coadread"}
}

# 2) Convert to ensembl gene id
# 2.1) Liu 2019
# Convert gene symbols to ensembl gene id
Liu_SKCM <- merge(Liu_SKCM,ensembl_gene_symbol, by.x = "row.names", by.y = "gene_name", all.x = TRUE, sort = FALSE)
Liu_SKCM$Row.names <- NULL
Liu_SKCM <- na.omit(Liu_SKCM)
rownames(Liu_SKCM) <- Liu_SKCM$gene_id

# 2.2) Braun 2020
# Braun_RCC <- merge(Braun_RCC,ensembl_gene_symbol, by.x = "row.names", by.y = "gene_name", all.x = TRUE, sort = FALSE)
# Braun_RCC$Row.names <- NULL
# Braun_RCC <- na.omit(Braun_RCC)
# rownames(Braun_RCC) <- Braun_RCC$gene_id

# 2.3) Carrol 2023
# Carrol_EAC <- merge(Carrol_EAC,ensembl_gene_symbol, by.x = "row.names", by.y = "gene_name", all.x = TRUE, sort = FALSE)
# Carrol_EAC$Row.names <- NULL
# Carrol_EAC <- na.omit(Carrol_EAC)
# rownames(Carrol_EAC) <- Carrol_EAC$gene_id

# 2.4) Motzer 2020
Motzer_RCC <- merge(Motzer_RCC,ensembl_gene_symbol, by.x = "row.names", by.y = "gene_name", all.x = TRUE, sort = FALSE)
Motzer_RCC$Row.names <- NULL
Motzer_RCC <- na.omit(Motzer_RCC)
rownames(Motzer_RCC) <- Motzer_RCC$gene_id

# 2.5) Riaz 2017
# Extract all the transcript ids
txid = keys(TxDb.Hsapiens.UCSC.hg19.knownGene, "TXID")
# and get their corresponding Entrez gene ids
TxGeneID = select(TxDb.Hsapiens.UCSC.hg19.knownGene, txid, "GENEID", "TXID")
# Gene Symbol
TxGeneNameID <- select(org.Hs.eg.db, unique(TxGeneID$GENEID), c("SYMBOL", "GENENAME","ENSEMBL"))
# Add Ensembl Gene ID (ours)
TxGeneNameID <- merge(TxGeneNameID,ensembl_gene_symbol, by.x = "SYMBOL", by.y = "gene_name", all.x = TRUE)
Riaz_SKCM <- merge(Riaz_SKCM,TxGeneNameID, by.x = "X", by.y = "ENTREZID", all.x = TRUE)
Riaz_SKCM <- Riaz_SKCM[!duplicated(Riaz_SKCM$gene_id),]
Riaz_SKCM <- Riaz_SKCM[!is.na(Riaz_SKCM$gene_id),]
rownames(Riaz_SKCM) <- Riaz_SKCM$gene_id
Riaz_SKCM[,c("X","SYMBOL","GENENAME","ENSEMBL")] <- NULL

# Normalize counts to Log2(TPM+1)
Riaz_SKCM <- merge(Riaz_SKCM,Carrol_EAC_raw[,"Length",drop=FALSE], by.x = "row.names", by.y  = "row.names", all.x=TRUE)
gene_lengths <- Riaz_SKCM$Length
row.names(Riaz_SKCM) <- Riaz_SKCM$Row.names
Riaz_SKCM[,c("Row.names","Length","gene_id")] <- NULL
Riaz_SKCM <- convertCounts(countsMatrix = as.matrix(Riaz_SKCM), geneLength = gene_lengths, 
                              unit = "TPM", prior.count = 1, log = TRUE )
# Remove NAs
to_remove <- which(is.na(Riaz_SKCM[,1]))
Riaz_SKCM <- Riaz_SKCM[-to_remove,]

# 3) Shared genes

# Liu_SKCM_genes <- colnames(Liu_SKCM)
# TCGA_genes <- rownames(RNAseq_TCGA_TPM_all)
# GTEx_genes <- rownames(RNAseq_GTEx_TPM_all)
# all_shared_genes <- intersect(intersect(Liu_SKCM_genes, TCGA_genes), GTEx_genes)
# all_shared_genes <- intersect(Liu_SKCM_genes, TCGA_genes)

# 4) NMDeff High/Low vs Immunotherapy response

# Observation --> SKCM patients have low survival when NMD-eff is high
# Hypothesis --> Those patients with NMD-eff High might be worse responders to immunotherapy as well

############################
#### 4.1) TCGA --> SKCM ####
############################

TCGA_df <- TCGA_NMDeff[TCGA_NMDeff$cancer_type == "SKCM",]
TCGA_df <- merge(TCGA_df,TCGA_cancer_metadata_drug_therapy_all, by.x = "sample", by.y = "sample")

TCGA_therapy_counts <- TCGA_df %>%
      dplyr::group_by(sample) %>%
      dplyr::summarize(immunotherapy_counts = sum(pharmaceutical_therapy_type == "Immunotherapy"),
              chemotherapy_counts = sum(pharmaceutical_therapy_type == "Chemotherapy"),
              other_therapy_counts = sum(!pharmaceutical_therapy_type %in% c("Immunotherapy","Chemotherapy")))

samples_to_keep <- as.character(TCGA_therapy_counts[TCGA_therapy_counts$immunotherapy_counts >= 1,"sample"]$sample)
TCGA_df <- TCGA_df %>%
    filter(sample %in% samples_to_keep) %>%
    filter(pharmaceutical_therapy_type %in% c("Immunotherapy")) #%>%# grep("ivolumab",TCGA_df$pharmaceutical_therapy_drug_name)

# Each sample can be treated with >1 immunotherapy drug
# Let's go one by one and choose the best treatment
TCGA_df_treatment <- c()
for (TCGA_sample in unique(TCGA_df$sample)) {
  TCGA_df_tmp <- TCGA_df %>%
    filter(sample %in% TCGA_sample)
  response <- TCGA_df_tmp$treatment_best_response
  if (length(response) == 1) {
    if (length(TCGA_df_treatment) == 0) {
      TCGA_df_treatment <- TCGA_df_tmp
    } else {
      TCGA_df_treatment <- rbind(TCGA_df_treatment,TCGA_df_tmp)
    }
  } else {
    # If NA
    NA_sum <- sum(response %in% c("[Not Applicable]","[Not Available]","[Unknown]")) == 2
    response_sum <- sum(table(response))
    if (isTRUE(NA_sum) | response_sum == 2) {
      # Take first
      TCGA_df_tmp <- TCGA_df_tmp[1,]
    } else {
      TCGA_df_tmp <- TCGA_df_tmp[!TCGA_df_tmp$treatment_best_response %in% c("[Not Applicable]","[Not Available]","[Unknown]"),]
    }
    if (nrow(TCGA_df_tmp) == 1) {
      TCGA_df_treatment <- rbind(TCGA_df_treatment,TCGA_df_tmp)
    } else if (nrow(TCGA_df_tmp) == 2) {
      TCGA_df_tmp <- TCGA_df_tmp[TCGA_df_tmp$treatment_best_response %in% c("Complete Response","Partial Response","Stable Disease"),]
      if (nrow(TCGA_df_tmp) == 1) {
        TCGA_df_treatment <- rbind(TCGA_df_treatment,TCGA_df_tmp)
      } else {
        print(TCGA_sample)
      }
    }
  }
}

# ipilimumab == anti CTLA-4 --> metastatic melanoma
# pembrolizumab and nivolumab (anti PD-1) --> melanoma and KIRC, LUSC

# Add metadata
TCGA_df_treatment$treatment_response <- TCGA_df_treatment$treatment_best_response 
TCGA_df_treatment$treatment_response <- ifelse(TCGA_df_treatment$treatment_response == "Clinical Progressive Disease" | TCGA_df_treatment$treatment_response == "Stable Disease","Non-Responders",TCGA_df_treatment$treatment_response)
TCGA_df_treatment$treatment_response <- ifelse(TCGA_df_treatment$treatment_response == "Partial Response" | TCGA_df_treatment$treatment_response == "Complete Response","Responders",TCGA_df_treatment$treatment_response)
TCGA_df_treatment$treatment_response <- ifelse(TCGA_df_treatment$treatment_response == "[Not Applicable]" | TCGA_df_treatment$treatment_response == "[Not Available]" | TCGA_df_treatment$treatment_response == "[Unknown]",
                                    NA,TCGA_df_treatment$treatment_response)
table(TCGA_df_treatment$treatment_response)

# out_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/phenotypes/immunotherapy/Immuno_NMDeff_TCGA_SKCM_results.txt")
# Immuno_NMDeff_all <- read.table(file = out_path, header = TRUE, sep = "\t")
# write.table(TCGA_df_treatment, file = out_path, sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)

# A) NMDeff vs Treatment Response

# TCGA_df$immunotherapy_response <- ifelse(TCGA_df$treatment_best_response == "Disease",1,0)
TCGA_df_stacked <- stack(TCGA_df_treatment[,c("endogenous_NMD_Consensus","ASE_PTC_NMD_triggering_0.2")])
# List of columns to repeat
columns_to_repeat <- c("treatment_response", "status", "TMB", "days_to_last_follow_up", "age_at_diagnosis", 
                        "gender", "ethnicity","CNV_burden", "cancer_subtype")
# Automatically repeat each column from TCGA_df_treatment and assign to TCGA_df_stacked
for (col in columns_to_repeat) {
  TCGA_df_stacked[[col]] <- rep(TCGA_df_treatment[[col]],2)
}
colnames(TCGA_df_stacked) <- c("NMDeff","NMD_method",columns_to_repeat)

# TCGA_df_stacked$treatment_response <- rep(TCGA_df_treatment$treatment_response,2)
# TCGA_df_stacked$status <- rep(TCGA_df_treatment$status,2)
# TCGA_df_stacked$TMB <- rep(TCGA_df_treatment$TMB,2)
# TCGA_df_stacked$days_to_last_follow_up <- rep(TCGA_df_treatment$days_to_last_follow_up,2)

TCGA_SKCM_imm_resp_model <- TCGA_df_stacked
# N of samples with iNMDeff
sum(!is.na(TCGA_SKCM_imm_resp_model[TCGA_SKCM_imm_resp_model$NMD_method == "endogenous_NMD_Consensus","NMDeff"]))
sum(!is.na(TCGA_SKCM_imm_resp_model[TCGA_SKCM_imm_resp_model$NMD_method == "ASE_PTC_NMD_triggering_0.2","NMDeff"]))

combinations <- combn(names(table(TCGA_df_treatment$treatment_response)), 2, simplify = FALSE)

p <- TCGA_df_stacked %>%  
        ggplot(aes(x = treatment_response, y = NMDeff, fill = factor(treatment_response))) +
          labs(x = "", y = "NMDeff", fill = "") +
          ggtitle("TCGA - SKCM") +
          theme_classic(base_size = 25)+
          facet_wrap(NMD_method ~ ., scale = "free_y") +
          geom_violin() + scale_fill_brewer(palette = "Dark2") + xlab("") +
          geom_boxplot(aes(fill = factor(treatment_response)), width=0.2, color="black", alpha=0.2) +
          geom_jitter(aes(fill = factor(treatment_response)), 
                      position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0, dodge.width = 0.9),
                      alpha = 0.25, size = 1) +        theme_bw(base_size = 25) +
          theme(
              axis.text.x = element_text(size = 18,hjust = 1), 
              legend.position = "top",
              plot.title = element_text(size = 35, hjust = 0.5)) +
          stat_compare_means(comparisons = combinations, size = 8,
                            #label.y = c(1,1.5,2),
                            label = "p.format", method = "wilcox.test", hide.ns = TRUE)# +

NMDeff_plot <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/phenotypes/immunotherapy/TCGA/TCGA_SKCM_NMDeff_vs_treatment_response_boxplots.png")
png(NMDeff_plot, width = 3500, height = 2500, res = 300)
print(p)
dev.off()

#a <- TCGA_SKCM_metadata[grep("[Pp]embrolizumab|[Nn]ivolumab|[iI]pilimumab|[Tt]remelimumab",TCGA_SKCM_metadata$pharmaceutical_therapy_drug_name),]

# B) Scatterplot NMDeff vs PFS

p <- TCGA_df_stacked %>%  
      ggplot(aes(x = days_to_last_follow_up, y = NMDeff)) +
        labs(x = "days_to_last_follow_up", y = "NMDeff") +
        ggtitle("TCGA - SKCM") +
        theme_classic(base_size = 25)+
        facet_wrap(NMD_method ~ ., scale = "free_y") +
          geom_point()+
          geom_smooth(method = "lm", se = TRUE, size = 1) +
          stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, size = 5) +
        theme(
            axis.text.x = element_text(size = 18,hjust = 1), 
            legend.position = "top",
            plot.title = element_text(size = 35, hjust = 0.5))
output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/phenotypes/immunotherapy/TCGA/TCGA_SKCM_NMDeff_vs_PFS_scatterplot.png")
png(output_path, width = 3500, height = 2000, res = 300)
print(p)
dev.off()

# C) Survival curve

TCGA_df_stacked$dataset <- "TCGA"
KM_plots <- survival_curve(Immuno_NMDeff = TCGA_df_stacked, dataset_to_predict = "TCGA", percentile = 30 )
KM_plot_final <- plot_grid(KM_plots[[1]],KM_plots[[2]])
KM_plot <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/phenotypes/immunotherapy/TCGA/TCGA_SKCM_NMDeff_survival_curve.png")
png(KM_plot, width = 3500, height = 1500, res = 300)
print(KM_plot_final)
dev.off()

####################################
#### 4.2) Liu 2019 --> Melanoma ####
####################################

Immuno_NMDeff_all <- c()
for (dataset in c("TCGA","GTEx")) {
  for (NMD_method in c("endogenous_NMD_Consensus","ASE_PTC_NMD_triggering_0.2")) {
      Immuno_NMDeff <- NMDeff_gene_model_validation(dataset = dataset, model = "LASSO", 
                          target_var = NMD_method, dataset_to_predict = "Liu_SKCM")
      Immuno_NMDeff$dataset <- dataset
      Immuno_NMDeff$NMD_method <- NMD_method
      if (length(Immuno_NMDeff_all) == 0) {
        Immuno_NMDeff_all <- Immuno_NMDeff
      } else {
        Immuno_NMDeff_all <- rbind(Immuno_NMDeff_all,Immuno_NMDeff)
      }
  }
}
Immuno_NMDeff_all$treatment_response <- factor(Immuno_NMDeff_all$treatment_response, levels = c("CR","PR","PD","SD"))
out_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/phenotypes/immunotherapy/Immuno_NMDeff_Liu_SKCM_results.txt")
Immuno_NMDeff_all <- read.table(file = out_path, header = TRUE, sep = "\t")
# write.table(Immuno_NMDeff_all, file = out_path, sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)

# Immuno_NMDeff_all %>%
#       group_by(dataset, NMD_method) %>%
#         summarize(
#           cor_coef = cor.test(ploidy, NMDeff)$estimate,
#           p_value = cor.test(ploidy, NMDeff)$p.value
#         )

# Scale
# Immuno_NMDeff_all <- Immuno_NMDeff_all %>%
#                 group_by(dataset,NMD_method) %>%
#                 mutate(NMDeff = scale(NMDeff))
# Mean of four models
Immuno_NMDeff_all <- Immuno_NMDeff_all %>%
                group_by(Row.names, NMD_method) %>%
                mutate(NMDeff_mean = mean(NMDeff, na.rm = TRUE))
# Duplicated (mean)
rows <- which(duplicated(Immuno_NMDeff_all[,c("Row.names","NMDeff_mean","NMD_method")]))
Immuno_NMDeff_all <- Immuno_NMDeff_all[-rows,]

# N of samples with iNMDeff
sum(!is.na(Immuno_NMDeff_all[Immuno_NMDeff_all$NMD_method == "endogenous_NMD_Consensus","NMDeff"]))
sum(!is.na(Immuno_NMDeff_all[Immuno_NMDeff_all$NMD_method == "ASE_PTC_NMD_triggering_0.2","NMDeff"]))

# A) Add iNMDeff binary predictions by Random Forest (multiclassPairs model)

# Test with different models (TCGA or GTex)
proxy_model_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/proxy_iNMDeff_models/tissues"
path <- paste0(proxy_model_path,"/SKCM_iNMDeff_RF_classifier.RData")
# path <- paste0(proxy_model_path,"/Skin_Not_Sun_Exposed_Suprapubic_iNMDeff_RF_classifier.RData")
# path <- paste0(proxy_model_path,"/Skin_Sun_Exposed_Lower_leg_iNMDeff_RF_classifier.RData")

iNMDeff_RF_classifier <- readRDS(path)

# Apply on validation datasets (gene expression matrix)
Immuno_dataset <- Liu_SKCM
results <- predict_RF(classifier = iNMDeff_RF_classifier, 
                      Data = Immuno_dataset,
                      impute = TRUE) # can handle missed genes by imputation
results$predictions
iNMDeff_class_res <- data.frame(results$predictions)
colnames(iNMDeff_class_res) <- c("NMD_High_prob","NMD_Low_prob")
#Just to check which is the best model
iNMDeff_class_res$NMD_classification <- case_when(
  iNMDeff_class_res$NMD_High_prob >= 0.65 ~ "NMD-High",
  iNMDeff_class_res$NMD_High_prob <= 0.35 ~ "NMD-Low",
  TRUE ~ "NMD-Mid"
) 
table(iNMDeff_class_res$NMD_classification)

Immuno_NMDeff_all <- merge(Immuno_NMDeff_all,iNMDeff_class_res, 
                by.x = "Row.names", by.y = "row.names",
                all.x = TRUE)

# B) NMDeff vs Treatment Response

# Immuno_NMDeff_all$treatment_response <- ifelse(Immuno_NMDeff_all$treatment_response %in% c("CR","PD"), as.character(Immuno_NMDeff_all$treatment_response),NA)
# Immuno_NMDeff_all$treatment_response <- ifelse(Immuno_NMDeff_all$treatment_response %in% c("SD"),NA,as.character(Immuno_NMDeff_all$treatment_response))
Immuno_NMDeff_all$treatment_response <- ifelse(Immuno_NMDeff_all$treatment_response %in% c("MR"),NA,as.character(Immuno_NMDeff_all$treatment_response))
Immuno_NMDeff_all$treatment_response <- ifelse(Immuno_NMDeff_all$treatment_response %in% c("CR","PR"), "Responders",as.character(Immuno_NMDeff_all$treatment_response))
Immuno_NMDeff_all$treatment_response <- ifelse(Immuno_NMDeff_all$treatment_response %in% c("PD","SD"), "Non-responders",as.character(Immuno_NMDeff_all$treatment_response))
Immuno_NMDeff_all_filt <- Immuno_NMDeff_all[!is.na(Immuno_NMDeff_all$treatment_response),]
table(Immuno_NMDeff_all_filt$treatment_response)

Liu_SKCM_imm_resp_model <- Immuno_NMDeff_all_filt

# Immuno_NMDeff_all$treatment_response <- ifelse(Immuno_NMDeff_all$treatment_response %in% c("CR","PR","SD"), "Responders","Non-responders")
# Immuno_NMDeff_all$treatment_response <- ifelse(Immuno_NMDeff_all$treatment_response %in% c("CR","PD"), as.character(Immuno_NMDeff_all$treatment_response),NA)
# Immuno_NMDeff_all$treatment_response <- ifelse(Immuno_NMDeff_all$treatment_response %in% c("CR"), "Responders","Non-responders")
combinations <- combn(names(table(Immuno_NMDeff_all_filt$treatment_response)), 2, simplify = FALSE)

# 2) Plot of iNMDeff vs treatment response
p <- Immuno_NMDeff_all_filt %>%  
      # filter(total_muts <= 1500 & CNA_prop <= 0.55 & purity >= 0.3 & ploidy <= 3 & Primary_Type != "occult") %>%
      filter(!is.na(treatment_response)) %>%
      filter(treatment_response != "MR") %>%
        ggplot(aes(x = treatment_response, y = NMDeff_mean, fill = factor(treatment_response))) +
          labs(x = "", y = "NMDeff", fill = "") +
          ggtitle("Liu 2019 - Melanoma") +
          theme_classic(base_size = 25)+
          facet_wrap(NMD_method ~ ., scale = "free_y") +
          geom_violin() + scale_fill_brewer(palette = "Dark2") + xlab("") +
          geom_boxplot(aes(fill = factor(treatment_response)), width=0.2, color="black", alpha=0.2) +
          geom_jitter(aes(fill = factor(treatment_response)), 
                      position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0, dodge.width = 0.9),
                      alpha = 0.25, size = 1) +        theme_bw(base_size = 25) +
          theme(
              axis.text.x = element_text(size = 18,hjust = 1), 
              legend.position = "top",
              plot.title = element_text(size = 35, hjust = 0.5)) +
          stat_compare_means(comparisons = combinations, size = 8,
                            #label.y = c(1,1.5,2),
                            label = "p.format", method = "wilcox.test", hide.ns = TRUE)# +
output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/phenotypes/immunotherapy/validation_datasets/Liu_SKCM/Liu_SKCM_NMDeff_vs_treatment_response_boxplots.png")
png(output_path, width = 4500, height = 3500, res = 300)
print(p)
dev.off()

# C) Scatterplot NMDeff vs PFS
p <- Immuno_NMDeff_all %>%  
      ggplot(aes(x = PFS, y = NMDeff_mean)) +
        labs(x = "PFS", y = "NMDeff") +
        ggtitle("Liu 2019 - Melanoma") +
        theme_classic(base_size = 25)+
        # facet_wrap(dataset + NMD_method ~ ., scale = "free_y") +
        facet_wrap(NMD_method ~ ., scale = "free_y") +
          geom_point()+
          geom_smooth(method = "lm", se = TRUE, size = 1) +
          stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, size = 5) +
        theme(
            axis.text.x = element_text(size = 18,hjust = 1), 
            legend.position = "top",
            plot.title = element_text(size = 35, hjust = 0.5))
output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/phenotypes/immunotherapy/validation_datasets/Liu_SKCM/Liu_SKCM_NMDeff_vs_PFS_scatterplot.png")
png(output_path, width = 4500, height = 3500, res = 300)
print(p)
dev.off()

# D) Survival Curve
Immuno_NMDeff_all$PFS <- as.numeric(Immuno_NMDeff_all$PFS)/365.25

# Continous iNMDeff classification by percentiles

for (percentile in c(5,10,15,20,25,30,35,40,45,50)) {
  Immuno_NMDeff_filt <- Immuno_NMDeff_all %>% filter(numPriorTherapies <= 1)
  KM_plots <- survival_curve(Immuno_NMDeff = Immuno_NMDeff_filt, dataset_to_predict = "Liu_SKCM", percentile = percentile,
                      KM_output = NULL, survival_var = "PFS", iNMDeff_model_type = "continous", binary_probability = NULL)
  # KM_plot_final <- plot_grid(KM_plots[[1]],KM_plots[[2]],KM_plots[[3]],KM_plots[[4]])
  KM_plot_final <- plot_grid(KM_plots[[1]],KM_plots[[2]])
  KM_plot <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/phenotypes/immunotherapy/validation_datasets/Liu_SKCM/Liu_SKCM_NMDeff_PFS_survival_curve_percentile_",percentile,".png")
  png(KM_plot, width = 3500, height = 2500, res = 300)
  print(KM_plot_final)
  dev.off()
  KM_plots <- survival_curve(Immuno_NMDeff = Immuno_NMDeff_filt, dataset_to_predict = "Liu_SKCM", percentile = percentile,
                      KM_output = NULL, survival_var = "OS", iNMDeff_model_type = "continous", binary_probability = NULL)
  # KM_plot_final <- plot_grid(KM_plots[[1]],KM_plots[[2]],KM_plots[[3]],KM_plots[[4]])
  KM_plot_final <- plot_grid(KM_plots[[1]],KM_plots[[2]])
  KM_plot <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/phenotypes/immunotherapy/validation_datasets/Liu_SKCM/Liu_SKCM_NMDeff_OS_survival_curve_percentile_",percentile,".png")
  png(KM_plot, width = 3500, height = 2500, res = 300)
  print(KM_plot_final)
  dev.off()
}

# Binary iNMDeff classification by probabilities

for (prob in c(0.5,0.55,0.6,0.65,0.7)) {
  Immuno_NMDeff_filt <- Immuno_NMDeff_all %>% filter(numPriorTherapies <= 1)
  KM_plots <- survival_curve(Immuno_NMDeff = Immuno_NMDeff_filt, dataset_to_predict = "Liu_SKCM", percentile = NULL,
                      KM_output = NULL, survival_var = "PFS", iNMDeff_model_type = "binary", binary_probability = prob)
  KM_plot_final <- plot_grid(KM_plots[[2]])
  KM_plot <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/phenotypes/immunotherapy/validation_datasets/Liu_SKCM/Liu_SKCM_NMDeff_PFS_survival_curve_probability_",prob,".png")
  png(KM_plot, width = 2500, height = 2500, res = 300)
  print(KM_plot_final)
  dev.off()
  KM_plots <- survival_curve(Immuno_NMDeff = Immuno_NMDeff_filt, dataset_to_predict = "Liu_SKCM", percentile = NULL,
                      KM_output = NULL, survival_var = "OS", iNMDeff_model_type = "binary", binary_probability = prob)
  KM_plot_final <- plot_grid(KM_plots[[2]])
  KM_plot <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/phenotypes/immunotherapy/validation_datasets/Liu_SKCM/Liu_SKCM_NMDeff_OS_survival_curve_probability_",prob,".png")
  png(KM_plot, width = 2500, height = 2500, res = 300)
  print(KM_plot_final)
  dev.off()
}

######## OUTPUT SAVE ########
Immuno_NMDeff_filt <- Immuno_NMDeff_all %>% filter(numPriorTherapies <= 1)
Liu_SKCM_PFS_object <- survival_curve(Immuno_NMDeff = Immuno_NMDeff_filt, dataset_to_predict = "Liu_SKCM", percentile = 20, KM_output = "yes", survival_var = "PFS")
Liu_SKCM_PFS_object_ASE <- survival_curve(Immuno_NMDeff = Immuno_NMDeff_filt, dataset_to_predict = "Liu_SKCM", percentile = 40, KM_output = "yes", survival_var = "PFS")
saveRDS(Liu_SKCM_PFS_object, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig6/Fig6D.RData")

# Cox regression
Liu_cox_reg_res <- c()
binary_probability <- 0.65

for (iNMDeff_model_type in c("continous","binary")) {

  for (NMD_method in c("endogenous_NMD_Consensus","ASE_PTC_NMD_triggering_0.2")) {

    if (NMD_method == "endogenous_NMD_Consensus") {
      percentile <- 20
      NMD_method_char <- "ETG"
    } else if (NMD_method == "ASE_PTC_NMD_triggering_0.2") {
      percentile <- 40
      NMD_method_char <- "ASE"
    }
    survival_var_char <- "PFS"
    status_variable <- "progressed"
    Immuno_NMDeff_filt <- Immuno_NMDeff_all %>% filter(numPriorTherapies <= 1)
    Immuno_NMDeff_filt <- data.frame(Immuno_NMDeff_filt[Immuno_NMDeff_filt$NMD_method == NMD_method,])

    if (iNMDeff_model_type == "continous") {
      NMD_classification <- "NMDeff_mean"
      Immuno_NMDeff_filt$NMD_classification <- NA
      percentiles <- quantile(Immuno_NMDeff_filt[,NMD_classification],probs = seq(0, 1, 0.01), na.rm = TRUE)
      thres_low <- as.numeric(percentiles[paste0(percentile,"%")])
      thres_high <- as.numeric(percentiles[paste0(100-percentile,"%")])
      Immuno_NMDeff_filt[which(Immuno_NMDeff_filt[,NMD_classification] <= thres_low),"NMD_classification"] <- "NMD-Low"
      Immuno_NMDeff_filt[which(Immuno_NMDeff_filt[,NMD_classification] >= thres_high),"NMD_classification"] <- "NMD-High"
    } else if (iNMDeff_model_type == "binary") {
      Immuno_NMDeff_filt$NMD_classification <- case_when(
        Immuno_NMDeff_filt$NMD_High_prob >= binary_probability ~ "NMD-High",
        Immuno_NMDeff_filt$NMD_High_prob <= (1-binary_probability) ~ "NMD-Low",
        TRUE ~ "NMD-Mid"
      ) 
      Immuno_NMDeff_filt <- Immuno_NMDeff_filt %>%
        filter(NMD_classification != "NMD-Mid") %>%
        mutate(NMD_classification = factor(NMD_classification))
      table(Immuno_NMDeff_filt$NMD_classification)
    }

    surv_obj <- with(Immuno_NMDeff_filt, Surv(as.numeric(eval(parse(text = survival_var_char))), as.numeric(eval(parse(text = status_variable)))))
    # km_fit <- survfit(surv_obj ~ NMD_classification, data = Immuno_NMDeff_filt)

    df_res <- data.frame(exp_HR = NA, p_value = NA, CI_2.5 = NA, CI_97.5 = NA)
    Immuno_NMDeff_filt$NMD_classification <- relevel(factor(Immuno_NMDeff_filt$NMD_classification), ref = "NMD-Low")
    Immuno_NMDeff_filt$CNV_burden_quartiles <- as.factor(ntile(Immuno_NMDeff_filt$CNA_prop, 4))
    # coxph(surv_obj ~  NMD_classification + age_quartiles + gender + cancer_type_strat + endogenous_purity + CNV_burden_quartiles + sample_lib_size, data = Immuno_NMDeff_filt)
    coxmodel <- coxph(surv_obj ~  NMD_classification + log(nonsyn_muts) + factor(gender..Male.1..Female.0.) + as.numeric(purity), data = Immuno_NMDeff_filt)
    #coxmodel <- coxph(surv_obj ~  NMD_classification + factor(gender..Male.1..Female.0.) + as.numeric(purity) + CNV_burden_quartiles, data = Immuno_NMDeff_filt)

    coxmodel_res <- summary(coxmodel)
    CI_lower_95 <- coxmodel_res$conf.int[1,"lower .95"]
    CI_upper_95 <- coxmodel_res$conf.int[1,"upper .95"]
    hasTestVarreg <- grep("NMD_classification",rownames(coxmodel_res$coefficients))
    if (length(hasTestVarreg) != 0) {
        df_res[,"exp_HR"] <- coxmodel_res$coefficients[hasTestVarreg,"exp(coef)"]
        df_res[,"p_value"] <- coxmodel_res$coefficients[hasTestVarreg,"Pr(>|z|)"]
        df_res[,"CI_2.5"] <- as.numeric(CI_lower_95)
        df_res[,"CI_97.5"] <- as.numeric(CI_upper_95)
    } 
    df_res$percentile <- percentile
    df_res$NMD_method <- NMD_method_char
    df_res$binary_probability <- binary_probability
    df_res$iNMDeff_model_type <- iNMDeff_model_type
    if (length(Liu_cox_reg_res) == 0) {
      Liu_cox_reg_res <- df_res
    } else {
      Liu_cox_reg_res <- rbind(Liu_cox_reg_res,df_res)
    }
  }
}
Liu_cox_reg_res$cancer_type <- "SKCM"
Liu_cox_reg_res$dataset_name <- "Liu-SKCM"
Liu_cox_reg_res

#################################
#### 4.3) Braun 2020 --> RCC ####
#################################

# Immuno_NMDeff_all <- c()
# for (dataset in c("TCGA","GTEx")) {
#   for (NMD_method in c("endogenous_NMD_Consensus","ASE_PTC_NMD_triggering_0.2")) {
#       Immuno_NMDeff <- NMDeff_gene_model_validation(dataset = dataset, model = "LASSO", 
#                           target_var = NMD_method, dataset_to_predict = "Braun_RCC")
#       Immuno_NMDeff$dataset <- dataset
#       Immuno_NMDeff$NMD_method <- NMD_method
#       if (length(Immuno_NMDeff_all) == 0) {
#         Immuno_NMDeff_all <- Immuno_NMDeff
#       } else {
#         Immuno_NMDeff_all <- rbind(Immuno_NMDeff_all,Immuno_NMDeff)
#       }
#   }
# }
# Immuno_NMDeff_all$treatment_response <- factor(Immuno_NMDeff_all$treatment_response, levels = c("CRPR","CR","PR","PD","SD"))
# out_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/phenotypes/immunotherapy/Immuno_NMDeff_Braun_RCC_results.txt")
# Immuno_NMDeff_all <- read.table(file = out_path, header = TRUE, sep = "\t")
# # write.table(Immuno_NMDeff_all, file = out_path, sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)

# # Plot

# p <- df %>%  
#     # mutate(TMB_Counts = as.numeric(TMB_Counts)) %>%
#     # filter(TMB_Counts >= 80) %>%
#     # filter(Arm == "NIVOLUMAB") %>%
#     filter(treatment_response != "CR") %>%
#     # filter(Tumor_Sample_Primary_or_Metastasis == "PRIMARY") %>%
#     # filter(Number_of_Prior_Therapies == "1") %>%
#     # filter(ImmunoPhenotype != "Infiltrated") %>%
#     # filter(PBRM1 == "WT") %>%
#     filter(WGII <= 0.25) %>%
#     filter(Purity >= 0.5) %>%
#     filter(Ploidy <= 2.5) %>%
#     # filter(Amplification_1q21.3 == "WT") %>%
#     # group_by(NMD_method) %>%
#     # mutate(NMDeff = scale(NMDeff)) %>%
#       ggplot(aes(x = treatment_response, y = NMDeff, fill = factor(treatment_response))) +
#       facet_wrap(dataset + NMD_method ~ ., scale = "free_y") +
#       geom_violin() + scale_fill_brewer(palette = "Dark2") + xlab("") +
#       #coord_cartesian(ylim = c(-3,2)) +
#       geom_boxplot(aes(fill = factor(treatment_response)), width=0.2, color="black", alpha=0.2) +
#       geom_jitter(aes(fill = factor(treatment_response)), 
#                   position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0, dodge.width = 0.9),
#                   alpha = 0.25, size = 1) +        theme_bw(base_size = 25) +
#       theme(axis.text.x = element_text(size = 18, angle = 45, hjust = 1)) #+
#     # annotate("text",
#     #         x = 1:length(table(Immuno_NMDeff_all$treatment_response)),
#     #         y = aggregate( NMDeff ~ treatment_response, Immuno_NMDeff_all, median)[ , 2],
#     #         label = table(Immuno_NMDeff_all$treatment_response),
#     #         col = "black",
#     #         vjust = - 1,
#     #         size = 5) #+
#     # stat_compare_means(comparisons = combinations, size = 12,
#     #                   # label.y = c(2.5,2.25,2.5,2.75,3),
#     #                   label = "p.format", method = "wilcox.test", hide.ns = TRUE)# +
#   output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/phenotypes/immunotherapy/validation_datasets/All_Immuno_Braun_RCC_NMDeff_predictions.png")
#   png(output_path, width = 4500, height = 3000, res = 300)
#   print(p)
#   dev.off()

####################################
#### # 4.4) Carrol 2023 --> EAC ####
####################################

Immuno_NMDeff_all <- c()
for (dataset in c("TCGA","GTEx")) {
  for (NMD_method in c("endogenous_NMD_Consensus","ASE_PTC_NMD_triggering_0.2")) {
      Immuno_NMDeff <- NMDeff_gene_model_validation(dataset = dataset, model = "LASSO", 
                          target_var = NMD_method, dataset_to_predict = "Carrol_EAC")
      Immuno_NMDeff$dataset <- dataset
      Immuno_NMDeff$NMD_method <- NMD_method
      if (length(Immuno_NMDeff_all) == 0) {
        Immuno_NMDeff_all <- Immuno_NMDeff
      } else {
        Immuno_NMDeff_all <- rbind(Immuno_NMDeff_all,Immuno_NMDeff)
      }
  }
}
# Immuno_NMDeff_all$treatment_response <- factor(Immuno_NMDeff_all$treatment_response, levels = c("CR","PR","PD","SD"))
out_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/phenotypes/immunotherapy/Immuno_NMDeff_Carrol_EAC_results.txt")
Immuno_NMDeff_all <- read.table(file = out_path, header = TRUE, sep = "\t")
# write.table(Immuno_NMDeff_all, file = out_path, sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)

# Immuno_NMDeff_all %>%
#       group_by(dataset, NMD_method) %>%
#         summarize(
#           cor_coef = cor.test(ICI.4W.Tumor_INCITE_change, NMDeff)$estimate,
#           p_value = cor.test(ICI.4W.Tumor_INCITE_change, NMDeff)$p.value
#         )
# aggregate( NMDeff ~ Stage, Immuno_NMDeff_all, median)
# Scale
# Immuno_NMDeff_all <- Immuno_NMDeff_all %>%
#                 group_by(dataset,NMD_method) %>%
#                 mutate(NMDeff = scale(NMDeff))

# Mean of two models (GTex and TCGA)
Immuno_NMDeff_all <- Immuno_NMDeff_all %>%
                group_by(sampleID,sample_type,NMD_method) %>%
                mutate(NMDeff_mean = mean(NMDeff, na.rm = TRUE))
# Duplicated (mean)
rows <- which(duplicated(Immuno_NMDeff_all[,c("sampleID","sample_type","NMDeff_mean","NMD_method")]))
Immuno_NMDeff_all <- Immuno_NMDeff_all[-rows,]

# Immuno_NMDeff_all$treatment_response <- ifelse(Immuno_NMDeff_all$treatment_response %in% "CB", "Responders","Non-responders")
# Immuno_NMDeff_all$Best_response <- ifelse(Immuno_NMDeff_all$Best_response %in% c("Progressive Disease (irPD)","Stable Disease (irSD)"), "Non-responders","Responders")
# Immuno_NMDeff_all$treatment_response <- Immuno_NMDeff_all$Best_response
# table(Immuno_NMDeff_all$Best_response,Immuno_NMDeff_all$treatment_response)
# Response_binary <--- this!

# A) Add iNMDeff binary predictions by Random Forest (multiclassPairs model)

# Test with different models (TCGA or GTex)
proxy_model_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/proxy_iNMDeff_models/tissues"
path <- paste0(proxy_model_path,"/ESCA_iNMDeff_RF_classifier.RData")
# path <- paste0(proxy_model_path,"/Esophagus_Mucosa_iNMDeff_RF_classifier.RData")
# path <- paste0(proxy_model_path,"/Esophagus_Muscularis_iNMDeff_RF_classifier.RData")
# path <- paste0(proxy_model_path,"/Esophagus_Gastroesophageal_Junction_iNMDeff_RF_classifier.RData")

iNMDeff_RF_classifier <- readRDS(path)

# Apply on validation datasets (gene expression matrix)
Immuno_dataset <- Carrol_EAC
results <- predict_RF(classifier = iNMDeff_RF_classifier, 
                      Data = Immuno_dataset,
                      impute = TRUE) # can handle missed genes by imputation
results$predictions
iNMDeff_class_res <- data.frame(results$predictions)
colnames(iNMDeff_class_res) <- c("NMD_High_prob","NMD_Low_prob")
#Just to check which is the best model
iNMDeff_class_res$NMD_classification <- case_when(
  iNMDeff_class_res$NMD_High_prob >= 0.6 ~ "NMD-High",
  iNMDeff_class_res$NMD_High_prob <= 0.4 ~ "NMD-Low",
  TRUE ~ "NMD-Mid"
) 
table(iNMDeff_class_res$NMD_classification)

rownames(iNMDeff_class_res) <- gsub("\\.","-",rownames(iNMDeff_class_res))
iNMDeff_class_res$sample_type <- substr(rownames(iNMDeff_class_res),10,50)
iNMDeff_class_res$sampleID <- substr(rownames(iNMDeff_class_res),1,8)
Immuno_NMDeff_all <- merge(Immuno_NMDeff_all,iNMDeff_class_res, 
                by.x = c("sample_type","sampleID"), by.y = c("sample_type","sampleID"),
                all.x = TRUE)

# B) NMDeff vs Treatment Response

# Immuno_NMDeff_all$treatment_response <- ifelse(Immuno_NMDeff_all$treatment_response %in% "CB", "Responders","Non-responders")
# Immuno_NMDeff_all$Best_response <- ifelse(Immuno_NMDeff_all$Best_response %in% c("Progressive Disease (irPD)","Stable Disease (irSD)"), "Non-responders","Responders")
Immuno_NMDeff_all$treatment_response <- Immuno_NMDeff_all$Response_binary
table(Immuno_NMDeff_all$treatment_response)

Carrol_EAC_imm_resp_model <- Immuno_NMDeff_all

# 2) Plots

for (sample_type in c("PostTx_Tumor","PreTx_Tumor","ICI-4W_Tumor")) {

  Immuno_NMDeff_all_filt <- Immuno_NMDeff_all[grep(sample_type,Immuno_NMDeff_all$sample_type),]
  combinations <- combn(names(table(Immuno_NMDeff_all_filt$Response_binary)), 2, simplify = FALSE)

  p <- Immuno_NMDeff_all_filt %>%  
        # filter(total_muts <= 1500 & CNA_prop <= 0.55 & purity >= 0.3 & ploidy <= 3 & Primary_Type != "occult") %>%
        filter(!is.na(Response_binary)) %>%
          ggplot(aes(x = Response_binary, y = NMDeff_mean, fill = factor(Response_binary))) +
            labs(x = "", y = "NMDeff", fill = "") +
            ggtitle("Carrol 2023 - EAC") +
            theme_classic(base_size = 25)+
            # facet_wrap(dataset + NMD_method ~ ., scale = "free_y") +
            facet_wrap(NMD_method ~ ., scale = "free_y") +
            geom_violin() + scale_fill_brewer(palette = "Dark2") + xlab("") +
            geom_boxplot(aes(fill = factor(Response_binary)), width=0.2, color="black", alpha=0.2) +
            geom_jitter(aes(fill = factor(Response_binary)), 
                        position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0, dodge.width = 0.9),
                        alpha = 0.25, size = 1) +        theme_bw(base_size = 25) +
            theme(
                axis.text.x = element_text(size = 14,hjust = 1, angle = 45), 
                legend.position = "top",
                plot.title = element_text(size = 35, hjust = 0.5)) +
            stat_compare_means(comparisons = combinations, size = 8,
                              #label.y = c(1,1.5,2),
                              label = "p.format", method = "wilcox.test", hide.ns = TRUE)# +

  output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/phenotypes/immunotherapy/validation_datasets/Carrol_EAC/Carrol_EAC_NMDeff_vs_treatment_response_boxplots_",sample_type,".png")
  png(output_path, width = 4500, height = 4500, res = 300)
  print(p)
  dev.off()

}

# C) Scatterplot NMDeff vs PFS

for (sample_type in c("PostTx_Tumor","PreTx_Tumor","ICI-4W_Tumor")) {

  Immuno_NMDeff_all_filt <- Immuno_NMDeff_all[grep(sample_type,Immuno_NMDeff_all$sample_type),]
  p <- Immuno_NMDeff_all_filt %>%  
        ggplot(aes(x = PFS, y = NMDeff_mean)) +
          labs(x = "PFS", y = "NMDeff") +
          ggtitle("Carrol 2023 - EAC") +
          theme_classic(base_size = 25)+
          # facet_wrap(dataset + NMD_method ~ ., scale = "free_y") +
          facet_wrap(NMD_method ~ ., scale = "free_y") +
            geom_point()+
            geom_smooth(method = "lm", se = TRUE, size = 1) +
            stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, size = 5) +
          theme(
              axis.text.x = element_text(size = 18,hjust = 1), 
              legend.position = "top",
              plot.title = element_text(size = 35, hjust = 0.5))
  # output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/phenotypes/immunotherapy/validation_datasets/Carrol_EAC/Carrol_EAC_NMDeff_vs_PFS_scatterplot.png")
  output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/phenotypes/immunotherapy/validation_datasets/Carrol_EAC/Carrol_EAC_NMDeff_vs_PFS_scatterplot_",sample_type,".png")
  png(output_path, width = 4500, height = 3500, res = 300)
  print(p)
  dev.off()
}

# D) Survival Curve

table(Immuno_NMDeff_all$sample_type)
# best percentile --> 
for (sample_type in c("PostTx_Tumor","PreTx_Tumor","ICI-4W_Tumor")) {
  for (percentile in c(5,10,15,20,25,30,35,40,45,50)) {
    Immuno_NMDeff_all_filt <- Immuno_NMDeff_all[grep(sample_type,Immuno_NMDeff_all$sample_type),]
    KM_plots <- survival_curve(Immuno_NMDeff = Immuno_NMDeff_all_filt, dataset_to_predict = "Carrol_EAC", percentile = percentile, 
                      KM_output = NULL, survival_var = "PFS", iNMDeff_model_type = "continous", binary_probability = NULL)
    # KM_plot_final <- plot_grid(KM_plots[[1]],KM_plots[[2]],KM_plots[[3]],KM_plots[[4]])
    KM_plot_final <- plot_grid(KM_plots[[1]],KM_plots[[2]])
    KM_plot <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/phenotypes/immunotherapy/validation_datasets/Carrol_EAC/Carrol_EAC_NMDeff_PFS_survival_curve_",sample_type,"_percentile_",percentile,".png")
    png(KM_plot, width = 3500, height = 2500, res = 300)
    print(KM_plot_final)
    dev.off()
    KM_plots <- survival_curve(Immuno_NMDeff = Immuno_NMDeff_all_filt, dataset_to_predict = "Carrol_EAC", percentile = percentile, 
                        KM_output = NULL, survival_var = "OS", iNMDeff_model_type = "continous", binary_probability = NULL)
    # KM_plot_final <- plot_grid(KM_plots[[1]],KM_plots[[2]],KM_plots[[3]],KM_plots[[4]])
    KM_plot_final <- plot_grid(KM_plots[[1]],KM_plots[[2]])
    KM_plot <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/phenotypes/immunotherapy/validation_datasets/Carrol_EAC/Carrol_EAC_NMDeff_OS_survival_curve_",sample_type,"_percentile_",percentile,".png")
    png(KM_plot, width = 3500, height = 2500, res = 300)
    print(KM_plot_final)
    dev.off()
  }
}

for (sample_type in c("PostTx_Tumor","PreTx_Tumor","ICI-4W_Tumor")) {
  for (prob in c(0.5,0.55,0.6,0.65,0.7)) {
      Immuno_NMDeff_all_filt <- Immuno_NMDeff_all[grep(sample_type,Immuno_NMDeff_all$sample_type),]
      KM_plots <- survival_curve(Immuno_NMDeff = Immuno_NMDeff_all_filt, dataset_to_predict = "Carrol_EAC", percentile = NULL, 
                        KM_output = NULL, survival_var = "PFS", iNMDeff_model_type = "binary", binary_probability = prob)
      if (length(KM_plots) != 0) {
        KM_plot_final <- plot_grid(KM_plots[[2]])
        KM_plot <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/phenotypes/immunotherapy/validation_datasets/Carrol_EAC/Carrol_EAC_NMDeff_PFS_survival_curve_",sample_type,"_probability_",prob,".png")
        png(KM_plot, width = 3500, height = 2500, res = 300)
        print(KM_plot_final)
        dev.off()
      }
      if (length(KM_plots) != 0) {
        KM_plots <- survival_curve(Immuno_NMDeff = Immuno_NMDeff_all_filt, dataset_to_predict = "Carrol_EAC", percentile = NULL, 
                            KM_output = NULL, survival_var = "OS", iNMDeff_model_type = "binary", binary_probability = prob)
        KM_plot_final <- plot_grid(KM_plots[[2]])
        KM_plot <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/phenotypes/immunotherapy/validation_datasets/Carrol_EAC/Carrol_EAC_NMDeff_OS_survival_curve_",sample_type,"_probability_",prob,".png")
        png(KM_plot, width = 3500, height = 2500, res = 300)
        print(KM_plot_final)
        dev.off()
      }
  }
}

######## OUTPUT SAVE ########
Immuno_NMDeff_all_filt <- Immuno_NMDeff_all[grep("PreTx_Tumor",Immuno_NMDeff_all$sample_type),]
# N of samples with iNMDeff
sum(!is.na(Immuno_NMDeff_all_filt[Immuno_NMDeff_all_filt$NMD_method == "endogenous_NMD_Consensus","NMDeff"]))
sum(!is.na(Immuno_NMDeff_all_filt[Immuno_NMDeff_all_filt$NMD_method == "ASE_PTC_NMD_triggering_0.2","NMDeff"]))
unique(Immuno_NMDeff_all_filt$sampleID)

Carrol_EAC_PFS_object <- survival_curve(Immuno_NMDeff = Immuno_NMDeff_all_filt, dataset_to_predict = "Carrol_EAC", percentile = 25, KM_output = "yes", survival_var = "PFS")
Carrol_EAC_PFS_object_ASE <- Carrol_EAC_PFS_object
saveRDS(Carrol_EAC_PFS_object, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig6/Fig6E.RData")

# Cox regression
Carrol_cox_reg_res <- c()
binary_probability <- 0.6
for (iNMDeff_model_type in c("binary","continous")) {

  for (NMD_method in c("endogenous_NMD_Consensus","ASE_PTC_NMD_triggering_0.2")) {

    if (NMD_method == "endogenous_NMD_Consensus") {
      percentile <- 25
      NMD_method_char <- "ETG"
    } else if (NMD_method == "ASE_PTC_NMD_triggering_0.2") {
      percentile <- 25
      NMD_method_char <- "ASE"
    }
    survival_var_char <- "PFS"
    status_variable <- "Progressed"

    Immuno_NMDeff_filt <- Immuno_NMDeff_all[grep("PreTx_Tumor",Immuno_NMDeff_all$sample_type),]
    Immuno_NMDeff_filt <- data.frame(Immuno_NMDeff_filt[Immuno_NMDeff_filt$NMD_method == NMD_method,])

    if (iNMDeff_model_type == "continous") {
      NMD_classification <- "NMDeff_mean"
      Immuno_NMDeff_filt$NMD_classification <- NA
      percentiles <- quantile(Immuno_NMDeff_filt[,NMD_classification],probs = seq(0, 1, 0.01), na.rm = TRUE)
      thres_low <- as.numeric(percentiles[paste0(percentile,"%")])
      thres_high <- as.numeric(percentiles[paste0(100-percentile,"%")])
      Immuno_NMDeff_filt[which(Immuno_NMDeff_filt[,NMD_classification] <= thres_low),"NMD_classification"] <- "NMD-Low"
      Immuno_NMDeff_filt[which(Immuno_NMDeff_filt[,NMD_classification] >= thres_high),"NMD_classification"] <- "NMD-High"
    } else if (iNMDeff_model_type == "binary") {
      Immuno_NMDeff_filt$NMD_classification <- case_when(
        Immuno_NMDeff_filt$NMD_High_prob >= binary_probability ~ "NMD-High",
        Immuno_NMDeff_filt$NMD_High_prob <= (1-binary_probability) ~ "NMD-Low",
        TRUE ~ "NMD-Mid"
      ) 
      Immuno_NMDeff_filt <- Immuno_NMDeff_filt %>%
        filter(NMD_classification != "NMD-Mid") %>%
        mutate(NMD_classification = factor(NMD_classification))
      table(Immuno_NMDeff_filt$NMD_classification)
    }

    surv_obj <- with(Immuno_NMDeff_filt, Surv(as.numeric(eval(parse(text = survival_var_char))), as.numeric(eval(parse(text = status_variable)))))
    km_fit <- survfit(surv_obj ~ NMD_classification, data = Immuno_NMDeff_filt)
    # km_fit

    df_res <- data.frame(exp_HR = NA, p_value = NA, CI_2.5 = NA, CI_97.5 = NA)
    Immuno_NMDeff_filt$NMD_classification <- relevel(factor(Immuno_NMDeff_filt$NMD_classification), ref = "NMD-Low")
    # coxph(surv_obj ~  NMD_classification + age_quartiles + gender + cancer_type_strat + endogenous_purity + CNV_burden_quartiles + sample_lib_size, data = Immuno_NMDeff_filt)
    coxmodel <- coxph(surv_obj ~ NMD_classification + as.numeric(Age) + factor(Sex) + log(PreTx.Tumor_Mutational_Load), data = Immuno_NMDeff_filt)
    #coxmodel <- coxph(surv_obj ~ NMD_classification + as.numeric(Age) + factor(Sex), data = Immuno_NMDeff_filt)
    # surv_obj <- with(Immuno_NMDeff_filt, Surv(as.numeric(eval(parse(text = survival_var_char))), as.numeric(eval(parse(text = status_variable)))))
    # coxmodel <- coxph(surv_obj ~  log(PreTx.Tumor_Mutational_Load) + as.numeric(Age) + factor(Sex), data = Immuno_NMDeff_filt)
    # coxmodel <- coxph(surv_obj ~  log(exp(NMDeff)) + as.numeric(Age) + factor(Sex), data = Immuno_NMDeff_filt)

    coxmodel_res <- summary(coxmodel)
    CI_lower_95 <- coxmodel_res$conf.int[1,"lower .95"]
    CI_upper_95 <- coxmodel_res$conf.int[1,"upper .95"]
    hasTestVarreg <- grep("NMD_classification",rownames(coxmodel_res$coefficients))
    if (length(hasTestVarreg) != 0) {
        df_res[,"exp_HR"] <- coxmodel_res$coefficients[hasTestVarreg,"exp(coef)"]
        df_res[,"p_value"] <- coxmodel_res$coefficients[hasTestVarreg,"Pr(>|z|)"]
        df_res[,"CI_2.5"] <- as.numeric(CI_lower_95)
        df_res[,"CI_97.5"] <- as.numeric(CI_upper_95)
    } 
    df_res$percentile <- percentile
    df_res$NMD_method <- NMD_method_char
    df_res$binary_probability <- binary_probability
    df_res$iNMDeff_model_type <- iNMDeff_model_type
    if (length(Carrol_cox_reg_res) == 0) {
      Carrol_cox_reg_res <- df_res
    } else {
      Carrol_cox_reg_res <- rbind(Carrol_cox_reg_res,df_res)
    }
  }
}
Carrol_cox_reg_res$cancer_type <- "EAC"
Carrol_cox_reg_res$dataset_name <- "Carrol-EAC"
Carrol_cox_reg_res

####################################
##### 4.5) Motzer 2020 --> RCC #####
####################################

# This data does not have treatment response, I personally asked the authors:
"PFS is in months. PD1FL is a flag for whether the sample met the criteria for PD-L1 IHC positivity. 
Unfortunately, the RECIST response data was not released with the manuscript, you can use the PFS and censoring data for survival analysis."

Immuno_NMDeff_all <- c()
for (dataset in c("TCGA","GTEx")) {
  for (NMD_method in c("endogenous_NMD_Consensus","ASE_PTC_NMD_triggering_0.2")) {
      Immuno_NMDeff <- NMDeff_gene_model_validation(dataset = dataset, model = "LASSO", 
                          target_var = NMD_method, dataset_to_predict = "Motzer_RCC")
      Immuno_NMDeff$dataset <- dataset
      Immuno_NMDeff$NMD_method <- NMD_method
      if (length(Immuno_NMDeff_all) == 0) {
        Immuno_NMDeff_all <- Immuno_NMDeff
      } else {
        Immuno_NMDeff_all <- rbind(Immuno_NMDeff_all,Immuno_NMDeff)
      }
  }
}
# Immuno_NMDeff_all$treatment_response <- factor(Immuno_NMDeff_all$treatment_response, levels = c("CR","PR","PD","SD"))
out_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/phenotypes/immunotherapy/Immuno_NMDeff_Motzer_RCC_results.txt")
Immuno_NMDeff_all <- read.table(file = out_path, header = TRUE, sep = "\t")
# write.table(Immuno_NMDeff_all, file = out_path, sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
colnames(Immuno_NMDeff_all)[colnames(Immuno_NMDeff_all) %in% "treatment_response"] <- "PFS_P_CNSR"

# Mean of four models
Immuno_NMDeff_all <- Immuno_NMDeff_all %>%
                group_by(Row.names,NMD_method) %>%
                mutate(NMDeff_mean = mean(NMDeff, na.rm = TRUE))
# Duplicated (mean)
rows <- which(duplicated(Immuno_NMDeff_all[,c("Row.names","NMDeff_mean","NMD_method")]))
Immuno_NMDeff_all <- Immuno_NMDeff_all[-rows,]

# N of samples with iNMDeff
# sum(!is.na(Immuno_NMDeff_all[Immuno_NMDeff_all$NMD_method == "endogenous_NMD_Consensus","NMDeff"]))
# sum(!is.na(Immuno_NMDeff_all[Immuno_NMDeff_all$NMD_method == "ASE_PTC_NMD_triggering_0.2","NMDeff"]))
# unique(Immuno_NMDeff_all$Row.names)

# A) Add iNMDeff binary predictions by Random Forest (multiclassPairs model)

# Test with different models (TCGA or GTex)
proxy_model_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/proxy_iNMDeff_models/tissues"
path <- paste0(proxy_model_path,"/Kidney_Cortex_iNMDeff_RF_classifier.RData")
# path <- paste0(proxy_model_path,"/KIRC_iNMDeff_RF_classifier.RData")
# path <- paste0(proxy_model_path,"/KIRP_iNMDeff_RF_classifier.RData")

iNMDeff_RF_classifier <- readRDS(path)

# Apply on validation datasets (gene expression matrix)
Immuno_dataset <- Motzer_RCC
Immuno_dataset$gene_id <- NULL
results <- predict_RF(classifier = iNMDeff_RF_classifier, 
                      Data = Immuno_dataset,
                      impute = TRUE) # can handle missed genes by imputation
results$predictions
iNMDeff_class_res <- data.frame(results$predictions)
colnames(iNMDeff_class_res) <- c("NMD_High_prob","NMD_Low_prob")
#Just to check which is the best model
iNMDeff_class_res$NMD_classification <- case_when(
  iNMDeff_class_res$NMD_High_prob >= 0.7 ~ "NMD-High",
  iNMDeff_class_res$NMD_High_prob <= 0.3 ~ "NMD-Low",
  TRUE ~ "NMD-Mid"
) 
table(iNMDeff_class_res$NMD_classification)

Immuno_NMDeff_all <- merge(Immuno_NMDeff_all,iNMDeff_class_res, 
                by.x = c("Row.names"), by.y = c("row.names"),
                all.x = TRUE)

# B) Scatterplot NMDeff vs PFS

p <- Immuno_NMDeff_all %>%  
      ggplot(aes(x = PFS_P, y = NMDeff_mean, fill = factor(TRT01P))) +
        labs(x = "PFS", y = "NMDeff", fill = "Treatment") +
        ggtitle("Motzer 2020 - RCC") +
        theme_classic(base_size = 25)+
        # facet_wrap(dataset + NMD_method ~ ., scale = "free_y") +
        facet_wrap(NMD_method ~ ., scale = "free_y") +
          geom_point()+
          geom_smooth(method = "lm", se = TRUE, size = 1) +
          stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, size = 5) +
        theme(
            axis.text.x = element_text(size = 18,hjust = 1), 
            legend.position = "top",
            plot.title = element_text(size = 35, hjust = 0.5))
output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/phenotypes/immunotherapy/validation_datasets/Motzer_RCC/Motzer_RCC_NMDeff_vs_PFS_scatterplot.png")
png(output_path, width = 4500, height = 3500, res = 300)
print(p)
dev.off()

# C) Survival curve
Immuno_NMDeff_all$PFS_P <- as.numeric(Immuno_NMDeff_all$PFS_P)
# best percentile --> 30
for (percentile in c(5,10,15,20,25,30,35,40,45,50)) {
  KM_plots <- survival_curve(Immuno_NMDeff = Immuno_NMDeff_all, dataset_to_predict = "Motzer_RCC", percentile = percentile, 
              KM_output = NULL, survival_var = "PFS", iNMDeff_model_type = "continous", binary_probability = NULL)
  # KM_plot_final <- plot_grid(KM_plots[[1]],KM_plots[[2]],KM_plots[[3]],KM_plots[[4]])
  KM_plot_final <- plot_grid(KM_plots[[1]],KM_plots[[2]])
  KM_plot <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/phenotypes/immunotherapy/validation_datasets/Motzer_RCC/Motzer_RCC_NMDeff_PFS_survival_curve_percentile_",percentile,".png")
  png(KM_plot, width = 3500, height = 2500, res = 300)
  print(KM_plot_final)
  dev.off()
}

for (prob in c(0.5,0.55,0.6,0.65,0.7)) {
    KM_plots <- survival_curve(Immuno_NMDeff = Immuno_NMDeff_all, dataset_to_predict = "Motzer_RCC", percentile = NULL, 
                      KM_output = NULL, survival_var = "PFS", iNMDeff_model_type = "binary", binary_probability = prob)
    if (length(KM_plots) != 0) {
      KM_plot_final <- plot_grid(KM_plots[[2]])
      KM_plot <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/phenotypes/immunotherapy/validation_datasets/Motzer_RCC/Motzer_RCC_NMDeff_PFS_survival_curve_probability_",prob,".png")
      png(KM_plot, width = 3500, height = 2500, res = 300)
      print(KM_plot_final)
      dev.off()
    }
}


######## OUTPUT SAVE ########
Motzer_RCC_PFS_object <- survival_curve(Immuno_NMDeff = Immuno_NMDeff_all, dataset_to_predict = "Motzer_RCC", percentile = 30, KM_output = "yes", survival_var = "PFS")
Motzer_RCC_PFS_object_ASE <- survival_curve(Immuno_NMDeff = Immuno_NMDeff_all, dataset_to_predict = "Motzer_RCC", percentile = 45, KM_output = "yes", survival_var = "PFS")
saveRDS(Motzer_RCC_PFS_object, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig6/Fig6F.RData")

# Cox regression
Motzer_cox_reg_res <- c()
binary_probability <- 0.7

for (iNMDeff_model_type in c("continous","binary")) {

  for (NMD_method in c("endogenous_NMD_Consensus","ASE_PTC_NMD_triggering_0.2")) {

    if (NMD_method == "endogenous_NMD_Consensus") {
      percentile <- 30
      NMD_method_char <- "ETG"
    } else if (NMD_method == "ASE_PTC_NMD_triggering_0.2") {
      percentile <- 45
      NMD_method_char <- "ASE"
    }
    
    survival_var_char <- "PFS_P"
    status_variable <- "PFS_P_CNSR"
    # status_variable <- "treatment_response" # It's not the treatment response...!
    # Immuno_NMDeff_filt <- data.frame(Immuno_NMDeff_all[Immuno_NMDeff_all$NMD_method == NMD_method,])
    Immuno_NMDeff_filt <- data.frame(Immuno_NMDeff_all[Immuno_NMDeff_all$NMD_method == NMD_method & Immuno_NMDeff_all$TRT01P != "Sunitinib",])
    # N of samples with iNMDeff
    sum(!is.na(Immuno_NMDeff_filt[Immuno_NMDeff_filt$NMD_method == "endogenous_NMD_Consensus","NMDeff"]))
    sum(!is.na(Immuno_NMDeff_filt[Immuno_NMDeff_filt$NMD_method == "ASE_PTC_NMD_triggering_0.2","NMDeff"]))
    unique(Immuno_NMDeff_filt$Row.names)

    if (iNMDeff_model_type == "continous") {
      NMD_classification <- "NMDeff_mean"
      Immuno_NMDeff_filt$NMD_classification <- NA
      percentiles <- quantile(Immuno_NMDeff_filt[,NMD_classification],probs = seq(0, 1, 0.01), na.rm = TRUE)
      thres_low <- as.numeric(percentiles[paste0(percentile,"%")])
      thres_high <- as.numeric(percentiles[paste0(100-percentile,"%")])
      Immuno_NMDeff_filt[which(Immuno_NMDeff_filt[,NMD_classification] <= thres_low),"NMD_classification"] <- "NMD-Low"
      Immuno_NMDeff_filt[which(Immuno_NMDeff_filt[,NMD_classification] >= thres_high),"NMD_classification"] <- "NMD-High"
    } else if (iNMDeff_model_type == "binary") {
      Immuno_NMDeff_filt$NMD_classification <- case_when(
        Immuno_NMDeff_filt$NMD_High_prob >= binary_probability ~ "NMD-High",
        Immuno_NMDeff_filt$NMD_High_prob <= (1-binary_probability) ~ "NMD-Low",
        TRUE ~ "NMD-Mid"
      ) 
      Immuno_NMDeff_filt <- Immuno_NMDeff_filt %>%
        filter(NMD_classification != "NMD-Mid") %>%
        mutate(NMD_classification = factor(NMD_classification))
      table(Immuno_NMDeff_filt$NMD_classification)
    }

    surv_obj <- with(Immuno_NMDeff_filt, Surv(as.numeric(eval(parse(text = survival_var_char))), as.numeric(eval(parse(text = status_variable)))))
    km_fit <- survfit(surv_obj ~ NMD_classification, data = Immuno_NMDeff_filt)

    df_res <- data.frame(exp_HR = NA, p_value = NA, CI_2.5 = NA, CI_97.5 = NA)
    Immuno_NMDeff_filt$NMD_classification <- relevel(factor(Immuno_NMDeff_filt$NMD_classification), ref = "NMD-Low")
    # coxph(surv_obj ~  NMD_classification + age_quartiles + gender + cancer_type_strat + endogenous_purity + CNV_burden_quartiles + sample_lib_size, data = Immuno_NMDeff_filt)
    coxmodel <- coxph(surv_obj ~ NMD_classification + as.numeric(AGE) + factor(SEX), data = Immuno_NMDeff_filt)

    coxmodel_res <- summary(coxmodel)
    CI_lower_95 <- coxmodel_res$conf.int[1,"lower .95"]
    CI_upper_95 <- coxmodel_res$conf.int[1,"upper .95"]
    hasTestVarreg <- grep("NMD_classification",rownames(coxmodel_res$coefficients))
    if (length(hasTestVarreg) != 0) {
        df_res[,"exp_HR"] <- coxmodel_res$coefficients[hasTestVarreg,"exp(coef)"]
        df_res[,"p_value"] <- coxmodel_res$coefficients[hasTestVarreg,"Pr(>|z|)"]
        df_res[,"CI_2.5"] <- as.numeric(CI_lower_95)
        df_res[,"CI_97.5"] <- as.numeric(CI_upper_95)
    } 
    df_res$percentile <- percentile
    df_res$NMD_method <- NMD_method_char
    df_res$iNMDeff_model_type <- iNMDeff_model_type
    df_res$binary_probability <- binary_probability
    if (length(Motzer_cox_reg_res) == 0) {
      Motzer_cox_reg_res <- df_res
    } else {
      Motzer_cox_reg_res <- rbind(Motzer_cox_reg_res,df_res)
    }
  }

}
Motzer_cox_reg_res$cancer_type <- "RCC"
Motzer_cox_reg_res$dataset_name <- "Motzer-RCC"
Motzer_cox_reg_res

#######################################
##### 4.6) Riaz 2017 --> Melanoma #####
#######################################

Immuno_NMDeff_all <- c()
for (dataset in c("TCGA","GTEx")) {
  for (NMD_method in c("endogenous_NMD_Consensus","ASE_PTC_NMD_triggering_0.2")) {
      Immuno_NMDeff <- NMDeff_gene_model_validation(dataset = dataset, model = "LASSO", 
                          target_var = NMD_method, dataset_to_predict = "Riaz_SKCM")
      Immuno_NMDeff$dataset <- dataset
      Immuno_NMDeff$NMD_method <- NMD_method
      if (length(Immuno_NMDeff_all) == 0) {
        Immuno_NMDeff_all <- Immuno_NMDeff
      } else {
        Immuno_NMDeff_all <- rbind(Immuno_NMDeff_all,Immuno_NMDeff)
      }
  }
}
# Immuno_NMDeff_all$treatment_response <- factor(Immuno_NMDeff_all$treatment_response, levels = c("CR","PR","PD","SD"))
out_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/phenotypes/immunotherapy/Immuno_NMDeff_Riaz_SKCM_results.txt")
Immuno_NMDeff_all <- read.table(file = out_path, header = TRUE, sep = "\t")
# df <- read.csv(file = out_path, header = TRUE, sep = "\t")
# write.table(Immuno_NMDeff_all, file = out_path, sep = "\t", quote = TRUE, col.names = TRUE, row.names = TRUE)

# Mean of four models
Immuno_NMDeff_all <- Immuno_NMDeff_all %>%
                group_by(patientID,NMD_method) %>%
                mutate(NMDeff_mean = mean(NMDeff, na.rm = TRUE))
# Duplicated (mean)
rows <- which(duplicated(Immuno_NMDeff_all[,c("patientID","NMDeff_mean","NMD_method")]))
Immuno_NMDeff_all <- Immuno_NMDeff_all[-rows,]

# N of samples with iNMDeff
# sum(!is.na(Immuno_NMDeff_filt[Immuno_NMDeff_filt$NMD_method == "endogenous_NMD_Consensus","NMDeff"]))
# sum(!is.na(Immuno_NMDeff_filt[Immuno_NMDeff_filt$NMD_method == "ASE_PTC_NMD_triggering_0.2","NMDeff"]))
# unique(Immuno_NMDeff_filt$Row.names)

# names
colnames(Immuno_NMDeff_all)[colnames(Immuno_NMDeff_all) %in% c("Dead.Alive..Dead...True.","Time.to.Death..weeks.")] <- c("status","time_to_death_weeks")
Immuno_NMDeff_all$treatment_response <- ifelse(Immuno_NMDeff_all$treatment_response %in% c("CR","PR"), "Responders","Non-Responders")
Riaz_imm_resp_model <- Immuno_NMDeff_all

# A) NMDeff vs Treatment Response
for (treatment in c("Pre","On")) {

  Immuno_NMDeff_all_filt <- Immuno_NMDeff_all %>% filter(grepl(treatment,patient_full_ID))
  Immuno_NMDeff_all_filt$treatment_response <- ifelse(Immuno_NMDeff_all_filt$treatment_response %in% c("CR","PR"), "Responders","Non-Responders")
  # Immuno_NMDeff_all$treatment_response <- ifelse(Immuno_NMDeff_all$treatment_response %in% c("CR","PR"), "Responders","Non-responders")
  combinations <- combn(names(table(Immuno_NMDeff_all_filt$treatment_response)), 2, simplify = FALSE)

  p <- Immuno_NMDeff_all_filt %>%  
        # filter(total_muts <= 1500 & CNA_prop <= 0.55 & purity >= 0.3 & ploidy <= 3 & Primary_Type != "occult") %>%
        filter(!is.na(treatment_response)) %>%
        filter(treatment_response != "NE") %>%
          ggplot(aes(x = treatment_response, y = NMDeff_mean, fill = factor(treatment_response))) +
            labs(x = "", y = "NMDeff", fill = "") +
            ggtitle("Riaz 2017 - Melanoma") +
            theme_classic(base_size = 25)+
            # facet_wrap(dataset + NMD_method ~ ., scale = "free_y") +
            facet_wrap(NMD_method ~ ., scale = "free_y") +
            geom_violin() + scale_fill_brewer(palette = "Dark2") + xlab("") +
            geom_boxplot(aes(fill = factor(treatment_response)), width=0.2, color="black", alpha=0.2) +
            geom_jitter(aes(fill = factor(treatment_response)), 
                        position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0, dodge.width = 0.9),
                        alpha = 0.25, size = 1) +        theme_bw(base_size = 25) +
            theme(
                axis.text.x = element_text(size = 18,hjust = 1), 
                legend.position = "top",
                plot.title = element_text(size = 35, hjust = 0.5)) +
            stat_compare_means(comparisons = combinations, size = 6,
                              #label.y = c(1,1.5,2),
                              label = "p.format", method = "wilcox.test", hide.ns = TRUE)
  output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/phenotypes/immunotherapy/validation_datasets/Riaz_SKCM/Riaz_SKCM_NMDeff_vs_treatment_response_boxplots_",treatment,"_treatment.png")
  png(output_path, width = 4500, height = 3500, res = 300)
  print(p)
  dev.off()

}

# B) Scatterplot NMDeff vs PFS
for (treatment in c("Pre","On")) {
  Immuno_NMDeff_all_filt <- Immuno_NMDeff_all %>% filter(grepl(treatment,patient_full_ID))
  p <- Immuno_NMDeff_all_filt %>%  
        ggplot(aes(x = time_to_death_weeks, y = NMDeff_mean)) +
          labs(x = "time_to_death_weeks", y = "NMDeff") +
          ggtitle("Riaz 2017 - Melanoma") +
          theme_classic(base_size = 25)+
          facet_wrap(dataset + NMD_method ~ ., scale = "free_y") +
            geom_point()+
            geom_smooth(method = "lm", se = TRUE, size = 1) +
            stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, size = 5) +
          theme(
              axis.text.x = element_text(size = 18,hjust = 1), 
              legend.position = "top",
              plot.title = element_text(size = 35, hjust = 0.5))
  # output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/phenotypes/immunotherapy/validation_datasets/Riaz_SKCM/Riaz_SKCM_NMDeff_vs_PFS_scatterplot.png")
  output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/phenotypes/immunotherapy/validation_datasets/Riaz_SKCM/Riaz_SKCM_NMDeff_vs_PFS_scatterplot_",treatment,"_treatment.png")
  png(output_path, width = 4500, height = 3500, res = 300)
  print(p)
  dev.off()
}

# C) Survival Curves

for (treatment in c("Pre","On")) {
  for (percentile in c(5,10,15,20,25,30,35,40,45,50)) {
    Immuno_NMDeff_all_filt <- Immuno_NMDeff_all %>% filter(grepl(treatment,patient_full_ID))
    # KM_plots <- survival_curve(Immuno_NMDeff = Immuno_NMDeff_all_filt, dataset_to_predict = "Riaz_SKCM", percentile = percentile,KM_output = NULL, survival_var = "PFS")
    # KM_plot_final <- plot_grid(KM_plots[[1]],KM_plots[[2]])
    # KM_plot <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/phenotypes/immunotherapy/validation_datasets/Riaz_SKCM/Riaz_SKCM_NMDeff_PFS_survival_curve_",treatment,"_treatment_percentile_",percentile,".png")
    # png(KM_plot, width = 3500, height = 2500, res = 300)
    # print(KM_plot_final)
    # dev.off()
    KM_plots <- survival_curve(Immuno_NMDeff = Immuno_NMDeff_all_filt, dataset_to_predict = "Riaz_SKCM", percentile = percentile, KM_output = NULL, survival_var = "OS")
    KM_plot_final <- plot_grid(KM_plots[[1]],KM_plots[[2]])
    KM_plot <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/phenotypes/immunotherapy/validation_datasets/Riaz_SKCM/Riaz_SKCM_NMDeff_OS_survival_curve_",treatment,"_treatment_percentile_",percentile,".png")
    png(KM_plot, width = 3500, height = 2500, res = 300)
    print(KM_plot_final)
    dev.off()
  }
}

# N of samples with iNMDeff
Immuno_NMDeff_all_filt <- Immuno_NMDeff_all %>% filter(grepl("On",patient_full_ID))
sum(!is.na(Immuno_NMDeff_all_filt[Immuno_NMDeff_all_filt$NMD_method == "endogenous_NMD_Consensus","NMDeff"]))
sum(!is.na(Immuno_NMDeff_all_filt[Immuno_NMDeff_all_filt$NMD_method == "ASE_PTC_NMD_triggering_0.2","NMDeff"]))
unique(Immuno_NMDeff_all_filt$patientID)

######## OUTPUT SAVE ########
# Riaz_SKCM_OS_object <- survival_curve(Immuno_NMDeff = Immuno_NMDeff_all_filt, dataset_to_predict = "Riaz_SKCM", percentile = percentile, KM_output = "yes", survival_var = "OS")

#######################################
#### # 4.7) Hartwig --> Pan-cancer ####
#######################################

Immuno_NMDeff_all <- c()
for (dataset in c("TCGA","GTEx")) {
  print(dataset)
  for (NMD_method in c("endogenous_NMD_Consensus","ASE_PTC_NMD_triggering_0.2")) {
    print(NMD_method)
      Immuno_NMDeff <- NMDeff_gene_model_validation(dataset = dataset, model = "LASSO", 
                          target_var = NMD_method, dataset_to_predict = "Hartwig")
      Immuno_NMDeff$dataset <- dataset
      Immuno_NMDeff$NMD_method <- NMD_method
      if (length(Immuno_NMDeff_all) == 0) {
        Immuno_NMDeff_all <- Immuno_NMDeff
      } else {
        Immuno_NMDeff_all <- rbind(Immuno_NMDeff_all,Immuno_NMDeff)
      }
  }
}
# Immuno_NMDeff_all$treatment_response <- factor(Immuno_NMDeff_all$treatment_response, levels = c("CR","PR","PD","SD"))
out_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/phenotypes/immunotherapy/Immuno_NMDeff_Hartwig_results.txt")
Immuno_NMDeff_all <- read.table(file = out_path, header = TRUE, sep = "\t")
# write.table(Immuno_NMDeff_all, file = out_path, sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
Hartwig_metadata_2_filt <- Hartwig_metadata_2[,!colnames(Hartwig_metadata_2) %in% c("gender","primaryTumorLocation","biopsySite","biopsyLocation")]

# Add TMB
table(Immuno_NMDeff_all$hmfPatientId %in% Hartwig_metadata_2_filt$patientId)
Immuno_NMDeff_all <- merge(Immuno_NMDeff_all, Hartwig_metadata_2_filt,
                    by.x = "hmfPatientId", by.y = "patientId", all.x = TRUE)
#Some Individuals have more than 1 sample...
#Remove
samples_remove <- names(table(Immuno_NMDeff_all$hmfPatientId)[table(Immuno_NMDeff_all$hmfPatientId) != 4])
Immuno_NMDeff_all <- Immuno_NMDeff_all[-which(Immuno_NMDeff_all$hmfPatientId %in% samples_remove),]
# Mean of four models
Immuno_NMDeff_all <- Immuno_NMDeff_all %>%
                group_by(hmfPatientId,NMD_method) %>%
                mutate(NMDeff_mean = mean(NMDeff, na.rm = TRUE))

# Immuno_NMDeff_all$NMDeff_mean <- Immuno_NMDeff_all$NMDeff
# Immuno_NMDeff_all <- Immuno_NMDeff_all[Immuno_NMDeff_all$dataset == "TCGA",]
# Duplicated (mean)
rows <- which(duplicated(Immuno_NMDeff_all[,c("hmfPatientId","NMDeff_mean","NMD_method")]))
Immuno_NMDeff_all <- data.frame(Immuno_NMDeff_all[-rows,])

# Filters
#treatment
Immuno_NMDeff_all_filt <- Immuno_NMDeff_all[Immuno_NMDeff_all$consolidatedTreatmentType == "Immunotherapy",]
Immuno_NMDeff_all_filt$status <- ifelse(Immuno_NMDeff_all_filt$deathDate == "null", 0, 1)

# A) NMDeff vs Treatment Response

Immuno_NMDeff_all_filt$treatment_response <- ifelse(Immuno_NMDeff_all_filt$treatment_response %in% c("ND","MR"),NA,as.character(Immuno_NMDeff_all_filt$treatment_response))
Immuno_NMDeff_all_filt$treatment_response <- ifelse(Immuno_NMDeff_all_filt$treatment_response %in% c("CR","PR"), "Responders",as.character(Immuno_NMDeff_all_filt$treatment_response))
Immuno_NMDeff_all_filt$treatment_response <- ifelse(Immuno_NMDeff_all_filt$treatment_response %in% c("PD","Non-CR/Non-PD","Clinical progression","SD"), "Non-Responders",as.character(Immuno_NMDeff_all_filt$treatment_response))
Immuno_NMDeff_all_filt <- Immuno_NMDeff_all_filt[!is.na(Immuno_NMDeff_all_filt$treatment_response),]
Immuno_NMDeff_all_filt$cancer <- paste0(Immuno_NMDeff_all_filt$cancerType,"-",Immuno_NMDeff_all_filt$primaryTumorType)
Immuno_NMDeff_all_filt$cancer <- ifelse(Immuno_NMDeff_all_filt$cancer %in% c("NA-Melanoma","Skin-Melanoma","Other-Melanoma"), "Melanoma", Immuno_NMDeff_all_filt$cancer)
Immuno_NMDeff_all_filt$cancer <- gsub("NA-","",Immuno_NMDeff_all_filt$cancer)
Immuno_NMDeff_all_filt$cancer <- gsub("-","",Immuno_NMDeff_all_filt$cancer)
Immuno_NMDeff_all_filt$cancer <- factor(Immuno_NMDeff_all_filt$cancer)

# Immuno_NMDeff_all_filt$cancer <- factor(Immuno_NMDeff_all_filt$primaryTumorType)

table(Immuno_NMDeff_all_filt$cancer)
table(Immuno_NMDeff_all_filt$treatment_response)

Hartwig_imm_resp_model <- Immuno_NMDeff_all_filt

# N of samples with iNMDeff
sum(!is.na(Hartwig_imm_resp_model[Hartwig_imm_resp_model$NMD_method == "endogenous_NMD_Consensus","NMDeff"]))
sum(!is.na(Hartwig_imm_resp_model[Hartwig_imm_resp_model$NMD_method == "ASE_PTC_NMD_triggering_0.2","NMDeff"]))
unique(Hartwig_imm_resp_model$Row.names)

###############################################################################
############################## FINAL OUTPUT SAVE ##############################
###############################################################################

datasets_pvalues <- c(0.16,0.27,0.17)
final_meta_pvalue <- survcomp::combine.test(datasets_pvalues, method = c("fisher"), hetero = FALSE, na.rm = FALSE)
datasets_pvalues <- c(0.55,0.7,0.44)
final_meta_pvalue <- survcomp::combine.test(datasets_pvalues, method = c("fisher"), hetero = FALSE, na.rm = FALSE)
final_meta_pvalue

## Supplementary Fig 21 ##

surv_objects <- list("Motzer_RCC_PFS" = Motzer_RCC_PFS_object, "Carrol_EAC_PFS" = Carrol_EAC_PFS_object, "Liu_SKCM_PFS" = Liu_SKCM_PFS_object)
saveRDS(surv_objects, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig6/Fig6D_F.RData")
surv_objects_ASE <- list("Motzer_RCC_PFS" = Motzer_RCC_PFS_object_ASE, 
                        "Carrol_EAC_PFS" = Carrol_EAC_PFS_object_ASE, 
                        "Liu_SKCM_PFS" = Liu_SKCM_PFS_object_ASE)
saveRDS(surv_objects_ASE, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig21/SuppFig21D_F.RData")

## Supplementary Table S7 ##
# conda activate /home/gpalou/.conda/envs/survival_analysis
library(survcomp)

# Liu --> 0.2233831
# survcomp::combine.test(c(0.2311370,0.2512966), method = c("fisher"), hetero = FALSE, na.rm = FALSE)
Liu_cox_reg_res$meta_p_value <- NA
Liu_cox_reg_res[Liu_cox_reg_res$NMD_method == "ETG","meta_p_value"] <- 0.2233831
# Carrol --> 0.2727334
# survcomp::combine.test(c(0.5316254,0.1436005), method = c("fisher"), hetero = FALSE, na.rm = FALSE)
Carrol_cox_reg_res$meta_p_value <- NA
Carrol_cox_reg_res[Carrol_cox_reg_res$NMD_method == "ETG","meta_p_value"] <- 0.2727334
# Motzer --> 0.1308121
# survcomp::combine.test(c(0.1678148,0.1713617), method = c("fisher"), hetero = FALSE, na.rm = FALSE)
Motzer_cox_reg_res$meta_p_value <- NA
Motzer_cox_reg_res[Motzer_cox_reg_res$NMD_method == "ETG","meta_p_value"] <- 0.1308121

# Merge
PFS_cox_reg_all <- rbind(rbind(Liu_cox_reg_res,Motzer_cox_reg_res),Carrol_cox_reg_res)
PFS_cox_reg_all$outcome <- "PFS"
PFS_cox_reg_all$treatment <- "immunotherapy"
PFS_cox_reg_all_SuppTableS7 <- PFS_cox_reg_all %>%
    dplyr::select(dataset_name,cancer_type,percentile,NMD_method,outcome,treatment,exp_HR,CI_2.5,CI_97.5,p_value,meta_p_value) %>%
    arrange(cancer_type,NMD_method)
colnames(PFS_cox_reg_all_SuppTableS7) <- c("dataset_name","cancer_type","percentile","NMD_method","outcome","treatment","exp_HR",
                              "CI_low","CI_high","p_value","meta_p_value")
PFS_cox_reg_all_SuppTableS7
output_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/Tables/SuppTableS7.txt"
write.table(PFS_cox_reg_all_SuppTableS7, file = output_path, 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)   

############################## IMMUNOTHERAPY PREDICTIONS ##############################


### 6) Immunotherapy response prediction by logistic regression
# library(pROC)

# Joint model using 4 datasets: TCGA-SKCM, Liu, Carrol, Motzer (no data), hartwig
WES_size <- 36 # 36Mb
# TCGA
TCGA_SKCM_imm_resp_model_df <- TCGA_SKCM_imm_resp_model[,c("NMDeff","NMD_method","treatment_response","TMB","age_at_diagnosis","gender")]
colnames(TCGA_SKCM_imm_resp_model_df)[colnames(TCGA_SKCM_imm_resp_model_df) 
                  %in% c("NMDeff","age_at_diagnosis","gender")] <- c("NMDeff_mean","age","sex")
table(TCGA_SKCM_imm_resp_model_df$treatment_response)
TCGA_SKCM_imm_resp_model_df$cancer_type <- "skin"
# Liu
Liu_SKCM_imm_resp_model_df <- data.frame(Liu_SKCM_imm_resp_model[,c("NMDeff_mean","NMD_method","treatment_response","total_muts","Primary_Type","gender..Male.1..Female.0.")])
colnames(Liu_SKCM_imm_resp_model_df)[colnames(Liu_SKCM_imm_resp_model_df) 
                  %in% c("total_muts","Primary_Type","gender..Male.1..Female.0.")] <- c("TMB","cancer_type","sex")
Liu_SKCM_imm_resp_model_df$treatment_response <- ifelse(Liu_SKCM_imm_resp_model_df$treatment_response == "Non-responders","Non-Responders","Responders")
Liu_SKCM_imm_resp_model_df$TMB <- Liu_SKCM_imm_resp_model_df$TMB/WES_size
table(Liu_SKCM_imm_resp_model_df$treatment_response)
# Carrol
Carrol_EAC_imm_resp_model_df <- data.frame(Carrol_EAC_imm_resp_model[,c("NMDeff_mean","sample_type","NMD_method","treatment_response","PreTx.Tumor_Mutational_Load","Age","Sex")])
Carrol_EAC_imm_resp_model_df <- Carrol_EAC_imm_resp_model_df[grep("PreTx_Tumor",Carrol_EAC_imm_resp_model_df$sample_type),]
colnames(Carrol_EAC_imm_resp_model_df)[colnames(Carrol_EAC_imm_resp_model_df) 
                  %in% c("PreTx.Tumor_Mutational_Load","Age","Sex")] <- c("TMB","age","sex")
Carrol_EAC_imm_resp_model_df$sample_type <- NULL
Carrol_EAC_imm_resp_model_df$treatment_response <- ifelse(Carrol_EAC_imm_resp_model_df$treatment_response == "Non-Responder","Non-Responders","Responders")
Carrol_EAC_imm_resp_model_df$cancer_type <- "esophagous"
table(Carrol_EAC_imm_resp_model_df$treatment_response)
# Hartwig
Hartwig_imm_resp_model_df <- data.frame(Hartwig_imm_resp_model[,c("NMDeff_mean","NMD_method","treatment_response","TOTAL_SNV","cancerType","primaryTumorType","ageAtBiopsy","gender")])
colnames(Hartwig_imm_resp_model_df)[colnames(Hartwig_imm_resp_model_df) 
                  %in% c("TOTAL_SNV","cancerType","ageAtBiopsy","gender")] <- c("TMB","cancer_type","age","sex")
Hartwig_imm_resp_model_df$TMB <- Hartwig_imm_resp_model_df$TMB/WES_size
Hartwig_imm_resp_model_df$cancer_type <- factor(tolower(Hartwig_imm_resp_model_df$cancer_type))
table(Hartwig_imm_resp_model_df$treatment_response)
# Split by cancer type
Hartwig_Lung_imm_resp_model_df <- Hartwig_imm_resp_model_df[which(Hartwig_imm_resp_model_df$cancer_type == "lung"),]
Hartwig_UT_imm_resp_model_df <- Hartwig_imm_resp_model_df[which(Hartwig_imm_resp_model_df$cancer_type == "urinary tract"),]
Hartwig_SKCM_imm_resp_model_df <- Hartwig_imm_resp_model_df[which(Hartwig_imm_resp_model_df$cancer_type == "skin" | Hartwig_imm_resp_model_df$primaryTumorType == "Melanoma"),]
Hartwig_all_imm_resp_model_df <- Hartwig_imm_resp_model_df
# Riaz
Riaz_imm_resp_model_df <- data.frame(Riaz_imm_resp_model[,c("NMDeff_mean","NMD_method","treatment_response","Mutation.Load","Subtype","patient_full_ID")])
treatment <- "On"
Riaz_imm_resp_model_df <- Riaz_imm_resp_model_df %>% filter(grepl(treatment,patient_full_ID))
colnames(Riaz_imm_resp_model_df)[colnames(Riaz_imm_resp_model_df) 
                  %in% c("Mutation.Load","Subtype")] <- c("TMB","cancer_subtype")
Riaz_imm_resp_model_df$TMB <- Riaz_imm_resp_model_df$TMB/WES_size
Riaz_imm_resp_model_df$cancer_type <- "skin"
table(Riaz_imm_resp_model_df$treatment_response)
# Merge
# joint_imm_resp_model <- rbind(rbind(TCGA_SKCM_imm_resp_model_df,Liu_SKCM_imm_resp_model_df),Carrol_EAC_imm_resp_model_df,Hartwig_imm_resp_model_df)
# table(joint_imm_resp_model$treatment_response)



# Try cox regression






# Prepare data
# set.seed(111) ## buena
# set.seed(222) ##
# set.seed(333) ## 

all_df_results <- c()

#all_df_results$predicted_responders_perc <- round( (all_df_results$predicted_responders / all_df_results$responders)*100 ,2)

for (dataset in c("TCGA_SKCM","Liu_SKCM","Carrol_EAC","Hartwig_Lung","Riaz_SKCM",
                  "Hartwig_all","Hartwig_SKCM","Hartwig_UT")) {
# for (dataset in c("TCGA_SKCM","Liu_SKCM","Carrol_EAC","Hartwig","all")) {

  if (dataset == "TCGA_SKCM") {
    Immuno_NMDeff_all <- TCGA_SKCM_imm_resp_model_df
  } else if (dataset == "Liu_SKCM") {
    Immuno_NMDeff_all <- Liu_SKCM_imm_resp_model_df
  } else if (dataset == "Carrol_EAC") {
    Immuno_NMDeff_all <- Carrol_EAC_imm_resp_model_df
  } else if (dataset == "Hartwig_all") {
    Immuno_NMDeff_all <- Hartwig_all_imm_resp_model_df
  } else if (dataset == "Hartwig_Lung") {
    Immuno_NMDeff_all <- Hartwig_Lung_imm_resp_model_df
  } else if (dataset == "Hartwig_SKCM") {
    Immuno_NMDeff_all <- Hartwig_SKCM_imm_resp_model_df
  } else if (dataset == "Hartwig_UT") {
    Immuno_NMDeff_all <- Hartwig_UT_imm_resp_model_df
  } else if (dataset == "Riaz_SKCM") {
    Immuno_NMDeff_all <- Riaz_imm_resp_model_df
  } else if (dataset == "all") {
    Immuno_NMDeff_all <- joint_imm_resp_model
  }

  for (NMD_method in c("endogenous_NMD_Consensus", "ASE_PTC_NMD_triggering_0.2")) {
    df_results <- data.frame(dataset = NA, NMD_method = NA, sample_size = NA, accuracy = NA, 
                               AUC = NA, AUC_CI_low = NA, AUC_CI_high = NA, AUC_p_value = NA,
                               specificity = NA, sensitivity = NA, precision = NA,
                              model = NA, responders = NA, non_responders = NA, coefficient = NA, p_value = NA)
    df_results$dataset <- dataset
    if (NMD_method == "endogenous_NMD_Consensus") {
      NMD_method_type <- "ETG"
    } else {
      NMD_method_type <- "ASE"
    }
    df_results$NMD_method <- NMD_method_type
    Immuno_NMDeff_all <- Immuno_NMDeff_all[!is.na(Immuno_NMDeff_all$treatment_response),]

    ####
    # Immuno_NMDeff_all_filt <- unstack(Immuno_NMDeff_all[,c("NMDeff_mean","NMD_method")])
    # cols <- colnames(Immuno_NMDeff_all)[!colnames(Immuno_NMDeff_all) %in% c("NMDeff_mean","NMD_method")]
    # n <- nrow(Immuno_NMDeff_all_filt)
    # for (col in cols) {
    #   Immuno_NMDeff_all_filt[,col] <- Immuno_NMDeff_all[,col][1:n]
    # }
    
    ####

    Immuno_NMDeff_all_filt <- Immuno_NMDeff_all[Immuno_NMDeff_all$NMD_method == NMD_method,]
    training_samples <- Immuno_NMDeff_all_filt$treatment_response %>% 
            createDataPartition(p = 1, list = FALSE)
    train_data  <- Immuno_NMDeff_all_filt[training_samples, ]
    print(dim(train_data))
    df_results$sample_size <- nrow(train_data)
    df_results$responders <- as.numeric(table(train_data$treatment_response=="Responders")[2])
    df_results$non_responders <- as.numeric(table(train_data$treatment_response=="Responders")[1])
    # test_data <- Immuno_NMDeff_all_filt[-training_samples, ]
    # print(dim(test_data))
    contrasts(factor(train_data$treatment_response))
    # # Make sure that the predictor variables are normally distributed. If not, you can use log, root, Box-Cox transformation.
    # Remove highly correlated predictors to minimize overfitting. The presence of highly correlated predictors might lead to an unstable model solution
    for ( model in c("TMB","iNMDeff","both")) {
      if (model == "TMB") {
        # Fit the model
        if (dataset %in% c("Hartwig_all")) {
          model_res <- glm( factor(treatment_response) ~ scale(log(TMB)) + factor(cancer_type) + scale(age) + factor(sex), data = train_data, family = binomial)
        } else if (dataset %in% c("TCGA_SKCM","Carrol_EAC","Hartwig_Lung","Hartwig_SKCM","Hartwig_UT")) {
          model_res <- glm( factor(treatment_response) ~ scale(log(TMB)) + scale(age) + factor(sex), data = train_data, family = binomial)
        } else if (dataset %in% c("Liu_SKCM")) {
          model_res <- glm( factor(treatment_response) ~ scale(log(TMB)) + factor(cancer_type) + factor(sex), data = train_data, family = binomial)
        } else if (dataset %in% c("Riaz_SKCM")) {
          model_res <- glm( factor(treatment_response) ~ scale(log(TMB)) + factor(cancer_subtype), data = train_data, family = binomial)
        }
        model_type <- "TMB"
        df_results <- model_results(model = model_res, df_results = df_results, train_data = train_data)
      } else if (model == "iNMDeff") {
        if (dataset %in% c("Hartwig_all")) {
          model_res <- glm( factor(treatment_response) ~ scale(log(exp(NMDeff_mean))) + factor(cancer_type) + scale(age) + factor(sex), data = train_data, family = binomial)
          # model_res <- glm( factor(treatment_response) ~ scale(log(exp(endogenous_NMD_Consensus))) + scale(log(exp(ASE_PTC_NMD_triggering_0.2))) + factor(cancer_type) + scale(age) + factor(sex), data = train_data, family = binomial)
        } else if (dataset %in% c("TCGA_SKCM","Carrol_EAC","Hartwig_Lung","Hartwig_SKCM","Hartwig_UT")) {
          model_res <- glm( factor(treatment_response) ~ scale(log(exp(NMDeff_mean))) + scale(age) + factor(sex), data = train_data, family = binomial)
          # model_res <- glm( factor(treatment_response) ~ scale(log(exp(endogenous_NMD_Consensus))) + scale(log(exp(ASE_PTC_NMD_triggering_0.2))) + scale(age) + factor(sex), data = train_data, family = binomial)
        } else if (dataset %in% c("Liu_SKCM")) {
          model_res <- glm( factor(treatment_response) ~ scale(log(exp(NMDeff_mean))) + factor(cancer_type) + factor(sex), data = train_data, family = binomial)
          # model_res <- glm( factor(treatment_response) ~ scale(log(exp(endogenous_NMD_Consensus))) + scale(log(exp(ASE_PTC_NMD_triggering_0.2))) + factor(cancer_type) + factor(sex), data = train_data, family = binomial)
        } else if (dataset %in% c("Riaz_SKCM")) {
          model_res <- glm( factor(treatment_response) ~ scale(log(exp(NMDeff_mean))) + factor(cancer_subtype), data = train_data, family = binomial)
          # model_res <- glm( factor(treatment_response) ~ scale(log(exp(endogenous_NMD_Consensus))) + scale(log(exp(ASE_PTC_NMD_triggering_0.2))) + factor(cancer_subtype), data = train_data, family = binomial)
        }
        model_type <- "iNMDeff"
        df_results <- model_results(model = model_res, df_results = df_results, train_data = train_data)
      } else if (model == "both") {
        if (dataset %in% c("Hartwig_all")) {
          model_res <- glm( factor(treatment_response) ~ scale(log(TMB)) + scale(log(exp(NMDeff_mean))) + scale(log(exp(NMDeff_mean))*log(TMB)) + 
                                  factor(cancer_type) + scale(age) + factor(sex), data = train_data, family = binomial)
          # model_res <- glm( factor(treatment_response) ~ scale(log(TMB)) + scale(log(exp(endogenous_NMD_Consensus))) + scale(log(exp(ASE_PTC_NMD_triggering_0.2))) + 
          #                         scale(log(exp(endogenous_NMD_Consensus))*log(TMB)) + scale(log(exp(ASE_PTC_NMD_triggering_0.2))*log(TMB)) + 
          #                         factor(cancer_type) + scale(age) + factor(sex), data = train_data, family = binomial)
        } else if (dataset %in% c("TCGA_SKCM","Carrol_EAC","Hartwig_Lung","Hartwig_SKCM","Hartwig_UT")) {
          model_res <- glm( factor(treatment_response) ~ scale(log(TMB)) + scale(log(exp(NMDeff_mean))) + scale(log(exp(NMDeff_mean))*log(TMB)) + 
                                   scale(age) + factor(sex), data = train_data, family = binomial)
          # model_res <- glm( factor(treatment_response) ~ scale(log(TMB)) + scale(log(exp(endogenous_NMD_Consensus))) + scale(log(exp(ASE_PTC_NMD_triggering_0.2))) + 
          #                           scale(log(exp(endogenous_NMD_Consensus))*log(TMB)) + scale(log(exp(ASE_PTC_NMD_triggering_0.2))*log(TMB)) + 
          #                          scale(age) + factor(sex), data = train_data, family = binomial)
        } else if (dataset %in% c("Liu_SKCM")) {
          model_res <- glm( factor(treatment_response) ~ scale(log(TMB)) + scale(log(exp(NMDeff_mean))) + scale(log(exp(NMDeff_mean))*log(TMB)) +
                                    factor(cancer_type) + factor(sex), data = train_data, family = binomial)
          # model_res <- glm( factor(treatment_response) ~ scale(log(TMB)) + scale(log(exp(endogenous_NMD_Consensus))) + scale(log(exp(ASE_PTC_NMD_triggering_0.2))) + 
          #                           scale(log(exp(endogenous_NMD_Consensus))*log(TMB)) + scale(log(exp(ASE_PTC_NMD_triggering_0.2))*log(TMB)) + 
          #                           factor(cancer_type) + factor(sex), data = train_data, family = binomial)
        } else if (dataset %in% c("Riaz_SKCM")) {
          model_res <- glm( factor(treatment_response) ~ scale(log(TMB)) + scale(log(exp(NMDeff_mean))) + scale(log(exp(NMDeff_mean))*log(TMB)) +
                                    factor(cancer_subtype), data = train_data, family = binomial)
          # model_res <- glm( factor(treatment_response) ~ scale(log(TMB)) + scale(log(exp(endogenous_NMD_Consensus))) + scale(log(exp(ASE_PTC_NMD_triggering_0.2))) + 
          #                           scale(log(exp(endogenous_NMD_Consensus))*log(TMB)) + scale(log(exp(ASE_PTC_NMD_triggering_0.2))*log(TMB)) + 
          #                           factor(cancer_subtype), data = train_data, family = binomial)
        }
        model_type <- "TMB,iNMDeff,TMB*iNMDeff"
        df_results <- model_results(model = model_res, df_results = df_results, train_data = train_data)
      }
      # Save results
      if (length(df_results) == 0) {
        all_df_results <- df_results
      } else {
        all_df_results <- rbind(all_df_results,df_results)
      }
    }
  } # NMD_method
}

model_results <- function(model, df_results, train_data) {
  # Summarize the model
  glm_res <- summary(model)
  remove_rows <- grep("cancer|sex|age",names(glm_res$coefficients[,"Estimate"]))
  df_results$p_value <- paste0(round(as.numeric(paste0(glm_res$coefficients[,"Pr(>|z|)"][-c(1,remove_rows)])),3),collapse=",")
  df_results$coefficient <- paste0(round(as.numeric(paste0(glm_res$coefficients[,"Estimate"][-c(1,remove_rows)])),3),collapse=",")
  df_results$model <- model_type
  # Make predictions
  predicted_probabilities <- model %>% predict(train_data, type = "response")
  predicted_classes <- ifelse(predicted_probabilities > 0.5, "Responders", "Non-Responders")
  table(predicted_classes)
  n_responders <- table(predicted_classes)["Responders"]
  if (is.na(n_responders)) { n_responders <- 0}
  df_results$predicted_responders <- n_responders
  n_non_responders <- table(predicted_classes)["Non-Responders"]
  if (is.na(n_non_responders)) { n_non_responders <- 0}
  df_results$predicted_nonresponders <- n_non_responders
  # 1) Model accuracy
  observed_classes <- train_data$treatment_response
  model_accuracy <- mean(predicted_classes == observed_classes,na.rm=TRUE)
  df_results$accuracy <- round(model_accuracy,2)
  # 2) Sensitivity (Recall) and Specificity
  conf_matrix <- table(predicted_classes,observed_classes)
  check <- (ncol(conf_matrix) == 2) & (nrow(conf_matrix) == 2)
  if (isTRUE(check)) {
    model_sensitivity <- caret::sensitivity(conf_matrix)
    df_results$sensitivity <- model_sensitivity
    model_specificity <- caret::specificity(conf_matrix)
    df_results$specificity <- model_specificity
    # 3) Precision
    model_precision <- caret::precision(conf_matrix)
    df_results$precision <- model_precision
  }
  ##### 4) Observed ROC curve with observed AUC ####
  # Create a ROC object and obtain AUC
  roc_obj <- roc(observed_classes, predicted_probabilities)
  AUC_observed <- auc(roc_obj)
  AUC_CI <- ci.auc(roc_obj)
  df_results$AUC <- round(AUC_observed,2)
  df_results$AUC_CI_low <- round(AUC_CI[1],2)
  df_results$AUC_CI_high <- round(AUC_CI[3],2)
  ##### Random ROC curve for an AUC of 0.5 ####
  # Number of instances
  # n <- length(train_data$treatment_response) # Assuming 100 instances for demonstration
  n <- 100
  # Generate random predictions
  random_predictions <- runif(n, min = 0, max = 1)
  # Generate actual outcomes with an equal number of 0s and 1s
  actual_outcomes <- c(rep(0, n/2), rep(1, n/2))
  # Shuffle the actual outcomes to ensure they are not ordered
  set.seed(555) # Setting seed for reproducibility
  actual_outcomes <- sample(actual_outcomes)
  # Calculate AUC
  roc_obj <- roc(actual_outcomes, random_predictions)
  AUC_random <- auc(roc_obj)
  # Print AUC value
  # print(paste("AUC:", AUC_random))
  ##### Test of Observed vs Random ROC curve ####
  AUC_obs_vs_random_test <- roc.test(AUC_observed, AUC_random, alternative="greater")
  df_results$AUC_p_value <- AUC_obs_vs_random_test$p.value
  # Return
  return(df_results)
}

# Plot

all_df_results$dataset_n <- paste0(all_df_results$dataset, " \n(n=",all_df_results$sample_size,")")
plot <- all_df_results %>%
            filter(!dataset %in% "Hartwig_all") %>%
            # mutate(coefficient = as.numeric(coefficient)) %>%
            ggplot(aes(x = factor(dataset_n), y = AUC, fill = factor(model))) + 
                geom_bar( stat = 'identity', position = position_dodge(width=0.9) ) +
                # Add error bars for confidence intervals
                geom_errorbar(aes(ymin = AUC_CI_low, ymax = AUC_CI_high), 
                               position = position_dodge(width=0.9), width = 0.2) +  # Adjust width as needed
                facet_wrap(NMD_method ~.) + ylim(c(0,1))+
                geom_hline(yintercept=0.5, linetype="dashed", color = "black", size=1) +
                labs(title = "", x = "", y = "AUC", fill = "Model") +                 
                theme_classic(base_size = 35) + scale_fill_brewer(palette = "OrRd", direction = -1) +
                theme(legend.position = "top", 
                      axis.text.x = element_text(angle = 25, hjust = 0.5, vjust = 0.75, size = 17),
                        legend.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
                geom_text(aes(label = ifelse(AUC_p_value < 0.05, "*", "")), 
                                        position = position_dodge(width=0.9), size = 13, color = "black", hjust = 0.5, vjust = 0.5)

output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/phenotypes/immunotherapy/validation_datasets/Immunotherapy_response_logistic_model_all_datasets.png")
# output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/phenotypes/immunotherapy/validation_datasets/ROC_curve.png")
png(output_path, width = 6500, height = 3000, res = 300)
print(plot)
dev.off()
saveRDS(all_df_results, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig6/Fig6G.RData")
output_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/Tables/SuppTableS8.txt"
all_df_results$dataset_n <- NULL
write.table(all_df_results, file = output_path, 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE) 

# ROC curves
# data(ExampleData)
# # specify column number of the outcome variable
# cOutcome <- 2

# # fit logistic regression models
# # all steps needed to construct a logistic regression model are written in a function
# # called 'ExampleModels', which is described on page 4-5
# riskmodel1 <- ExampleModels()$riskModel1
# riskmodel2 <- ExampleModels()$riskModel2

# # obtain predicted risks
# predRisk1 <- predRisk(riskmodel1)
# predRisk2 <- predRisk(riskmodel2)

# # specify label of the ROC curve
# labels <- c("without genetic factors", "with genetic factors")

# # produce ROC curve
# plotROC(data=ExampleData, cOutcome=cOutcome, 
# predrisk=cbind(predRisk1,predRisk2), labels=labels)

# train_data$treatment_response <- ifelse(train_data$treatment_response == "Responders",1,0)
# plot <- plotROC(train_data, 3, model_res$fitted.values) 



# Plot of coefficients
all_df_results[all_df_results$model == "TMB,iNMDeff,TMB*iNMDeff",]
all_df_results_coeffs <- all_df_results
# df <- all_df_results[all_df_results$NMD_method == "ASE" & all_df_results$dataset != "Hartwig_all",]
# df

all_df_results_coeffs <- data.frame(all_df_results_coeffs %>%
  separate_rows(coefficient, p_value, sep = ",") %>%
  mutate(
    coefficient = as.numeric(coefficient),
    p_value = as.numeric(p_value)
  ) %>%
  mutate(term = model)
  )


sequence_vector <- c()

for (start in seq(3, 78, by = 5)) {
  sequence_vector <- c(sequence_vector, seq(start, start+2))
}

all_df_results_coeffs[sequence_vector,"term"] <- rep(c("TMB","iNMDeff","TMB*iNMDeff"),length(sequence_vector)/3)

# x = interaction(model, dataset), 

plot_1 <- all_df_results_coeffs %>%
            filter(NMD_method == "ASE") %>%
            filter(abs(coefficient) < 10) %>%
            filter(!dataset %in% "Hartwig_all") %>%
            # ggplot(aes(x = factor(term), y = coefficient, color = factor(model), fill = factor(dataset) )) + 
            ggplot(aes(x = term, y = coefficient, fill = factor(dataset) )) + 
                geom_bar(stat = 'identity', position = position_dodge(width=.9), size = 1.5) +
                facet_wrap(NMD_method + model ~., scale = "free_y") +# coord_cartesian(ylim = c(-10,10))+
                geom_hline(yintercept=00, linetype="dashed", color = "black", size=1) +
                labs(title = "", x = "Variables in the model", y = "coefficients") +                 
                theme_classic(base_size = 35) + scale_fill_brewer(palette = "Spectral", direction = -1) +
                theme(legend.position = "top", 
                      axis.text.x = element_text(angle = 25, hjust = 0.5, vjust = 0.75, size = 17),
                        legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
                        legend.box = "horizontal",
                        legend.title = element_text(size = 16), # Decrease size of legend title
                        legend.text = element_text(size = 15)) +
                guides(fill = guide_legend(title = "dataset", 
                                          title.position = "top", 
                                          title.theme = element_text(size = 16), 
                                          label.theme = element_text(size = 15)),
                      color = guide_legend(title = "model", 
                                            title.position = "top", 
                                            title.theme = element_text(size = 16), 
                                            label.theme = element_text(size = 15))) +
                geom_text(aes(label = ifelse(p_value <= 0.25 & p_value > 0.1, "*", "")),
                            stat = 'identity', position = position_dodge(width=.9), size = 9, color = "black", hjust = 0.5) +
                geom_text(aes(label = ifelse(p_value <= 0.1 & p_value > 0.05, "**", "")), 
                                        position = position_dodge(width=.9), size = 9, color = "black", hjust = 0.5) +
                geom_text(aes(label = ifelse(p_value <= 0.05, "***", "")), 
                                        position = position_dodge(width=.9), size = 9, color = "black", hjust = 0.5)

plot_2 <- all_df_results_coeffs %>%
            filter(NMD_method == "ETG") %>%
            filter(abs(coefficient) < 10) %>%
            filter(!dataset %in% "Hartwig_all") %>%
            # ggplot(aes(x = factor(term), y = coefficient, color = factor(model), fill = factor(dataset) )) + 
            ggplot(aes(x = term, y = coefficient, fill = factor(dataset) )) + 
                geom_bar(stat = 'identity', position = position_dodge(width=.9), size = 1.5) +
                facet_wrap(NMD_method + model ~., scale = "free_y") +# coord_cartesian(ylim = c(-10,10))+
                geom_hline(yintercept=00, linetype="dashed", color = "black", size=1) +
                labs(title = "", x = "Variables in the model", y = "coefficients") +                 
                theme_classic(base_size = 35) + scale_fill_brewer(palette = "Spectral", direction = -1) +
                theme(legend.position = "top", 
                      axis.text.x = element_text(angle = 25, hjust = 0.5, vjust = 0.75, size = 17),
                        legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
                        legend.box = "horizontal",
                        legend.title = element_text(size = 16), # Decrease size of legend title
                        legend.text = element_text(size = 15)) +
                guides(fill = guide_legend(title = "dataset", 
                                          title.position = "top", 
                                          title.theme = element_text(size = 16), 
                                          label.theme = element_text(size = 15)),
                      color = guide_legend(title = "model", 
                                            title.position = "top", 
                                            title.theme = element_text(size = 16), 
                                            label.theme = element_text(size = 15))) +
                geom_text(aes(label = ifelse(p_value <= 0.25 & p_value > 0.1, "*", "")),
                            stat = 'identity', position = position_dodge(width=.9), size = 9, color = "black", hjust = 0.5) +
                geom_text(aes(label = ifelse(p_value <= 0.1 & p_value > 0.05, "**", "")), 
                                        position = position_dodge(width=.9), size = 9, color = "black", hjust = 0.5) +
                geom_text(aes(label = ifelse(p_value <= 0.05, "***", "")), 
                                        position = position_dodge(width=.9), size = 8, color = "black", hjust = 0.5)

plot <- plot_grid(plot_1,plot_2,nrow = 2,  align = "hv",
                        rel_heights = c(1, 1), axis = "top")

output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/phenotypes/immunotherapy/validation_datasets/Immunotherapy_response_logistic_model_all_datasets_coefficients_pvals.png")
png(output_path, width = 6000, height = 5000, res = 300)
print(plot)
dev.off()
