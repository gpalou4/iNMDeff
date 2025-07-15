#conda activate /home/gpalou/.conda/envs/survival_analysis
library(survival)
library(survminer)
library(forestmodel)
library(dplyr)
library(ggpubr)
library(cowplot)
library(ggh4x)
library(survcomp)
library(qvalue)
library(RColorBrewer)
# library(poolr)
# options(scipen=999)

survival_analysis_plots <- function(TCGA_df, test_variable, percentile, survival_var = "OS", TCGA_cancer_metadata_drug_therapy_all = NULL,
                            cancer_type, randomization = "no", KM_curve = NULL, treatment_type = "all") {
   
  # Cancer type or pancancer
  if (cancer_type != "pancancer") {
    TCGA_df <- TCGA_df[which(TCGA_df$cancer_type_strat %in% cancer_type),]
    char <- "by_cancer"
  } else {
    char <- "pancancer"
  }

  if (survival_var == "OS") {
    survival_var_char <- "days_to_last_follow_up"
    status_var_char <- "status"
  } else if (survival_var == "PFS") {
    survival_var_char <- "PFS_months"
    status_var_char <- "PFS_status"
  }

  # df <- TCGA_df[which(TCGA_df$cancer_type_strat == "SKCM" & TCGA_df$pharmaceutical_therapy_type == "Immunotherapy"),]

  if (treatment_type != "all") {
    if (treatment_type == "no_treatment") {
      TCGA_df <- TCGA_df[!TCGA_df$sample %in% TCGA_cancer_metadata_drug_therapy_all$sample,]
    } else {
      TCGA_df <- merge(TCGA_df,TCGA_cancer_metadata_drug_therapy_all, by.x = "sample", by.y = "sample")
      TCGA_therapy_counts <- TCGA_df %>%
            dplyr::group_by(sample) %>%
            dplyr::summarize(immunotherapy_counts = sum(pharmaceutical_therapy_type == "Immunotherapy"),
                    chemotherapy_counts = sum(pharmaceutical_therapy_type == "Chemotherapy"),
                    other_therapy_counts = sum(!pharmaceutical_therapy_type %in% c("Immunotherapy","Chemotherapy")))
      if (treatment_type == "immunotherapy") {
        samples_to_keep <- as.character(TCGA_therapy_counts[TCGA_therapy_counts$immunotherapy_counts >= 1,"sample"]$sample)
        TCGA_df <- TCGA_df %>%
            filter(sample %in% samples_to_keep) %>%
            filter(pharmaceutical_therapy_type %in% c("Immunotherapy")) #%>%
            # filter(grepl("[Pp]embrolizumab|[Nn]ivolumab|[iI]pilimumab|[Tt]remelimumab|[Cc]emiplimab|[Aa]tezolizumab|[Aa]velumab|[Dd]urvalumab|[Dd]ostarlimab|[Rr]etifanlimab",pharmaceutical_therapy_drug_name))
      } else if (treatment_type == "chemotherapy") {
        samples_to_keep <- as.character(TCGA_therapy_counts[TCGA_therapy_counts$chemotherapy_counts >= 1 & TCGA_therapy_counts$immunotherapy_counts == 0,"sample"]$sample)
        TCGA_df <- TCGA_df %>%
            dplyr::filter(sample %in% samples_to_keep) %>%
            dplyr::filter(pharmaceutical_therapy_type %in% c("Chemotherapy")) %>%
            dplyr::filter(!grepl("[Pp]embrolizumab|[Nn]ivolumab|[iI]pilimumab|[Tt]remelimumab|[Cc]emiplimab|[Aa]tezolizumab|[Aa]velumab|[Dd]urvalumab|[Dd]ostarlimab|[Rr]etifanlimab",pharmaceutical_therapy_drug_name))
      } else if (treatment_type == "all_treatment") {
        TCGA_df <- TCGA_df
      }
      # Remove duplicate samples
      remove_cols <- colnames(TCGA_cancer_metadata_drug_therapy_all)[2:4]
      TCGA_df <- TCGA_df[,-which(colnames(TCGA_df) %in% remove_cols)]
      TCGA_df <- TCGA_df[!duplicated(TCGA_df),]
    }
  }
  # Create DF results
  df_res <- data.frame(exp_HR = NA, p_value = NA, CI_2.5 = NA, CI_97.5 = NA)
  if (nrow(TCGA_df) == 0) {return(df_res)}

  # Randomize NMDeff
  if (randomization == "yes") {
    TCGA_df[,test_variable] <- sample(TCGA_df[,test_variable])
  }

  # Create NMDeff variable high vs low
  TCGA_df$test_variable <- NA
  percentiles_num <- quantile(TCGA_df[,test_variable],probs = seq(0, 1, 0.01), na.rm = TRUE)
  thres_low <- as.numeric(percentiles_num[paste0(percentile,"%")])
  thres_high <- as.numeric(percentiles_num[paste0(100-percentile,"%")])
  TCGA_df[which(TCGA_df[,test_variable] <= thres_low),"test_variable"] <- "Low"
  TCGA_df[which(TCGA_df[,test_variable] >= thres_high),"test_variable"] <- "High"
  ### Fran test
  # thres_high <- as.numeric(percentiles_num[paste0(100-percentile,"%")])
  # TCGA_df[which(TCGA_df[,test_variable] < thres_high),"test_variable"] <- "Low"
  # TCGA_df[which(TCGA_df[,test_variable] >= thres_high),"test_variable"] <- "High"
  ### Fran test

  if ( ( sum(na.omit(TCGA_df$test_variable == "Low")) < 5 ) & ( sum(na.omit(TCGA_df$test_variable == "High")) < 5 ) ) {
    return(df_res)
  } else if (cancer_type == "DLBC") {
    return(df_res)
  }

  surv_obj <<- NULL
  surv_obj <<- with(TCGA_df, Surv(as.numeric(eval(parse(text = survival_var_char))), as.numeric(eval(parse(text = status_var_char)))))
  km_fit <<- survfit(surv_obj ~ test_variable, data = TCGA_df)

  # 1) Kaplan-Meier curves
  # p <- ggsurvplot(km_fit, data = TCGA_df, risk.table = TRUE,
  #           pval = TRUE, conf.int = FALSE,
  #           ylim = c(0, 1), #xlim = c(0, 22), 
  #           xlab = "Time (years)", ylab = "Survival probability",
  #           ggtheme = theme_bw(), surv.median.line = "v", legend.labs = c("High values", "Low values"), 
  #           legend.title = paste0(test_variable)) 
  # KM_plot <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/phenotypes/survival/",char,"/",cancer_type,"_",survival_var,"_Kaplan_Meier_curve_",test_variable,"_perc_",percentile,"_treatment_",treatment_type,".png")
  # # KM_plot <- paste0("/g/strcombio/fsupek_home/gpalou/Manuscript/figures_cowplot/test2.png")
  # png(KM_plot, width = 3500, height = 2500, res = 300)
  # print(p)
  # dev.off()

  if (!is.null(KM_curve)) {
    return(list(km_fit = km_fit, df = TCGA_df, surv_obj = surv_obj))
  }

  # 2) Cox Model
  TCGA_df$test_variable <- relevel(factor(TCGA_df$test_variable), ref = "Low")
  if (cancer_type == "pancancer") {
    TCGA_df$cancer_type_strat <- relevel(TCGA_df$cancer_type_strat, ref = "BRCA_Normal")
  }
  TCGA_df$age_quartiles <- as.factor(ntile(TCGA_df$age_at_diagnosis, 4))
  TCGA_df$CNV_burden_quartiles <- as.factor(ntile(TCGA_df$CNV_burden, 4))
  TCGA_df$cancer_subtype <- as.factor(as.character(TCGA_df$cancer_subtype))
  TCGA_df$stage <- as.factor(TCGA_df$stage)
  stageNA <- ( sum(is.na(TCGA_df$stage)) / nrow(TCGA_df) ) > 0.9
  isGenderUnique1 <- length(table(TCGA_df[!is.na(TCGA_df$test_variable),"gender"])) == 1
  isGenderUnique2 <- length(table(TCGA_df[((!is.na(TCGA_df$test_variable)) & (!is.na(TCGA_df$CNV_burden_quartiles))),"gender"])) == 1
  # isGenderUnique1 <- length(table(TCGA_df[!is.na(TCGA_df[,test_variable]),"gender"])) == 1
  # isGenderUnique2 <- length(table(TCGA_df[((!is.na(TCGA_df[,test_variable])) & (!is.na(TCGA_df$CNV_burden_quartiles))),"gender"])) == 1
  
  # a <- coxph(surv_obj ~  test_variable + age_quartiles + gender + endogenous_purity + 
  #           sample_lib_size + log(TMB) + CNV_burden, data = TCGA_df)
  # summary(a)

  # Cox regression
  coxmodel <- NULL
  error <- NULL
  last_warning <- NULL
  grep_match <- grep("BRCA|UCEC|CESC|UCS|TGCT|OV",cancer_type)
  tryCatch( { 
    coxmodel <- withCallingHandlers(
      # Remove stage variable from cancers where > 90% stages are NAs
      if (isTRUE(stageNA)) {
          if ( length(grep_match) == 1 ) {
            coxph(surv_obj ~ test_variable + age_quartiles + endogenous_purity + sample_lib_size + log(TMB) , data = TCGA_df, iter.max = 10000)
          } else if (cancer_type == "pancancer") {
            coxph(surv_obj ~  test_variable + age_quartiles + gender + cancer_type_strat + endogenous_purity + sample_lib_size + log(TMB), data = TCGA_df, iter.max = 10000)
          } else if (cancer_type %in% c("LAML", "THYM","DLBC")) {
            coxph(surv_obj ~  test_variable + age_quartiles + gender + sample_lib_size + log(TMB), data = TCGA_df, iter.max = 10000)
          } else if (isGenderUnique1 | isGenderUnique2) {
            coxph(surv_obj ~  test_variable + age_quartiles  + endogenous_purity + sample_lib_size + log(TMB), data = TCGA_df, iter.max = 10000)
          } else {
            coxph(surv_obj ~  test_variable + age_quartiles + gender + endogenous_purity + sample_lib_size + log(TMB), data = TCGA_df, iter.max = 10000)
          }
      } else {
        # I will skip "stage" for now
          if ( length(grep_match) == 1 ) {
            coxph(surv_obj ~ test_variable + age_quartiles + endogenous_purity + sample_lib_size + log(TMB) , data = TCGA_df, iter.max = 10000)
          } else if (cancer_type == "pancancer") {
            coxph(surv_obj ~  test_variable + age_quartiles + gender + cancer_type_strat + endogenous_purity + sample_lib_size + log(TMB), data = TCGA_df, iter.max = 10000)
          } else if (cancer_type %in% c("LAML", "THYM","DLBC")) {
            coxph(surv_obj ~  test_variable + age_quartiles + gender + sample_lib_size + log(TMB), data = TCGA_df, iter.max = 10000)
          } else if (isGenderUnique1 | isGenderUnique2) {
            coxph(surv_obj ~  test_variable + age_quartiles  + endogenous_purity + sample_lib_size + log(TMB), data = TCGA_df, iter.max = 10000)
          } else {
            coxph(surv_obj ~  test_variable + age_quartiles + gender + endogenous_purity + sample_lib_size + log(TMB), data = TCGA_df, iter.max = 10000)
          }
          # if ( length(grep_match) == 1 ) {
          #   coxph(surv_obj ~ test_variable + age_quartiles + endogenous_purity + sample_lib_size + stage + cancer_subtype, data = TCGA_df)
          # } else if (cancer_type == "pancancer") {
          #   coxph(surv_obj ~  test_variable + age_quartiles + gender + cancer_type_strat + endogenous_purity + cancer_subtype + sample_lib_size + stage, data = TCGA_df)
          # } else if (cancer_type %in% c("LAML", "THYM","DLBC")) {
          #   coxph(surv_obj ~  test_variable + age_quartiles + gender + sample_lib_size + stage + cancer_subtype, data = TCGA_df)
          # } else if (isGenderUnique1 | isGenderUnique2) {
          #   coxph(surv_obj ~  test_variable + age_quartiles  + endogenous_purity + sample_lib_size + stage + cancer_subtype, data = TCGA_df)
          # } else {
          #   coxph(surv_obj ~  test_variable + age_quartiles + gender + endogenous_purity + sample_lib_size + stage + cancer_subtype, data = TCGA_df)
          # }
        }

        # coxmodel <- coxph(surv_obj ~  test_variable + age_quartiles + gender + cancer_type_strat + endogenous_purity + 
        #         sample_lib_size + log(TMB), data = TCGA_df, iter.max = 10000)

        # coxmodel_res <- summary(coxmodel)
        # CI_lower_95 <- coxmodel_res$conf.int[1,"lower .95"]
        # CI_upper_95 <- coxmodel_res$conf.int[1,"upper .95"]
        # hasTestVarreg <- grep("test_variable",rownames(coxmodel_res$coefficients))
        # if (length(hasTestVarreg) != 0) {
        #     df_res[,"exp_HR"] <- coxmodel_res$coefficients[hasTestVarreg,"exp(coef)"]
        #     df_res[,"p_value"] <- coxmodel_res$coefficients[hasTestVarreg,"Pr(>|z|)"]
        #     df_res[,"CI_2.5"] <- as.numeric(CI_lower_95)
        #     df_res[,"CI_97.5"] <- as.numeric(CI_upper_95)
        # } 
        # df_res


        # table(TCGA_df$stage)
        # df <- subset(TCGA_df, stage %in% c("STAGE I/II (NOS)","STAGE I", "STAGE IA", "STAGE IB"))
        # df <- subset(TCGA_df, stage %in% c("STAGE II","STAGE IIA", "STAGE IIB", "STAGE IIC"))
        # df <- subset(TCGA_df, stage %in% c("STAGE III","STAGE IIIA","STAGE IIIB", "STAGE IIIC"))
        # df <- subset(TCGA_df, stage %in% c("STAGE IV","STAGE IVA","STAGE IVB"))

        # df <- subset(TCGA_df, stage %in% c("STAGE I/II (NOS)","STAGE I", "STAGE IA", "STAGE IB","STAGE II","STAGE IIA", "STAGE IIB", "STAGE IIC"))
        # df <- subset(TCGA_df, stage %in% c("STAGE III","STAGE IIIA","STAGE IIIB", "STAGE IIIC","STAGE IV","STAGE IVA","STAGE IVB"))

        # df$test_variable <- NA
        # percentiles_num <- quantile(df[,test_variable],probs = seq(0, 1, 0.01), na.rm = TRUE)
        # thres_low <- as.numeric(percentiles_num[paste0(percentile,"%")])
        # thres_high <- as.numeric(percentiles_num[paste0(100-percentile,"%")])
        # df[which(df[,test_variable] <= thres_low),"test_variable"] <- "Low"
        # df[which(df[,test_variable] >= thres_high),"test_variable"] <- "High"
        # df$test_variable <- relevel(factor(df$test_variable), ref = "Low")

        # surv_obj <<- with(df, Surv(as.numeric(eval(parse(text = survival_var_char))), as.numeric(eval(parse(text = status_var_char)))))
        # km_fit <<- survfit(surv_obj ~ test_variable, data = df)
        # coxph(surv_obj ~  test_variable + log(TMB) + gender + endogenous_purity + CNV_burden, data = df)


        ,
      warning = function(w) {
            if (grepl("Ran out of iterations and did not converge", conditionMessage(w))) {
              # Take action for non-convergent result directly here
              print("Non-convergence warning detected.")
              error <<- TRUE # Update the flag
            }
            # Optionally, use invokeRestart("muffleWarning") if you want to suppress the warning
          }
        )
      }, error = function(e) {
        print(paste("Error occurred:", conditionMessage(e)))
        error <<- TRUE # Update the flag to indicate an error occurred
      })

  if (isTRUE(error)) {return(df_res)}
  if (is.na(coxmodel$coefficients[1])) {
    if (length(table(is.na(coxmodel$coefficients))) == 1) {return(df_res)}
  }  
  print(summary(coxmodel))
  # Hazard Ratios plot
  # format <- forest_model_format_options(colour = "firebrick", shape = 15, text_size = 5, point_size = 5)
  # p <- forest_model(coxmodel, covariates = c("test_variable", "gender", "age_quartiles", "CNV_burden_quartiles"), format_options = format)
  # Cox_HR_plot <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/phenotypes/survival/",char,"/",cancer_type,"_Cox_reg_Hazard_Ratios_",test_variable,"_perc_",percentile,".png")
  # png(Cox_HR_plot, width = 3500, height = 2500, res = 300)
  # print(p)
  # dev.off()

  #3) Save coefficients
  coxmodel_res <- summary(coxmodel)
  CI_lower_95 <- coxmodel_res$conf.int[1,"lower .95"]
  CI_upper_95 <- coxmodel_res$conf.int[1,"upper .95"]
  hasTestVarreg <- grep("test_variable",rownames(coxmodel_res$coefficients))
  if (length(hasTestVarreg) != 0) {
      df_res[,"exp_HR"] <- coxmodel_res$coefficients[hasTestVarreg,"exp(coef)"]
      df_res[,"p_value"] <- coxmodel_res$coefficients[hasTestVarreg,"Pr(>|z|)"]
      df_res[,"CI_2.5"] <- as.numeric(CI_lower_95)
      df_res[,"CI_97.5"] <- as.numeric(CI_upper_95)
  } 
  return(df_res)
}

FDR_adjust <- function(cox_res) {

  # Choose the best percentile per cancer_type based on p-value
  cancer_best_percentile <- data.frame(cancer = cancer_types, percentile = NA)
  cox_res_final <- c()
  # outliers
  # outliers <- c("UVM","UCS","DLBC","CHOL","ACC","KICH","BRCA_Normal","UCEC_POLE")
  outliers <- c("DLBC","CHOL","ACC","BRCA_Normal","UCEC_POLE")
  cox_res_filt <- cox_res %>% filter(!cancer_type %in% outliers)

  for ( cancer in cancer_types ) {
    print(cancer)
    cox_res_tmp <- cox_res_filt %>% filter(cancer_type == cancer)
    if (nrow(cox_res_tmp) == 0) {next}
    # Q-value correction across percentiles
    #cox_res_tmp[cox_res_tmp$NMD_method == "ASE_PTC_NMD_triggering_0.2","meta_pvalue"] <- NA
    cox_res_tmp <- cox_res_tmp[order(cox_res_tmp$meta_pvalue),]
    # pvals <- as.numeric(na.omit(cox_res_tmp$meta_pvalue))
    # res <- qvalue(pvals, lambda = 0.1, pi0.method = "smoother")
    cox_res_tmp$meta_pvalue_FDR_adjust <- p.adjust(cox_res_tmp$meta_pvalue, method = "fdr")
    best_percentile <- unique(cox_res_tmp[which(cox_res_tmp$meta_pvalue == min(cox_res_tmp$meta_pvalue, na.rm = TRUE)),"percentile"])
    if (length(best_percentile) == 0) {next}
    row <- cancer_best_percentile$cancer %in% cancer
    cancer_best_percentile[row,"percentile"] <- best_percentile
    cox_res_tmp <- cox_res_tmp[cox_res_tmp$percentile == best_percentile,]
    # Save
    if (length(cox_res_final) == 0 ) {
      cox_res_final <- cox_res_tmp
    } else {
      cox_res_final <- rbind(cox_res_final,cox_res_tmp)
    }
  }
  cox_res_final <- cox_res_final[order(cox_res_final$meta_pvalue_FDR_adjust),]
  return(cox_res_final)
}

survival_analysis_by_cancer <- function(TCGA_NMDeff, randomization) {
  cox_res_df_all <- c()
  for (percentile in percentiles) {
    for ( cancer in cancer_types ) {
      print(paste0("----------------- CANCER TYPE -------------------> ",cancer))
      cox_res_tmp <- c()
      for (NMD_method in c("ASE_PTC_NMD_triggering_0.2","endogenous_NMD_Consensus")) {
        cox_res_df <- survival_analysis_plots(TCGA_df = TCGA_NMDeff, test_variable = NMD_method, percentile = percentile, cancer_type = cancer, 
                                              randomization = randomization, KM_curve = NULL, survival_var = "OS")
        cox_res_df$cancer_type <- cancer
        cox_res_df$percentile <- percentile
        cox_res_df$NMD_method <- NMD_method
        # Outlier removal to avoid a biased meta p-value
        # outlier_rows <- which(cox_res_df$CI_2.5 > 100 | cox_res_df$CI_97.5 > 100 | cox_res_df$exp_HR == Inf | cox_res_df$exp_HR > 100 | cox_res_df$exp_HR < 0.001)
        # outlier_rows <- which(cox_res_df$CI_2.5 > 5 | cox_res_df$CI_97.5 > 10 | cox_res_df$exp_HR == Inf | cox_res_df$exp_HR > 10 | cox_res_df$exp_HR < 0.01)
        # cox_res_df[outlier_rows, 1:4] <- NA
        if (length(cox_res_tmp) == 0 ) {
          cox_res_tmp <- cox_res_df
        } else {
          cox_res_tmp <- rbind(cox_res_tmp,cox_res_df)
        }
      }
      # Meta p-value between ASE & ETG for the same cancer type and percentile
      p_values <- cox_res_tmp$p_value
      exp_HR <- cox_res_tmp$exp_HR
      NA_pvals <- sum(is.na(p_values))
      if ( NA_pvals == 2) {
        final_meta_pvalue <- NA
      } else if ( NA_pvals == 1) {
        final_meta_pvalue <- as.numeric(na.omit(p_values))
      } else {
        # Check that direction of effect is the same in both NMD_methods
        HR_direction <- sum(exp_HR > 1)
        if (HR_direction %in% c(0,2)) {
          final_meta_pvalue <- combine.test(cox_res_tmp$p_value, method = c("fisher"), hetero = FALSE, na.rm = FALSE)
          # "the by-the-book way is to run, instead of the usual 2-tailed test, two separate 1-tailed test.  
          # Then meta-analysis left tails with left tails, and separately right tails with right tails."
        } else {
          final_meta_pvalue <- NA
        }
      }
      cox_res_tmp$meta_pvalue <- final_meta_pvalue
      print(cox_res_tmp)
      # Save results
      if (length(cox_res_df_all) == 0 ) {
        cox_res_df_all <- cox_res_tmp
      } else {
        cox_res_df_all <- rbind(cox_res_df_all,cox_res_tmp)
      }
    }
  }
  return(cox_res_df_all)
}

survival_analysis_by_cancer_stratified_by_treatment <- function(TCGA_NMDeff, randomization, treatment_type, survival_var) {
  cox_res_df_all <- c()
  for ( cancer in cancer_types ) {
    for (percentile in percentiles) {
      print(paste0("----------------- CANCER TYPE -------------------> ",cancer))
      cox_res_tmp <- c()
      for (NMD_method in c("ASE_PTC_NMD_triggering_0.2","endogenous_NMD_Consensus")) {
        cox_res_df <- survival_analysis_plots(TCGA_df = TCGA_NMDeff, test_variable = NMD_method, percentile = percentile, cancer_type = cancer, 
                                              randomization = randomization, KM_curve = NULL, survival_var = survival_var, treatment_type = treatment_type,
                                              TCGA_cancer_metadata_drug_therapy_all = TCGA_cancer_metadata_drug_therapy_all)
        cox_res_df$cancer_type <- cancer
        cox_res_df$percentile <- percentile
        cox_res_df$NMD_method <- NMD_method
        # Outlier removal to avoid a biased meta p-value
        # outlier_rows <- which(cox_res_df$CI_2.5 > 100 | cox_res_df$CI_97.5 > 100 | cox_res_df$exp_HR == Inf | cox_res_df$exp_HR > 100 | cox_res_df$exp_HR < 0.001)
        # outlier_rows <- which(cox_res_df$CI_2.5 > 10 | cox_res_df$CI_97.5 > 15 | cox_res_df$exp_HR == Inf | cox_res_df$exp_HR > 15 | cox_res_df$exp_HR < 0.01)
        # cox_res_df[outlier_rows, 1:4] <- NA
        if (length(cox_res_tmp) == 0 ) {
          cox_res_tmp <- cox_res_df
        } else {
          cox_res_tmp <- rbind(cox_res_tmp,cox_res_df)
        }
      }
      # Meta p-value between ASE & ETG for the same cancer type and percentile
      p_values <- cox_res_tmp$p_value
      exp_HR <- cox_res_tmp$exp_HR
      NA_pvals <- sum(is.na(p_values))
      if ( NA_pvals == 2) {
        final_meta_pvalue <- NA
      } else if ( NA_pvals == 1) {
        final_meta_pvalue <- as.numeric(na.omit(p_values))
      } else {
        # Check that direction of effect is the same in both NMD_methods
        HR_direction <- sum(exp_HR > 1)
        if (HR_direction %in% c(0,2)) {
          final_meta_pvalue <- combine.test(cox_res_tmp$p_value, method = c("fisher"), hetero = FALSE, na.rm = FALSE)
          # "the by-the-book way is to run, instead of the usual 2-tailed test, two separate 1-tailed test.  
          # Then meta-analysis left tails with left tails, and separately right tails with right tails."
        } else {
          final_meta_pvalue <- NA
        }
      }
      cox_res_tmp$meta_pvalue <- final_meta_pvalue
      print(cox_res_tmp)
      # Save results
      if (length(cox_res_df_all) == 0 ) {
        cox_res_df_all <- cox_res_tmp
      } else {
        cox_res_df_all <- rbind(cox_res_df_all,cox_res_tmp)
      }
    }
  }
  return(cox_res_df_all)
}

# Function to count number of p-values below each threshold
count_below_FDR_threshold <- function(FDR) {
  n_hits_obs <- sum(cox_res_non_rand_Fadj$meta_pvalue_FDR_adjust < FDR, na.rm = TRUE)
  n_size_obs <- sum(!is.na(cox_res_non_rand_Fadj$meta_pvalue_FDR_adjust))
  n_hits_rand <- sum(cox_res_rand_Fadj$meta_pvalue_FDR_adjust < FDR, na.rm = TRUE)
  n_size_rand <- sum(!is.na(cox_res_rand_Fadj$meta_pvalue_FDR_adjust))
  data.frame(n_hits_obs,n_size_obs, n_hits_rand, n_size_rand)
}

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

endogenous_NMD_genesets <-  c("endogenous_NMD_Colombo","endogenous_NMD_Karousis","endogenous_NMD_Tani","endogenous_NMD_Courtney","endogenous_NMD_ensembl",
                      "endogenous_NMD_all","endogenous_NMD_Consensus","endogenous_SMG6","endogenous_SMG7",
                      "endogenous_non_NMD_neg_control","endogenous_non_NMD_neg_control_with_NMD_features")
ASE_NMD_genesets <- c("ASE_PTC_NMD_triggering_0.01","ASE_PTC_NMD_evading_0.01","ASE_synonymous_0.01",
                      "ASE_PTC_NMD_triggering_0.2","ASE_PTC_NMD_evading_0.2","ASE_synonymous_0.2")

# 1) TCGA samples NMDeff
sample_NMD_efficiencies_TCGA_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/cancers/pancancer/NMD_efficiencies_TCGA.txt"
sample_NMD_efficiencies_TCGA <- read.table(file = sample_NMD_efficiencies_TCGA_path, header = TRUE, sep = "\t")
sample_NMD_efficiencies_TCGA <- modify_NMDeff_dataframe(sample_NMD_efficiencies_TCGA, dataset = "TCGA", scale = FALSE)
cancer_types_original <- as.character(unique(sample_NMD_efficiencies_TCGA$cancer_type))

# 2) TCGA metadata
TCGA_metadata <- read.table("/g/strcombio/fsupek_cancer1/gpalou/TCGA_metadata/TCGA_clinical_all.tsv", 
                            header = TRUE, sep = "\t", stringsAsFactors = FALSE)
TCGA_metadata <- TCGA_metadata[,colnames(TCGA_metadata) %in% c("submitter_id", "project_id", "gender", "race", "age_at_diagnosis", "vital_status", "days_to_death", "days_to_last_follow_up")]
colnames(TCGA_metadata) <- c("sample", "cancer_type","gender", "race", "age_at_diagnosis", "vital_status", "days_to_death", "days_to_last_follow_up")
# TCGA-SKCM immunotherapy data
# df <- GDCquery_clinic(project = "TCGA-SKCM", type = "clinical", save.csv = FALSE)
# df$treatments_radiation_treatment_type
# grep("Radiation Th",df)

TCGA_cancer_metadata_PFS_all <- c()
TCGA_cancer_metadata_drug_therapy_all <- c()

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
  # PFS metadata
  TCGA_cancer_metadata_PFS <- read.table(paste0("/g/strcombio/fsupek_cancer1/gpalou/TCGA_metadata/TCGA_cBioPortal/",cancer,"_tcga_pan_can_atlas_2018_clinical_data.tsv"), 
                              header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  TCGA_cancer_metadata_PFS <- TCGA_cancer_metadata_PFS[,c("Patient.ID","Overall.Survival..Months.","Overall.Survival.Status","Progress.Free.Survival..Months.","Progression.Free.Status","Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code")]
  colnames(TCGA_cancer_metadata_PFS) <- c("sample","OS_months","OS_status","PFS_months","PFS_status","stage")
  if (length(TCGA_cancer_metadata_PFS_all) == 0) {
    TCGA_cancer_metadata_PFS_all <- TCGA_cancer_metadata_PFS
  } else {
    TCGA_cancer_metadata_PFS_all <- rbind(TCGA_cancer_metadata_PFS_all,TCGA_cancer_metadata_PFS)
  }
}

TCGA_cancer_metadata_PFS_all <- TCGA_cancer_metadata_PFS_all[-which(duplicated(TCGA_cancer_metadata_PFS_all)),]
# TCGA_cancer_metadata_drug_therapy_all <- TCGA_cancer_metadata_drug_therapy_all[-which(duplicated(TCGA_cancer_metadata_drug_therapy_all)),]
# df <- merge(TCGA_metadata,TCGA_cancer_metadata_drug_therapy_all)

# 3) Merge and add survival variables
TCGA_metadata$cancer_type <- gsub("TCGA-","",TCGA_metadata$cancer_type)
sample_NMD_efficiencies_TCGA[,c("vital_status","days_to_death","days_to_last_follow_up")] <- NULL
TCGA_NMDeff <- merge(sample_NMD_efficiencies_TCGA, TCGA_metadata, by = c("sample","cancer_type"), all.x = TRUE)
TCGA_NMDeff <- merge(TCGA_NMDeff,TCGA_cancer_metadata_PFS_all, by = "sample", all.x = TRUE)
TCGA_NMDeff[which(TCGA_NMDeff$vital_status == "dead"),]$days_to_last_follow_up <- TCGA_NMDeff[which(TCGA_NMDeff$vital_status == "dead"),]$days_to_death
TCGA_NMDeff$status <- NA
TCGA_NMDeff[which(TCGA_NMDeff$vital_status == "dead"),]$status <- 1
TCGA_NMDeff[which(TCGA_NMDeff$vital_status == "alive"),]$status <- 0
TCGA_NMDeff$days_to_last_follow_up <- as.numeric(TCGA_NMDeff$days_to_last_follow_up)/365.25
TCGA_NMDeff$age_at_diagnosis <- as.numeric(TCGA_NMDeff$age_at_diagnosis)/365.25
# PFS
TCGA_NMDeff[which(TCGA_NMDeff$PFS_status == "1:PROGRESSION"),]$PFS_status <- 1
TCGA_NMDeff[which(TCGA_NMDeff$PFS_status == "0:CENSORED"),]$PFS_status <- 0
### Subtype_Selected ###
# tags <- as.character(TCGA_NMDeff[which(TCGA_NMDeff$cancer_type %in% c("GBM","LGG")),"Subtype_Selected"])
# tags <- gsub(".*\\.","",tags)
# add <- TCGA_NMDeff[which(TCGA_NMDeff$cancer_type %in% c("GBM","LGG")),"cancer_type"]
# tags <- paste0(add,"_",tags)
# TCGA_NMDeff$cancer_type_strat <- as.character(TCGA_NMDeff$cancer_type_strat)
# TCGA_NMDeff[which(TCGA_NMDeff$cancer_type %in% c("GBM","LGG")),"cancer_type_strat"] <- tags
# TCGA_NMDeff$cancer_type_strat <- factor(TCGA_NMDeff$cancer_type_strat)

### MSI, MSS and POLE ###
# tags <- TCGA_NMDeff[which(TCGA_NMDeff$cancer_type %in% c("TCGA-UCEC","TCGA-COAD","TCGA-STAD")),"Subtype_Selected"]
# tags <- gsub(".*\\.","",tags)
# add <- gsub("TCGA-","",TCGA_NMDeff[which(TCGA_NMDeff$cancer_type %in% c("TCGA-UCEC","TCGA-COAD","TCGA-STAD")),"cancer_type"])
# tags <- paste0(add,"_",tags)
# tags <- gsub("_CIN","",tags)
# tags <- gsub("_GS","",tags)
# tags <- gsub("_HM-SNV","",tags)
# tags <- gsub("_NA","",tags)
# tags <- gsub("_EBV","",tags)
# tags <- gsub("_CN_LOW","",tags)
# tags <- gsub("_CN_HIGH","",tags)

# TCGA_NMDeff[which(TCGA_NMDeff$cancer_type %in% c("TCGA-UCEC","TCGA-COAD","TCGA-STAD")),"cancer_type_strat"] <- tags
# table(TCGA_NMDeff$cancer_type_strat)
# table(TCGA_NMDeff[-which(is.na(TCGA_NMDeff$ASE_stopgain_0.2) | TCGA_NMDeff$ASE_num_PTCs_0.2 < 3),"cancer_type_strat"])

###################################################
######## SURVIVAL ANALYSIS ON ALL SAMPLES #########
###################################################

cancer_types <- c(levels(TCGA_NMDeff$cancer_type_strat),"pancancer") 
percentiles <- c(5,10,15,20,25,30,35,40,45,50)

# 4) Survival analysis for NMDeff
cox_res_non_rand_all <- survival_analysis_by_cancer(TCGA_NMDeff = TCGA_NMDeff, randomization = "no")
outlier_rows <- which(cox_res_non_rand_all$CI_2.5 > 20| cox_res_non_rand_all$CI_97.5 > 20 | 
                    cox_res_non_rand_all$exp_HR == Inf | cox_res_non_rand_all$exp_HR > 15 | cox_res_non_rand_all$exp_HR < 0.001)
cox_res_non_rand <- cox_res_non_rand_all[-outlier_rows,]
# Remove meta-pvalues that actually comes from only 1 method, that are already removed previously
rem_rows <- which(( is.na(cox_res_non_rand$exp_HR) & !is.na(cox_res_non_rand$meta_pvalue) ))
cox_res_non_rand <- cox_res_non_rand[-rem_rows,]

cox_res_rand_all <- c()
for (rand_iteration in 1:10) {
  cox_res_rand_iteration <- survival_analysis_by_cancer(TCGA_NMDeff = TCGA_NMDeff, randomization = "yes")
  cox_res_rand_iteration$rand_iteration <- rand_iteration
  if (length(cox_res_rand_all) == 0 ) {
    cox_res_rand_all <- cox_res_rand_iteration
  } else {
    cox_res_rand_all <- rbind(cox_res_rand_all,cox_res_rand_iteration)
  }
}
outlier_rows <- which(cox_res_rand_all$CI_2.5 > 20| cox_res_rand_all$CI_97.5 > 20 | 
                    cox_res_rand_all$exp_HR == Inf | cox_res_rand_all$exp_HR > 15 | cox_res_rand_all$exp_HR < 0.001)
cox_res_rand <- cox_res_rand_all[-outlier_rows,]
# Remove meta-pvalues that actually comes from only 1 method, that are already removed previously
rem_rows <- which(( is.na(cox_res_rand_all$exp_HR) & !is.na(cox_res_rand_all$meta_pvalue) ))
cox_res_rand_all <- cox_res_rand_all[-rem_rows,]

# cox_res_final_rand <- c()
# for (rand_iteration in 1:10) {
#   cox_res_rand_tmp <- cox_res_rand[cox_res_rand$rand_iteration == rand_iteration,]
#   cox_res_final_rand_tmp <- FDR_adjust(cox_res_rand_tmp)
#   cox_res_final_rand_tmp$meta_pvalue_pancancer_FDR_adjust <- p.adjust(cox_res_final_rand_tmp$meta_pvalue, method = "fdr")
#   if (length(cox_res_final_rand) == 0) {
#     cox_res_final_rand <- cox_res_final_rand_tmp 
#   } else {
#     cox_res_final_rand <- rbind(cox_res_final_rand, cox_res_final_rand_tmp)
#   }
# }
# cox_res_final_non_rand <- FDR_adjust(cox_res_non_rand)

# 5) Same as RVAS analysis
# Take significant hits for each cancer type at different % FDR thresholds, 
# then, calculate empirical FDR based on hits on randomization for each FDR threshold

cox_res_non_rand_Fadj <- cox_res_non_rand %>%
    group_by(cancer_type) %>%
    mutate(meta_pvalue_FDR_adjust = p.adjust(meta_pvalue, method = "fdr"))

cox_res_rand_Fadj <- cox_res_rand %>%
    group_by(cancer_type,rand_iteration) %>%
    mutate(meta_pvalue_FDR_adjust = p.adjust(meta_pvalue, method = "fdr"))

cox_res_non_rand_Fadj[cox_res_non_rand_Fadj$NMD_method == "ASE_PTC_NMD_triggering_0.2","meta_pvalue_FDR_adjust"] <- NA
cox_res_rand_Fadj[cox_res_rand_Fadj$NMD_method == "ASE_PTC_NMD_triggering_0.2","meta_pvalue_FDR_adjust"] <- NA
cox_res_non_rand_Fadj <- data.frame(cox_res_non_rand_Fadj[order(cox_res_non_rand_Fadj$meta_pvalue_FDR_adjust),])
df <- data.frame(cox_res_rand_Fadj[which(cox_res_rand_Fadj$meta_pvalue_FDR_adjust < 0.05),])
data.frame(cox_res_non_rand_Fadj[which(cox_res_non_rand_Fadj$meta_pvalue_FDR_adjust < 0.05),])

# Save
output_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/phenotypes/survival/cox_res_non_rand.txt"
# write.table(cox_res_non_rand_Fadj, file = output_path, 
#                 sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
cox_res_non_rand_Fadj <- read.table(file = output_path, header = TRUE, sep = "\t")
# Create Supp Table S3
cox_res_non_rand_Fadj_SuppTableS3 <- cox_res_non_rand_Fadj %>%
    filter( (cancer_type == "SKCM" & percentile %in% seq(15,45,5) ) |
    ( cancer_type == "ESCA_ac" & percentile %in% seq(30,45,5) ) |
    ( cancer_type == "PRAD" & percentile %in% seq(35,50,5) ) |
    ( cancer_type == "LUAD" & percentile %in% seq(5,10,5) ) ) %>%
    mutate(NMD_method = ifelse(NMD_method == "endogenous_NMD_Consensus","ETG","ASE"),
          outcome = "OS",
          treatment = "any") %>%
    select(cancer_type,percentile,NMD_method,outcome,treatment,exp_HR,CI_2.5,CI_97.5,p_value,meta_pvalue,meta_pvalue_FDR_adjust) %>%
    arrange(cancer_type,percentile,NMD_method)
colnames(cox_res_non_rand_Fadj_SuppTableS3) <- c("cancer_type","percentile","NMD_method","outcome","treatment","exp_HR",
                              "CI_low","CI_high","p_value","meta_p_value","FDR_meta_p_value")
output_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/Tables/SuppTableS3.txt"
write.table(cox_res_non_rand_Fadj_SuppTableS3, file = output_path, 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)   

output_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/phenotypes/survival/cox_res_rand.txt"
# write.table(cox_res_rand_Fadj, file = output_path, 
#                 sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
cox_res_rand_Fadj <- read.table(file = output_path, header = TRUE, sep = "\t")
# write.table(cox_res_non_rand_Fadj, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig6/Fig6A_B.txt", 
#                 sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
# saveRDS(cox_res_non_rand_Fadj, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig6/Fig6A_B.RData")

cancer_type <- "PRAD"
filter <- which(cox_res_non_rand_Fadj$cancer_type == cancer_type & cox_res_non_rand_Fadj$meta_pvalue_FDR_adjust < 0.05)
percentiles_to_check <- cox_res_non_rand_Fadj[filter,"percentile"]
percentiles_to_check <- percentiles_to_check
sort(percentiles_to_check)
df <- cox_res_non_rand_Fadj[which(cox_res_non_rand_Fadj$cancer_type == cancer_type & cox_res_non_rand_Fadj$percentile %in% percentiles_to_check),]
df <- data.frame(df[order(df$percentile),])
median(df$meta_pvalue_FDR_adjust,na.rm=TRUE)
median(df$exp_HR,na.rm=TRUE)

# Apply the function to each threshold
FDR_thresholds <- c(0.03,0.05, 0.1, 0.15, 0.2, 0.25)
n_hits <- sapply(FDR_thresholds, count_below_FDR_threshold)
# Convert to dataframe
n_hits_df <- data.frame(t(n_hits))
for (col in colnames(n_hits_df)) {
  n_hits_df[,col] <- unlist(n_hits_df[,col])

}
n_hits_df$FDR_thresholds <- FDR_thresholds*100
n_hits_df <- n_hits_df %>%
    mutate(obs_hits_ratio = n_hits_obs / n_size_obs,
          rand_hits_ratio = n_hits_rand / n_size_rand,
          FDR_empirical = rand_hits_ratio/obs_hits_ratio)
n_hits_df <- round(n_hits_df,2)
n_hits_stacked <- stack(n_hits_df[,c("obs_hits_ratio","rand_hits_ratio")])
n_hits_stacked$FDR_thresholds <- rep(n_hits_df$FDR_thresholds,2)
n_hits_stacked$FDR_empirical <- rep(n_hits_df$FDR_empirical,2)
colnames(n_hits_stacked) <- c("hits_ratio","type","FDR_thresholds","FDR_empirical")

sig_hits_plot <- n_hits_stacked %>%
            ggplot(aes(x = factor(FDR_thresholds), y = hits_ratio, fill = factor(type))) + 
                geom_bar( stat = 'identity', position = position_dodge() ) +
                geom_text(data = . %>% filter(type == "obs_hits_ratio"),size = 10,
                          aes(label = scales::percent(FDR_empirical, accuracy = 0.01), 
                                          y = hits_ratio + 0.005, # A small adjustment to place the text slightly above the bar
                                          group = factor(FDR_thresholds)), 
                                      position = position_dodge(width = 0.9),
                                      vjust = -0.5) + # Adjust this as needed to position the text correctly
                labs(title = "", x = "% FDR threshold", y = "% significat hits", fill = "") +                 
                theme_classic(base_size = 35) + scale_fill_brewer(palette = "OrRd", direction = -1) +
                theme(legend.position = "top", 
                        legend.margin = margin(t = 0, r = 0, b = 0, l = 0))

sig_hits_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/phenotypes/survival/cox_sig_hits_ratio.png")
png(sig_hits_path, width = 6500, height = 4500, res = 300)
print(sig_hits_plot)
dev.off()

# 5) Final plot for the Cox regressions

# Order
# filter <- cox_res_final$NMD_method == "endogenous_NMD_Consensus"
# tmp <-  cox_res_final[filter,]
# tissue_order <- tmp[order(tmp$exp_HR),"cancer_type"]
# cox_res_final$cancer_type <- factor(cox_res_final$cancer_type, levels = tissue_order)

# # Aggregate the data
# agg_data <- cox_res_final %>%
#   group_by(cancer_type) %>%
#   summarize(significant = any(meta_pvalue_FDR_adjust < 0.05),
#             marginal_sig = any(meta_pvalue_FDR_adjust >= 0.05 & meta_pvalue_FDR_adjust < 0.25),
#             exp_HR = max(exp_HR,na.rm=TRUE)) %>%
#   mutate(pval_sig = ifelse(significant == TRUE,"***",ifelse(marginal_sig == TRUE, "*","")))


# plot_cox <- cox_res_final %>%
#           mutate(pval_sig = ifelse(meta_pvalue_FDR_adjust < 0.05,"***",ifelse(meta_pvalue_FDR_adjust >= 0.05 & meta_pvalue_FDR_adjust < 0.25, "*",""))) %>%
#           ggplot(aes(y = exp_HR, x = cancer_type, group = factor(NMD_method), color = factor(NMD_method), fill = factor(pval_sig))) +
#             geom_point(size = 10, shape = 16, position = position_dodge(width = 0.3)) +
#             # facet_wrap( ~ NMD_method) +
#             geom_errorbar(aes(ymin = CI_2.5, ymax = CI_97.5), size = 2, width = 0.2, position = position_dodge(width = 0.3)) +
#             geom_hline(yintercept = 1, color = "red", linetype = "dashed", size = 0.8, alpha = 0.5) +
#             labs(x = "", y = "Harzard Ratio") +   
#             scale_color_brewer(palette = "Paired",
#                             breaks = c("ASE_PTC_NMD_triggering_0.2", "endogenous_NMD_Consensus"),
#                             labels = c("ASE", "ETG")) +
#             scale_fill_manual(values = c("red","black","white"),breaks = c("***", "*",""), labels = c("<0.05 (***)", "<0.25 (*)","")) +
#             theme_classic(base_size = 25) + 
#             theme(plot.title = element_text(hjust = 0.5, size = 55),
#                     legend.position = "top",
#                     axis.text.y = element_text(size = 35),
#                     axis.title.x = element_text(size = 45),
#                     axis.text.x = element_text(size= 35, angle = 45, hjust = 1, vjust = 0.85),
#                     strip.text = element_text(size = 45),
#                     panel.spacing = grid::unit(3, "lines"),
#                     legend.title = element_text(size = 45),
#                     legend.text = element_text(size = 40)) +
#             geom_text(data = agg_data, aes(x = cancer_type, y = exp_HR, 
#                                           label = ifelse(significant, "***", "")),
#                       hjust = 0.5, vjust = -0.5, size = 20, color = "red") +
#             geom_text(data = agg_data, aes(x = cancer_type, y = exp_HR, 
#                                           label = ifelse(marginal_sig, "*", "")), 
#                       hjust = 0.5, vjust = -0.5, size = 20, color = "black") +
#             guides(color = guide_legend(title = "NMD_method", override.aes = list(size = 16)),
#                   fill = guide_legend(title = "significant", override.aes = list(size = 16))) +
#             coord_cartesian(ylim=c(0,4))

# cox_fores_plot <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/phenotypes/survival/cox_final_plot_best_percentiles.png")
# png(cox_fores_plot, width = 13500, height = 6500, res = 300)
# print(plot_cox)
# dev.off()


##############################################################################
######## SURVIVAL ANALYSIS OF PFS ON STRAFIFIED SAMPLES BY TREATMENT #########
##############################################################################

# For KM plots
# change treatment to chemotherapy/immunotherapy
for (percentile in percentiles) {
  for ( cancer in cancer_types ) {
    print(paste0("----------------- CANCER TYPE -------------------> ",cancer))
    cox_res_tmp <- c()
    for (NMD_method in c("ASE_PTC_NMD_triggering_0.2","endogenous_NMD_Consensus")) {
      df <- survival_analysis_plots(TCGA_df = TCGA_NMDeff, test_variable = NMD_method, percentile = percentile, cancer_type = cancer, TCGA_cancer_metadata_drug_therapy_all = TCGA_cancer_metadata_drug_therapy_all,
                                            randomization = "no", KM_curve = "yes", survival_var = "PFS", treatment_type = "immunotherapy")
    }
  }
}

#############################################
######## TREATMENT --> CHEMOTHERAPY #########
#############################################

cox_res_non_rand_chemo_all <- survival_analysis_by_cancer_stratified_by_treatment(TCGA_NMDeff = TCGA_NMDeff, randomization = "no", 
                                            treatment_type = "chemotherapy", survival_var = "PFS")
outlier_rows <- which(cox_res_non_rand_chemo_all$CI_2.5 > 20| cox_res_non_rand_chemo_all$CI_97.5 > 20 | 
                    cox_res_non_rand_chemo_all$exp_HR == Inf | cox_res_non_rand_chemo_all$exp_HR > 15 | cox_res_non_rand_chemo_all$exp_HR < 0.001)
cox_res_non_rand_chemo <- cox_res_non_rand_chemo_all[-outlier_rows,]
# Remove meta-pvalues that actually comes from only 1 method, that are already removed previously
rem_rows <- which(( is.na(cox_res_non_rand_chemo$exp_HR) & !is.na(cox_res_non_rand_chemo$meta_pvalue) ))
cox_res_non_rand_chemo <- cox_res_non_rand_chemo[-rem_rows,]

### Randomization ###
# cox_res_rand_chemo <- c()
# for (rand_iteration in 1:10) {
#   cox_res_rand_iteration <- survival_analysis_by_cancer_stratified_by_treatment(TCGA_NMDeff = TCGA_NMDeff, randomization = "yes", treatment_type = "chemotherapy")
#   cox_res_rand_iteration$rand_iteration <- rand_iteration
#   if (length(cox_res_rand) == 0 ) {
#     cox_res_rand_chemo <- cox_res_rand_iteration
#   } else {
#     cox_res_rand_chemo <- rbind(cox_res_rand_chemo,cox_res_rand_iteration)
#   }
# }

# Same as RVAS analysis
# Take significant hits for each cancer type at different % FDR thresholds, 
# then, calculate empirical FDR based on hits on randomization for each FDR threshold

cox_res_non_rand_chemo_Fadj <- cox_res_non_rand_chemo %>%
    group_by(cancer_type) %>%
    mutate(meta_pvalue_FDR_adjust = p.adjust(meta_pvalue, method = "fdr"))

# cox_res_rand_chemo_Fadj <- cox_res_rand_chemo %>%
#     group_by(cancer_type,rand_iteration) %>%
#     mutate(meta_pvalue_FDR_adjust = p.adjust(meta_pvalue, method = "fdr"))

# data.frame(cox_res_non_rand_chemo_Fadj[cox_res_non_rand_chemo_Fadj$cancer_type == "SKCM",])

cox_res_non_rand_chemo_Fadj[cox_res_non_rand_chemo_Fadj$NMD_method == "ASE_PTC_NMD_triggering_0.2","meta_pvalue_FDR_adjust"] <- NA
# cox_res_rand_chemo_Fadj[cox_res_rand_chemo_Fadj$NMD_method == "ASE_PTC_NMD_triggering_0.2","meta_pvalue_FDR_adjust"] <- NA
cox_res_non_rand_chemo_Fadj <- data.frame(cox_res_non_rand_chemo_Fadj[order(cox_res_non_rand_chemo_Fadj$meta_pvalue_FDR_adjust),])
data.frame(cox_res_non_rand_chemo_Fadj[which(cox_res_non_rand_chemo_Fadj$meta_pvalue_FDR_adjust < 0.25),])

# Percentiles sig
# cancer_type <- "pancancer"
# filter <- which(cox_res_non_rand_chemo_Fadj$cancer_type == cancer_type & cox_res_non_rand_chemo_Fadj$meta_pvalue_FDR_adjust < 1)
# percentiles_to_check <- cox_res_non_rand_chemo_Fadj[filter,"percentile"]
# percentiles_to_check <- percentiles_to_check$percentile
# percentiles_to_check
# df <- cox_res_non_rand_chemo_Fadj[which(cox_res_non_rand_chemo_Fadj$cancer_type == cancer_type & cox_res_non_rand_chemo_Fadj$percentile %in% percentiles_to_check),]
# df <- data.frame(df[order(df$percentile),])
# median(df$meta_pvalue_FDR_adjust,na.rm=TRUE)
# median(df$exp_HR,na.rm=TRUE)

cancer_type <- "ESCA_scc"
filter <- which(cox_res_non_rand_chemo_Fadj$cancer_type == cancer_type & cox_res_non_rand_chemo_Fadj$meta_pvalue_FDR_adjust < 0.1)
percentiles_to_check <- cox_res_non_rand_chemo_Fadj[filter,"percentile"]
percentiles_to_check <- percentiles_to_check
sort(percentiles_to_check)
df <- cox_res_non_rand_chemo_Fadj[which(cox_res_non_rand_chemo_Fadj$cancer_type == cancer_type & cox_res_non_rand_chemo_Fadj$percentile %in% percentiles_to_check),]
df <- df[order(df$percentile),]
median(df$meta_pvalue_FDR_adjust,na.rm=TRUE)
median(df$exp_HR,na.rm=TRUE)

# Save
output_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/phenotypes/survival/cox_res_non_rand_strafified_by_chemotherapy.txt"
# write.table(cox_res_non_rand_chemo_Fadj, file = output_path, 
#                 sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
cox_res_non_rand_chemo_Fadj <- read.table(file = output_path, header = TRUE, sep = "\t")
# output_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/phenotypes/survival/cox_res_rand_strafified_by_chemotherapy.txt"
# write.table(cox_res_rand_chemo_Fadj, file = output_path, 
#                 sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
# cox_res_rand_chemo_Fadj <- read.table(file = output_path, header = TRUE, sep = "\t")

#############################################
######## TREATMENT --> IMMUNOTHERAPY ########
#############################################

cox_res_non_rand_immuno_all <- survival_analysis_by_cancer_stratified_by_treatment(TCGA_NMDeff = TCGA_NMDeff, randomization = "no", 
                                                                treatment_type = "immunotherapy", survival_var = "PFS")
outlier_rows <- which(cox_res_non_rand_immuno_all$CI_2.5 > 20 | cox_res_non_rand_immuno_all$CI_97.5 > 20 | 
                    cox_res_non_rand_immuno_all$exp_HR == Inf | cox_res_non_rand_immuno_all$exp_HR > 15 | cox_res_non_rand_immuno_all$exp_HR < 0.001)
cox_res_non_rand_immuno <- cox_res_non_rand_immuno_all[-outlier_rows,]
# Remove meta-pvalues that actually comes from only 1 method, that are already removed previously
rem_rows <- which(( is.na(cox_res_non_rand_immuno$exp_HR) & !is.na(cox_res_non_rand_immuno$meta_pvalue) ))
cox_res_non_rand_immuno <- cox_res_non_rand_immuno[-rem_rows,]

## Randomization ##
# cox_res_rand_immuno <- c()
# for (rand_iteration in 1:10) {
#   cox_res_rand_iteration <- survival_analysis_by_cancer_stratified_by_treatment(TCGA_NMDeff = TCGA_NMDeff, randomization = "yes", treatment_type = "immunotherapy")
#   cox_res_rand_iteration$rand_iteration <- rand_iteration
#   if (length(cox_res_rand) == 0 ) {
#     cox_res_rand_immuno <- cox_res_rand_iteration
#   } else {
#     cox_res_rand_immuno <- rbind(cox_res_rand_immuno,cox_res_rand_iteration)
#   }
# }

# Same as RVAS analysis
# Take significant hits for each cancer type at different % FDR thresholds, 
# then, calculate empirical FDR based on hits on randomization for each FDR threshold

cox_res_non_rand_immuno_Fadj <- cox_res_non_rand_immuno %>%
    group_by(cancer_type) %>%
    mutate(meta_pvalue_FDR_adjust = p.adjust(meta_pvalue, method = "fdr"))

# cox_res_rand_immuno_Fadj <- cox_res_rand_immuno %>%
#     group_by(cancer_type,rand_iteration) %>%
#     mutate(meta_pvalue_FDR_adjust = p.adjust(meta_pvalue, method = "fdr"))

cox_res_non_rand_immuno_Fadj[cox_res_non_rand_immuno_Fadj$NMD_method == "ASE_PTC_NMD_triggering_0.2","meta_pvalue_FDR_adjust"] <- NA
# cox_res_rand_immuno_Fadj[cox_res_rand_immuno_Fadj$NMD_method == "ASE_PTC_NMD_triggering_0.2","meta_pvalue_FDR_adjust"] <- NA
cox_res_non_rand_immuno_Fadj <- data.frame(cox_res_non_rand_immuno_Fadj[order(cox_res_non_rand_immuno_Fadj$meta_pvalue_FDR_adjust),])
data.frame(cox_res_non_rand_immuno_Fadj[which(cox_res_non_rand_immuno_Fadj$meta_pvalue_FDR_adjust < 1),])

cancer_type <- "pancancer"
filter <- which(cox_res_non_rand_immuno_Fadj$cancer_type == cancer_type & cox_res_non_rand_immuno_Fadj$meta_pvalue_FDR_adjust < 1)
percentiles_to_check <- unique(cox_res_non_rand_immuno_Fadj[filter,"percentile"])
percentiles_to_check <- percentiles_to_check
sort(percentiles_to_check)
df <- cox_res_non_rand_immuno_Fadj[which(cox_res_non_rand_immuno_Fadj$cancer_type == cancer_type & cox_res_non_rand_immuno_Fadj$percentile %in% percentiles_to_check),]
df <- df[order(df$percentile),]
median(df$meta_pvalue_FDR_adjust,na.rm=TRUE)
median(df$exp_HR,na.rm=TRUE)

cancer_type <- "SKCM"
filter <- which(cox_res_non_rand_immuno_Fadj$cancer_type == cancer_type & cox_res_non_rand_immuno_Fadj$meta_pvalue_FDR_adjust < 1)
percentiles_to_check <- unique(cox_res_non_rand_immuno_Fadj[filter,"percentile"])
percentiles_to_check <- percentiles_to_check
sort(percentiles_to_check)
df <- cox_res_non_rand_immuno_Fadj[which(cox_res_non_rand_immuno_Fadj$cancer_type == cancer_type & cox_res_non_rand_immuno_Fadj$percentile %in% percentiles_to_check),]
df <- df[order(df$percentile),]
median(df$meta_pvalue_FDR_adjust,na.rm=TRUE)
median(df$exp_HR,na.rm=TRUE)

# Save
output_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/phenotypes/survival/cox_res_non_rand_strafified_by_immunotherapy.txt"
# write.table(cox_res_non_rand_immuno_Fadj, file = output_path, 
                # sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
cox_res_non_rand_immuno_Fadj <- read.table(file = output_path, header = TRUE, sep = "\t")
# output_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/phenotypes/survival/cox_res_rand_strafified_by_immunotherapy.txt"
# write.table(cox_res_rand_immuno_Fadj, file = output_path, 
#                 sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
# cox_res_rand_immuno_Fadj <- read.table(file = output_path, header = TRUE, sep = "\t")



### Main Figures ###

sig_df <- data.frame(cox_res_non_rand_Fadj[which(cox_res_non_rand_Fadj$meta_pvalue_FDR_adjust < 0.05),])
percentile <- 20
# OS_surv_1 <- survival_analysis_plots(TCGA_df = TCGA_NMDeff, test_variable = "endogenous_NMD_Consensus", 
#             percentile = percentile, cancer_type = "SKCM", randomization = "no", TCGA_cancer_metadata_drug_therapy_all = TCGA_cancer_metadata_drug_therapy_all,
#             KM_curve = "yes", treatment_type = "no_treatment", survival_var = "OS")
PFS_surv_1 <- survival_analysis_plots(TCGA_df = TCGA_NMDeff, test_variable = "endogenous_NMD_Consensus", 
            percentile = percentile, cancer_type = "SKCM", randomization = "no", TCGA_cancer_metadata_drug_therapy_all = TCGA_cancer_metadata_drug_therapy_all,
            KM_curve = "yes", treatment_type = "no_treatment", survival_var = "PFS")
# OS_surv_2 <- survival_analysis_plots(TCGA_df = TCGA_NMDeff, test_variable = "endogenous_NMD_Consensus", 
#             percentile = percentile, cancer_type = "SKCM", randomization = "no", TCGA_cancer_metadata_drug_therapy_all = TCGA_cancer_metadata_drug_therapy_all,
#             KM_curve = "yes", treatment_type = "immunotherapy", survival_var = "OS")
PFS_surv_2 <- survival_analysis_plots(TCGA_df = TCGA_NMDeff, test_variable = "endogenous_NMD_Consensus", 
            percentile = percentile, cancer_type = "SKCM", randomization = "no", TCGA_cancer_metadata_drug_therapy_all = TCGA_cancer_metadata_drug_therapy_all,
            KM_curve = "yes", treatment_type = "immunotherapy", survival_var = "PFS")
# OS_surv_3 <- survival_analysis_plots(TCGA_df = TCGA_NMDeff, test_variable = "endogenous_NMD_Consensus", 
#             percentile = percentile, cancer_type = "SKCM", randomization = "no", TCGA_cancer_metadata_drug_therapy_all = TCGA_cancer_metadata_drug_therapy_all,
#             KM_curve = "yes", treatment_type = "chemotherapy", survival_var = "OS")
PFS_surv_3 <- survival_analysis_plots(TCGA_df = TCGA_NMDeff, test_variable = "endogenous_NMD_Consensus", 
            percentile = percentile, cancer_type = "SKCM", randomization = "no", TCGA_cancer_metadata_drug_therapy_all = TCGA_cancer_metadata_drug_therapy_all,
            KM_curve = "yes", treatment_type = "chemotherapy", survival_var = "PFS")
# OS_surv_4 <- survival_analysis_plots(TCGA_df = TCGA_NMDeff, test_variable = "endogenous_NMD_Consensus", 
#             percentile = percentile, cancer_type = "SKCM", randomization = "no", TCGA_cancer_metadata_drug_therapy_all = TCGA_cancer_metadata_drug_therapy_all,
#             KM_curve = "yes", treatment_type = "all_treatment", survival_var = "OS")
PFS_surv_4 <- survival_analysis_plots(TCGA_df = TCGA_NMDeff, test_variable = "endogenous_NMD_Consensus", 
            percentile = percentile, cancer_type = "SKCM", randomization = "no", TCGA_cancer_metadata_drug_therapy_all = TCGA_cancer_metadata_drug_therapy_all,
            KM_curve = "yes", treatment_type = "all_treatment", survival_var = "PFS")
# OS_surv_5 <- survival_analysis_plots(TCGA_df = TCGA_NMDeff, test_variable = "endogenous_NMD_Consensus", 
#             percentile = 35, cancer_type = "PRAD", randomization = "no", KM_curve = "yes")
# OS_surv_6 <- survival_analysis_plots(TCGA_df = TCGA_NMDeff, test_variable = "endogenous_NMD_Consensus", 
#             percentile = 15, cancer_type = "COAD", randomization = "no", KM_curve = "yes")

surv_objects <- list("SKCM_no_treatment_PFS" = PFS_surv_1, "SKCM_immunotherapy_PFS" = PFS_surv_2, 
                      "SKCM_chemotherapy_PFS" = PFS_surv_3, "SKCM_all_treatment_PFS" = PFS_surv_4)
saveRDS(surv_objects, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig6/Fig6A_C.RData")

### Supplementary Figures ###

percentile <- 20
OS_surv_1 <- survival_analysis_plots(TCGA_df = TCGA_NMDeff, test_variable = "endogenous_NMD_Consensus", 
            percentile = percentile, cancer_type = "SKCM", randomization = "no", TCGA_cancer_metadata_drug_therapy_all = TCGA_cancer_metadata_drug_therapy_all,
            KM_curve = "yes", treatment_type = "all", survival_var = "OS")
OS_surv_2 <- survival_analysis_plots(TCGA_df = TCGA_NMDeff, test_variable = "ASE_PTC_NMD_triggering_0.2", 
            percentile = percentile, cancer_type = "SKCM", randomization = "no", TCGA_cancer_metadata_drug_therapy_all = TCGA_cancer_metadata_drug_therapy_all,
            KM_curve = "yes", treatment_type = "all", survival_var = "OS")
OS_surv_3 <- survival_analysis_plots(TCGA_df = TCGA_NMDeff, test_variable = "endogenous_NMD_Consensus", 
            percentile = 35, cancer_type = "PRAD", randomization = "no", TCGA_cancer_metadata_drug_therapy_all = TCGA_cancer_metadata_drug_therapy_all,
            KM_curve = "yes", treatment_type = "all", survival_var = "OS")
OS_surv_4 <- survival_analysis_plots(TCGA_df = TCGA_NMDeff, test_variable = "ASE_PTC_NMD_triggering_0.2", 
            percentile = 35, cancer_type = "PRAD", randomization = "no", TCGA_cancer_metadata_drug_therapy_all = TCGA_cancer_metadata_drug_therapy_all,
            KM_curve = "yes", treatment_type = "all", survival_var = "OS")
OS_surv_5 <- survival_analysis_plots(TCGA_df = TCGA_NMDeff, test_variable = "endogenous_NMD_Consensus", 
            percentile = 30, cancer_type = "ESCA_ac", randomization = "no", TCGA_cancer_metadata_drug_therapy_all = TCGA_cancer_metadata_drug_therapy_all,
            KM_curve = "yes", treatment_type = "all", survival_var = "OS")
OS_surv_6 <- survival_analysis_plots(TCGA_df = TCGA_NMDeff, test_variable = "ASE_PTC_NMD_triggering_0.2", 
            percentile = 30, cancer_type = "ESCA_ac", randomization = "no", TCGA_cancer_metadata_drug_therapy_all = TCGA_cancer_metadata_drug_therapy_all,
            KM_curve = "yes", treatment_type = "all", survival_var = "OS")
OS_surv_7 <- survival_analysis_plots(TCGA_df = TCGA_NMDeff, test_variable = "endogenous_NMD_Consensus", 
            percentile = 10, cancer_type = "LUAD", randomization = "no", TCGA_cancer_metadata_drug_therapy_all = TCGA_cancer_metadata_drug_therapy_all,
            KM_curve = "yes", treatment_type = "all", survival_var = "OS")
OS_surv_8 <- survival_analysis_plots(TCGA_df = TCGA_NMDeff, test_variable = "ASE_PTC_NMD_triggering_0.2", 
            percentile = 10, cancer_type = "LUAD", randomization = "no", TCGA_cancer_metadata_drug_therapy_all = TCGA_cancer_metadata_drug_therapy_all,
            KM_curve = "yes", treatment_type = "all", survival_var = "OS")
surv_objects <- list("SKCM_all_OS_ETG" = OS_surv_1, "SKCM_all_OS_ASE" = OS_surv_2, 
                      "PRAD_all_OS_ETG" = OS_surv_3, "PRAD_all_OS_ASE" = OS_surv_4, 
                      "ESCA_ac_all_OS_ETG" = OS_surv_5, "ESCA_ac_all_OS_ASE" = OS_surv_6,
                      "LUAD_all_OS_ETG" = OS_surv_7, "LUAD_all_OS_ASE" = OS_surv_8)
saveRDS(surv_objects, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig20/SuppFig20A_H.RData")

percentile <- 50

PFS_surv_1 <- survival_analysis_plots(TCGA_df = TCGA_NMDeff, test_variable = "ASE_PTC_NMD_triggering_0.2", 
            percentile = percentile, cancer_type = "SKCM", randomization = "no", TCGA_cancer_metadata_drug_therapy_all = TCGA_cancer_metadata_drug_therapy_all,
            KM_curve = "yes", treatment_type = "no_treatment", survival_var = "PFS")
PFS_surv_2 <- survival_analysis_plots(TCGA_df = TCGA_NMDeff, test_variable = "ASE_PTC_NMD_triggering_0.2", 
            percentile = percentile, cancer_type = "SKCM", randomization = "no", TCGA_cancer_metadata_drug_therapy_all = TCGA_cancer_metadata_drug_therapy_all,
            KM_curve = "yes", treatment_type = "immunotherapy", survival_var = "PFS")
PFS_surv_3 <- survival_analysis_plots(TCGA_df = TCGA_NMDeff, test_variable = "ASE_PTC_NMD_triggering_0.2", 
            percentile = percentile, cancer_type = "SKCM", randomization = "no", TCGA_cancer_metadata_drug_therapy_all = TCGA_cancer_metadata_drug_therapy_all,
            KM_curve = "yes", treatment_type = "chemotherapy", survival_var = "PFS")
# PFS_surv_4 <- survival_analysis_plots(TCGA_df = TCGA_NMDeff, test_variable = "ASE_PTC_NMD_triggering_0.2", 
#             percentile = percentile, cancer_type = "SKCM", randomization = "no", TCGA_cancer_metadata_drug_therapy_all = TCGA_cancer_metadata_drug_therapy_all,
#             KM_curve = "yes", treatment_type = "all_treatment", survival_var = "PFS")

surv_objects <- list("SKCM_no_treatment_PFS_ASE" = PFS_surv_1, "SKCM_immunotherapy_PFS_ASE" = PFS_surv_2, 
                    "SKCM_chemotherapy_PFS_ASE" = PFS_surv_3)
                    # "SKCM_chemotherapy_PFS_ASE" = PFS_surv_3, "SKCM_all_treatment_PFS_ASE" = PFS_surv_4)
                    
saveRDS(surv_objects, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig21/SuppFig21A_C.RData")

### Supplementary Tables ###

# Create Supp Table S4 for PFS - immunotherapy
cox_res_non_rand_immuno_Fadj_SuppTableS4 <- cox_res_non_rand_immuno_Fadj %>%
    filter( (cancer_type == "SKCM" & percentile %in% seq(20,30,5) ) |
    ( cancer_type == "pancancer" & percentile %in% seq(20,25,5) ) ) %>%
    mutate(NMD_method = ifelse(NMD_method == "endogenous_NMD_Consensus","ETG","ASE"),
          outcome = "PFS",
          treatment = "immunotherapy") %>%
    select(cancer_type,percentile,NMD_method,outcome,treatment,exp_HR,CI_2.5,CI_97.5,p_value,meta_pvalue,meta_pvalue_FDR_adjust) %>%
    arrange(cancer_type,percentile,NMD_method)
colnames(cox_res_non_rand_immuno_Fadj_SuppTableS4) <- c("cancer_type","percentile","NMD_method","outcome","treatment","exp_HR",
                              "CI_low","CI_high","p_value","meta_p_value","FDR_meta_p_value")
#Save
output_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/Tables/SuppTableS4.txt"
write.table(cox_res_non_rand_immuno_Fadj_SuppTableS4, file = output_path, 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE) 

# Create Supp Table S5 for PFS - chemotherapy (FDR < 5%)
cox_res_non_rand_immuno_Fadj_SuppTableS5 <- cox_res_non_rand_chemo_Fadj %>%
    filter( (cancer_type == "PAAD" & percentile %in% seq(20,35,5) ) |
    ( cancer_type == "LGG" & percentile %in% seq(15,50,5) ) |
    ( cancer_type == "SKCM" & percentile %in% seq(40,50,5) ) |
    ( cancer_type == "BLCA_Basal_scc" & percentile %in% 50 ) |
    ( cancer_type == "ESCA_scc" & percentile %in% 50 ) ) %>%
    mutate(NMD_method = ifelse(NMD_method == "endogenous_NMD_Consensus","ETG","ASE"),
          outcome = "PFS",
          treatment = "chemotherapy") %>%
    select(cancer_type,percentile,NMD_method,outcome,treatment,exp_HR,CI_2.5,CI_97.5,p_value,meta_pvalue,meta_pvalue_FDR_adjust) %>%
    arrange(cancer_type,percentile,NMD_method)
colnames(cox_res_non_rand_immuno_Fadj_SuppTableS5) <- c("cancer_type","percentile","NMD_method","outcome","treatment","exp_HR",
                              "CI_low","CI_high","p_value","meta_p_value","FDR_meta_p_value")
#Save
output_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/Tables/SuppTableS5.txt"
write.table(cox_res_non_rand_immuno_Fadj_SuppTableS5, file = output_path, 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE) 

# Create Supp Table S6 for PFS - chemotherapy (FDR < 25%)
cox_res_non_rand_immuno_Fadj_SuppTableS6 <- cox_res_non_rand_chemo_Fadj %>%
    filter( (cancer_type == "PAAD" & percentile %in% seq(10,50,5) ) |
    ( cancer_type == "LGG" & percentile %in% seq(5,50,5) ) |
    ( cancer_type == "SKCM" & percentile %in% seq(40,50,5) ) |
    ( cancer_type == "BLCA_Basal_scc" & percentile %in% 50 ) |
    ( cancer_type == "BLCA_Lum_pap" & percentile %in% seq(25,50,5) ) |
    ( cancer_type == "CESC" & percentile %in% 10 ) |
    ( cancer_type == "STAD" & percentile %in% seq(45,50,5) ) |
    ( cancer_type == "STAD_MSI" & percentile %in% seq(45,50,5) ) |
    ( cancer_type == "ESCA_scc" & percentile %in% 50 ) ) %>%
    mutate(NMD_method = ifelse(NMD_method == "endogenous_NMD_Consensus","ETG","ASE"),
          outcome = "PFS",
          treatment = "chemotherapy") %>%
    select(cancer_type,percentile,NMD_method,outcome,treatment,exp_HR,CI_2.5,CI_97.5,p_value,meta_pvalue,meta_pvalue_FDR_adjust) %>%
    arrange(cancer_type,percentile,NMD_method)
colnames(cox_res_non_rand_immuno_Fadj_SuppTableS6) <- c("cancer_type","percentile","NMD_method","outcome","treatment","exp_HR",
                              "CI_low","CI_high","p_value","meta_p_value","FDR_meta_p_value")
#Save
output_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/Tables/SuppTableS6.txt"
write.table(cox_res_non_rand_immuno_Fadj_SuppTableS6, file = output_path, 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE) 





























##########################################################################################################33



# 5) Survival anaysis for CNA-PCs

for (percentile in c(5,10,15,20,25,30,35,40,45,50)) {
  for ( cancer in c(levels(TCGA_NMDeff$cancer_type),"pancancer") ) {
    print(paste0("----------------- CANCER TYPE -------------------> ",cancer))
    for (NMD_method in c("ASE_PTC_NMD_triggering_0.2","endogenous_NMD_Consensus")) {
      survival_analysis_plots(TCGA_NMDeff, test_variable = "ASE_PTC_NMD_triggering_0.2", percentile = percentile, cancer_type = cancer)
      cox_res_df <- survival_analysis_plots(TCGA_df = TCGA_NMDeff, test_variable = "endogenous_NMD_Consensus", percentile = percentile, cancer_type = cancer)
      cox_res_df$cancer_type <- cancer
      cox_res_df$percentile <- percentile
      cox_res_df$NMD_method <- NMD_method
      # survival_analysis_plots(TCGA_df, test_variable = "Dim.3", percentile = percentile, cancer_type = cancer)
      # survival_analysis_plots(TCGA_df, test_variable = "Dim.37", percentile = percentile, cancer_type = cancer)
      # survival_analysis_plots(TCGA_df, test_variable = "Dim.53", percentile = percentile, cancer_type = cancer)
      # survival_analysis_plots(TCGA_df, test_variable = "Dim.96", percentile = percentile, cancer_type = cancer)
      # survival_analysis_plots(TCGA_df, test_variable = "Dim.13", percentile = percentile, cancer_type = cancer)
    }
  }
}






