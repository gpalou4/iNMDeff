library("ggplot2")

# find_best_PC_old <- function(all_lambdas_filt,lambda_diff) {
#     best_PC <- NULL
#     for (PC in 1:nrow(all_lambdas_filt)) {
#         if ( is.null(best_PC) ) {
#             # If lambda is below 1 we will not correct
#             if (sum(all_lambdas_filt[PC,] < 1) == 2) {
#                 isBestPC <- 2
#             } else if (sum(all_lambdas_filt[PC,] < 1) == 1) {
#                 isBestPC <- 1
#             } else {
#                 # Check if both mutations have lambda 0.8 <= lambda >= 1.20, otherwise take the next PC
#                 isBestPC <- sum(all_lambdas_filt[PC,] >= (1-lambda_diff) | all_lambdas_filt[PC,] <= (1+lambda_diff))
#             }
#             if (ncol(all_lambdas_filt) > 1 ) { 
#                 if ( ( (mut == "CNV") | (mut == "miss_tr") ) && ( isBestPC == 2 ) ) {
#                     best_PC <- PC
#                     print(all_lambdas_filt[PC,])
#                 }
#             } else if (ncol(all_lambdas_filt) == 1 ) {
#                 if (isBestPC == 1) {
#                     best_PC <- PC
#                     print(all_lambdas_filt[PC,])
#                 }
#             }
#         } 
#     }
#     return(best_PC)
# }

find_best_PC <- function(all_lambdas_filt,lambda_diff) {
    best_PC <- NULL
    col_NA <- which( colSums(is.na(all_lambdas_filt)) > round(nrow(all_lambdas_filt)/2) )
    if (length(col_NA) != 0 ) {
        all_lambdas_filt <- all_lambdas_filt[,-col_NA, drop = FALSE]
    }
    for (PC in 1:nrow(all_lambdas_filt)) {
        all_lambdas_filt_PC <- all_lambdas_filt[PC,]
        if ( is.null(best_PC) ) {
            # Check if at least 1 mutation has lambda > 1.25
            to_correct <- sum(all_lambdas_filt_PC >= (1+lambda_diff) )
            if (to_correct == 0) {
                best_PC <- PC
                print(all_lambdas_filt_PC)
            }
        } 
    }
    return(best_PC)
}

TCGA_cancer_names_path <- "/g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/TCGA_projects_names.txt"
#TCGA_cancer_names_path <- "/g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/TCGA_cancers_strat.txt"
TCGA_cancers <- read.table(file = TCGA_cancer_names_path, stringsAsFactors = FALSE)$V1
TCGA_cancers <- TCGA_cancers[!TCGA_cancers%in%c("TCGA-LAML","TCGA-MESO","TCGA-SKCM")]

euclidean <- function(a, b) sqrt(sum((a - b)^2))
best_PC_cancers <- data.frame(TCGA_cancers = c(TCGA_cancers,"pancancer"), best_PC_som_CNV = NA, best_PC_som_mut_miss_tr = NA, best_PC_som_mut_synonymous = NA,
                                lambda_som_CNV_amp = NA, lambda_som_CNV_del = NA, lambda_som_miss = NA, lambda_som_tr = NA, lambda_som_syn = NA)

NMD_method <- "endogenous" # ASE / endogenous / PTC
if (NMD_method == "ASE") {
    NMD_method_VAF <- paste0(NMD_method,"_",VAF)
} else {
    NMD_method_VAF <- NMD_method
}
VAF <- 0.2
#TCGA_cancers <- "pancancer"
#TCGA_cancer <- "pancancer"

test <- "test_1"
robust_SPCA <- "no"
alpha <- "1e-04"
files_left <- list()

for (TCGA_cancer in TCGA_cancers) {

    print(TCGA_cancer)
    TCGA_cancer_original <- strsplit(TCGA_cancer, "_")[[1]][1]
    
    # 1) Number of PCs
    if (TCGA_cancer == "pancancer") {
        if (NMD_method == "ASE") {
            PCs <- 100
        } else {
            PCs <- 500
        }
        
    } else {
        TCGA_CNV_PCA <- read.table(file = paste0("/g/strcombio/fsupek_cancer1/gpalou/TCGA_CNV_PCA/",TCGA_cancer_original,"_sparse_PCA_ind_",alpha,"_robust_",robust_SPCA,".txt"),
                                            header = TRUE, stringsAsFactors = FALSE, sep = "\t")
        # Number of PCs
        PCs <- ncol(TCGA_CNV_PCA)
    }
    # 2) Obtain lambdas per PC
    all_lambdas <- c()
    PCs_left <- c()
    for (PC in 1:PCs) {
        lambda_path <- paste0("/home/gpalou/analysis_results/NMD_project/associations/somatic_associations/optimizing_PCs/lambda_values/",TCGA_cancer,"_",NMD_method_VAF,"_CGC_somatic_mut_CNV_PCs_",PC,"_lambda_",test,"_",alpha,"_robust_",robust_SPCA,".txt")
        if (file.exists(lambda_path)) {
            lambda <- read.table(lambda_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
        } else {
            PCs_left <- c(PCs_left,PC)
            all_lambdas <- rbind(all_lambdas,rep(NA,5))
            next
        }
        if ( is.null(all_lambdas) ) {
            all_lambdas <- lambda
        } else {
            all_lambdas <- rbind(all_lambdas,lambda)
        } 
    }
    files_left[[TCGA_cancer]] <- PCs_left

    print(head(all_lambdas))
    col_NA <- which( colSums(is.na(all_lambdas)) > round(nrow(all_lambdas)/2) )
    if ( length(col_NA) == 5 ) {
        print("Error, Lambdas NA or no files")
        next
    }
    
    for (mut in c("CNV","miss_tr","syn")) {
        print(mut)
        if (mut == "CNV") {
            all_lambdas_filt <- all_lambdas[,colnames(all_lambdas) %in% c("CNV_amp","CNV_del"), drop = FALSE]
            col <- "best_PC_som_CNV"
            col2 <- c("lambda_som_CNV_amp","lambda_som_CNV_del")
        } else if (mut == "miss_tr") {
            all_lambdas_filt <- all_lambdas[,colnames(all_lambdas) %in% c("mut_truncating","mut_missense"), drop = FALSE]
            col <- "best_PC_som_mut_miss_tr"
            col2 <- c("lambda_som_miss","lambda_som_tr")
        } else if (mut == "syn") {
            all_lambdas_filt <- all_lambdas[,colnames(all_lambdas) %in% c("mut_synonymous"), drop = FALSE]
            col <- "best_PC_som_mut_synonymous"
            col2 <- "lambda_som_syn"
        }
        best_PC <- NULL
        lambda_diff <- 0
        best_PC <- find_best_PC(all_lambdas_filt, lambda_diff = 0.25)
        # If don't find a best PC increase the lambda difference slightly until 1.5
        if (is.null(best_PC)) {
            while ( is.null(best_PC) ) {
                lambda_diff <- lambda_diff + 0.05
                best_PC <- find_best_PC(all_lambdas_filt, lambda_diff = 0.25 + lambda_diff)
            }  
        }
        if (!is.null(best_PC) & lambda_diff <= 0.5) {
            best_PC_cancers[best_PC_cancers$TCGA_cancer%in%TCGA_cancer,col] <- best_PC
            best_PC_cancers[best_PC_cancers$TCGA_cancer%in%TCGA_cancer,col2] <- as.numeric(all_lambdas_filt[best_PC,, drop = FALSE])
        }
    }

    # 3) Lambda plot
    if (length(col_NA) != 0 ) {
        all_lambdas <- all_lambdas[,-col_NA]
    }
    all_lambdas$e_dist <- NULL
    all_lambdas_df <- stack(all_lambdas)
    n <- 5 - length(col_NA)
    all_lambdas_df$PC <- rep(1:nrow(all_lambdas),n)

    lambda_plot <- paste0("/home/gpalou/analysis_results/NMD_project/associations/somatic_associations/optimizing_PCs/plots_lambda/",TCGA_cancer,"_lambda.png")
    png(lambda_plot, width = 3500, height = 2500, res = 300)
    p <- ggplot(all_lambdas_df, aes(y = values,x = PC, color = ind)) +
        geom_point(size = 4, alpha = 1) +
        geom_line() +
        geom_hline(yintercept = 1, size = 2) +
        theme_classic() + ylim(c(0,4)) +
        xlab("PCs") + ylab("Lambda") + ggtitle(paste0(TCGA_cancer," --> Lambda over sparse-PCs")) +
        theme(legend.position="top", plot.title = element_text(hjust = 0.5, size = 18),
            axis.text.x = element_text(color="black", size=18),
            axis.text.y = element_text(color="black", size=18),
            legend.title = element_text(colour="black", size=18, face="bold"))
    print(p)
    dev.off()

}

# Put manually the best PC for pancancer
if (NMD_method == "ASE") {
    # best_PC_cancers[best_PC_cancers$TCGA_cancers %in% "pancancer","best_PC_som_CNV"] <- 18
    # best_PC_cancers[best_PC_cancers$TCGA_cancers %in% "pancancer","best_PC_som_mut_miss_tr"] <- 1
    # best_PC_cancers[best_PC_cancers$TCGA_cancers %in% "pancancer","best_PC_som_mut_synonymous"] <- 1
} else if (NMD_method == "endogenous") {
    best_PC_cancers[best_PC_cancers$TCGA_cancers %in% "pancancer","best_PC_som_CNV"] <- 300
    best_PC_cancers[best_PC_cancers$TCGA_cancers %in% "pancancer","best_PC_som_mut_miss_tr"] <- 1
    best_PC_cancers[best_PC_cancers$TCGA_cancers %in% "pancancer","best_PC_som_mut_synonymous"] <- 1
}

write.table(best_PC_cancers, file = paste0("/home/gpalou/analysis_results/NMD_project/associations/somatic_associations/optimizing_PCs/",NMD_method_VAF,"_best_PC_all_TCGA_cancers.txt"), 
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
