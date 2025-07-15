library("GWASTools")
library("VennDiagram")
library("dplyr")
library("ggplot2")

RGWAS_df <- function(dataset_name, database, NMD_method, geneList, randomGenes, randomization, geneVariantsThres, VAF, tissues) {

    if (NMD_method == "ASE") {
        NMD_method_VAF <- paste0(NMD_method,"_",VAF)
    } else {
        NMD_method_VAF <- NMD_method
    }
    if (geneList == "yes") {
        genelist <- "withgenelist"
        if (randomGenes == "yes") {
            randomGenes_info <- "random"
        } else if (randomGenes == "no") {
            randomGenes_info <- "selected"
        }
    } else if (geneList == "no") {
        genelist <- "nogenelist"
        randomGenes_info <- "all"
    }

    dir_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/rare_germline_variants/",database,"/",genelist,"/",dataset_name)
    plots_path <- paste0(dir_path,"/plots")
    RGWAS_res_path <-paste0(dir_path,"/results")

    if (! dir.exists(plots_path) ) {
        dir.create(plots_path)
    }

    RGWAS_df_tissue_all <- c()
    for (tissue in tissues) {
        print(tissue)
        # RGWAS results
        RGWAS_res_name_pattern <- paste0(NMD_method_VAF,"_",gsub("TCGA-","",tissue),".*_2_genes_",randomGenes_info,"_randomization_",randomization,".txt")
        RGWAS_res_name <- grep(RGWAS_res_name_pattern, list.files(RGWAS_res_path), value = TRUE)
        #print(RGWAS_res_name)
        tryCatch({
            error <<- FALSE
            RGWAS_res_tissue <- read.table( file = paste0(RGWAS_res_path,"/",RGWAS_res_name),
                            sep = "\t", header = TRUE)
        },error = function(e) {
            print("Tissue RGWAS results not found...")
            error <<- TRUE
        })
        if (isTRUE(error)) {next}
        RGWAS_res_tissue$ranking <- 1:nrow(RGWAS_res_tissue)
        # Variant-gene threshold
        RGWAS_res <- RGWAS_res_tissue[which(RGWAS_res_tissue$n_carriers >= geneVariantsThres),]
        #print(dim(RGWAS_res))
        if (nrow(RGWAS_res) == 0 ) {next}
        if (length(RGWAS_df_tissue_all) == 0) {
            RGWAS_df_tissue_all <- RGWAS_res
        } else {
            RGWAS_df_tissue_all <- rbind(RGWAS_df_tissue_all,RGWAS_res)
        }
    }
    return(RGWAS_df_tissue_all)
}

QQplots_and_pval_histograms <- function(data,genelist,plots_type) {
            
            lambdas <- data.frame(database = NA, NMD_method = NA, dataset = NA, tissue = NA)
            # Obtain name info
            tissue_name <- as.character(unique(data$tissue))
            random_geneset_name <- as.character(unique(data$random_geneset))
            if (random_geneset_name == "yes") {
                random_geneset_name_char <- "Selected"
            } else if (random_geneset_name == "no") {
                random_geneset_name_char <- "Random"
            } else if (random_geneset_name == "all") {
                random_geneset_name_char <- "All"
            }
            dataset_name <- as.character(unique(data$dataset))
            NMD_method_name <- as.character(unique(data$NMD_method))

            if (plots_type == "dataset") {
                plots_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/rare_germline_variants/GTEx/",genelist,"/",dataset_name,"/plots")
            } else if (plots_type == "all") {
                plots_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/rare_germline_variants/GTEx/",genelist,"/QC_plots")
                dataset_name <- "all"
            }

            # 1) QQplots
            # Calculate lambda
            p_values_randomization <- data[data$randomization == "yes","pValue"] %>% pull()
            p_values_non_randomization <- data[data$randomization == "no","pValue"]  %>% pull()
            ylim_max <- max(c(-log10(p_values_randomization),log10(p_values_non_randomization)))
            lambda_randomization <- lambda(p_values_randomization)
            lambda_non_randomization <- lambda(p_values_non_randomization)

            # 1.1) Non randomization
            plot_title <- paste0(random_geneset_name_char," Genes --- Non-Randomization --> Lambda = ",lambda_non_randomization)
            df_non_randomization <- data[data$randomization == "no",]
            df_non_randomization <- df_non_randomization[order(df_non_randomization[,"pValue"]),]

            # Expected p-values
            n1 <- length(p_values_non_randomization)
            if (n1 != 0) {
                a <- 1:n1
                pvals_exp <- a/n1
                df_non_randomization[,"pValue_expected"] <- pvals_exp
                # QQplot       
                p1 <- ggplot(df_non_randomization, aes(x = -log10(eval(parse(text="pValue_expected"))), y = -log10(eval(parse(text="pValue"))) )) +
                    geom_point(size = 2) + ylim(c(0,ylim_max)) +
                    geom_abline(intercept = 0, slope = 1) +
                    ylab("-log10(P-values") + xlab("-log10(expected P-values") + ggtitle(plot_title) +
                    theme_classic(base_size = 16)
            }

            # 1.2) Randomization
            plot_title <- paste0(random_geneset_name_char," Genes --- Randomization --> Lambda = ",lambda_randomization)
            df_randomization <- data[data$randomization == "yes",]
            df_randomization <- df_randomization[order(df_randomization[,"pValue"]),]

            # Expected p-values
            n2 <- length(p_values_randomization)
            if (n2 != 0) {
                a <- 1:n2
                pvals_exp <- a/n2
                df_randomization[,"pValue_expected"] <- pvals_exp
                # QQplot       
                p2 <- ggplot(df_randomization, aes(x = -log10(eval(parse(text="pValue_expected"))), y = -log10(eval(parse(text="pValue"))) )) +
                    geom_point(size = 2) + ylim(c(0,ylim_max)) +
                    geom_abline(intercept = 0, slope = 1) +
                    ylab("-log10(P-values") + xlab("-log10(expected P-values") + ggtitle(plot_title) +
                    theme_classic(base_size = 16)
            }
            if ( (n1 != 0) & (n2 != 0) ) {
                # Final qqplot
                qqplot <- paste0(plots_path,"/",tissue_name,"_",NMD_method_name,"_qqplot_genes_",
                                    random_geneset_name_char,"_",dataset_name,".png")
                p <- cowplot::plot_grid(plotlist=list(p1,p2), labels = "AUTO", align = "v", ncol = 2, nrow = 1)
                ggsave(filename = qqplot, device = "png", plot = p, width = 5000, height = 2500, dpi = 300, units = "px")
            }
            # 2) Pvals histogram
            # 2.1) Non Randomization
            p1 <- ggplot(df_non_randomization, aes(x = eval(parse(text="pValue")))) +
                geom_histogram(aes(y=..density..), colour="black", fill="white", bins = 15) +
                geom_vline(aes(xintercept=median(eval(parse(text="pValue")), na.rm = TRUE)),color="red", linetype="dashed", size=1) + 
                xlab("p-values") + ylab("frequency") + ggtitle(paste0(random_geneset_name_char," Genes --- Non-Randomization")) +
                theme_classic(base_size = 16)
            # 2.2) Randomization
            p2 <- ggplot(df_randomization, aes(x = eval(parse(text="pValue")))) +
                geom_histogram(aes(y=..density..), colour="black", fill="white", bins = 15) +
                geom_vline(aes(xintercept=median(eval(parse(text="pValue")), na.rm = TRUE)),color="red", linetype="dashed", size=1) + 
                xlab("p-values") + ylab("frequency") +ggtitle(paste0(random_geneset_name_char," Genes --- Randomization")) +
                theme_classic(base_size = 16)
            # Final histogram
            histogram <- paste0(plots_path,"/",tissue_name,"_",NMD_method_name,"_histogram_genes_",
                                random_geneset_name_char,"_",dataset_name,".png")
            p <- cowplot::plot_grid(plotlist=list(p1,p2), labels = "AUTO", align = "v", ncol = 2, nrow = 1)
            ggsave(filename = histogram, device = "png", plot = p, width = 5000, height = 2500, dpi = 300, units = "px")
            # 3) Save lambda
            lambdas$dataset <- dataset_name
            lambdas$NMD_method <- NMD_method_name
            lambdas$lambda_non_randomization <- as.numeric(lambda_non_randomization)
            lambdas$lambda_randomization <- as.numeric(lambda_randomization)
            lambdas$tissue <- tissue_name
            print(lambdas)
            return(lambdas)
}

lambda <- function(p_values) {
    # Lambda calculation
    chisq <- qchisq(1-p_values,1)
    lambda = round(median(chisq,na.rm=TRUE)/qchisq(0.5,1),2)
    return(lambda)
}

GTEx_tissue_names_path <- "/g/strcombio/fsupek_cancer1/gpalou/GTEx/RNAseq/GTEx_tissues_clean.txt"
GTEx_tissues <- read.table(file = GTEx_tissue_names_path, stringsAsFactors = FALSE)$V1
NMD_genes <- read.table( file = "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/NMD_genes.txt",
                        sep = "\t", header = TRUE)
GTEx_tissues <- c(GTEx_tissues,"pantissue")
datasets <- c("PTV_Missense_CADD15_0.1perc","PTV_Missense_CADD25_0.1perc","PTV_0.1perc")

##################### 1) WITH GENE LIST #####################

# 1) Obtain datasets

RGWAS_replicated_hits <- list()
# GTEx
GTEx_RGWAS_df_all <- c()
for (randomization in c("yes","no")) {
    for (randomGenes in c("yes","no")) {
        for (dataset in datasets) {
            print(dataset)
            RWGAS_ASE <- RGWAS_df(dataset_name = dataset, NMD_method = "ASE", geneList = "yes", randomGenes = randomGenes,
                        randomization = randomization, geneVariantsThres = 2, VAF = 0.2, tissues = GTEx_tissues, database = "GTEx")
            RWGAS_ASE$NMD_method <- "ASE"
            RWGAS_end <- RGWAS_df(dataset_name = dataset, NMD_method = "endogenous", geneList = "yes", randomGenes = randomGenes,
                            randomization = randomization, geneVariantsThres = 2, VAF = 0.2, tissues = GTEx_tissues, database = "GTEx")
            RWGAS_end$NMD_method <- "Endogenous"
            if (is.null(nrow(RWGAS_end))) {
                RWGAS_dataset <- RWGAS_ASE
            } else {
                RWGAS_dataset <- rbind(RWGAS_ASE,RWGAS_end)
            }  
            RWGAS_dataset$dataset <- dataset
            RWGAS_dataset$random_geneset <- randomGenes
            RWGAS_dataset$randomization <- randomization
            if (length(GTEx_RGWAS_df_all) == 0) {
                GTEx_RGWAS_df_all <- RWGAS_dataset
            } else {
                GTEx_RGWAS_df_all <- rbind(GTEx_RGWAS_df_all,RWGAS_dataset)
            }
        }
    }
}
print(dim(GTEx_RGWAS_df_all))
table(duplicated(GTEx_RGWAS_df_all))
GTEx_RGWAS_df_all$database <- "GTEx"

# 2) QQ-plots
# For each NMD_method - tissue - dataset - Randomization vs non-randomization

GTEx_RGWAS_df_all %>%
        group_by(NMD_method,dataset,tissue,random_geneset) %>%
            do( QQplots_and_pval_histograms(.,"withgenelist") )

##################### 2) NO GENE LIST #####################

# 1) Obtain datasets

# GTEx
GTEx_RGWAS_df_all <- c()
for (randomization in c("yes","no")) {
    for (dataset in datasets) {
        print(dataset)
        RWGAS_ASE <- RGWAS_df(dataset_name = dataset, NMD_method = "ASE", geneList = "no", randomGenes = "all",
                    randomization = randomization, geneVariantsThres = 2, VAF = 0.2, tissues = GTEx_tissues, database = "GTEx")
        RWGAS_ASE$NMD_method <- "ASE"
        RWGAS_end <- RGWAS_df(dataset_name = dataset, NMD_method = "endogenous", geneList = "no", randomGenes = "all",
                        randomization = randomization, geneVariantsThres = 2, VAF = 0.2, tissues = GTEx_tissues, database = "GTEx")
        RWGAS_end$NMD_method <- "Endogenous"
        if (is.null(nrow(RWGAS_end))) {
            RWGAS_dataset <- RWGAS_ASE
        } else {
            RWGAS_dataset <- rbind(RWGAS_ASE,RWGAS_end)
        }  
        RWGAS_dataset$dataset <- dataset
        RWGAS_dataset$randomization <- randomization
        RWGAS_dataset$random_geneset <- "all"
        if (length(GTEx_RGWAS_df_all) == 0) {
            GTEx_RGWAS_df_all <- RWGAS_dataset
        } else {
            GTEx_RGWAS_df_all <- rbind(GTEx_RGWAS_df_all,RWGAS_dataset)
        }
    }
    
}
print(dim(GTEx_RGWAS_df_all))
GTEx_RGWAS_df_all$database <- "GTEx"

# 2) QQ-plots
# For each NMD_method - tissue - dataset - Randomization vs non-randomization
all_lambdas <- GTEx_RGWAS_df_all %>%
        group_by(NMD_method,dataset,tissue) %>%
            do( QQplots_and_pval_histograms(.,"nogenelist","dataset") )
all_lambdas$database <- "GTEx"
genelist <- "nogenelist"
output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/rare_germline_variants/GTEx/",genelist,"/GTEx_",genelist,"_all_lambdas.txt")
write.table(data.frame(all_lambdas), file = output_path, sep = "\t", quote = FALSE,
                      col.names = TRUE, row.names = FALSE)
