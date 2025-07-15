
################################################################################################

########################################## LIBRARIES ###########################################

################################################################################################

# conda activate RGWAS

library(tidyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(tclust)
library(metap)
library(MASS)
library(data.table)
library(boot)
library(broom)
library(ggrepel)
library(ggsignif)
library(ggpubr)
library(SKAT)
library(fastDummies)
library("FactoMineR")
library("factoextra")
library("broom.mixed")

################################################################################################

########################################## SCRIPT ##############################################

################################################################################################

args <- commandArgs(trailingOnly=TRUE)
tissue_type <- args[1]
rare_germline_variants_name <- args[2]
NMD_method <- args[3]
VAF <- args[4]
geneList <- args[5]
numberGenesThres <- as.numeric(args[6])
randomGenes <- args[7]
randomization <- args[8]
print("ARGUMENS --> ")
print(args)

####################################
# tissue_type <- "Adipose_Visceral_Omentum"
# tissue_type <- "pantissue"
# rare_germline_variants_name <- "PTV_Missense_CADD25_0.1perc"
# NMD_method <- "ASE"
# VAF <- "0.2"
# geneList <- "yes"
# numberGenesThres <- 2
# randomGenes <- "yes"
# randomization <- "no"
####################################

# Decide whether use gene list or no
if (geneList == "yes") {
    genelist <- TRUE
} else {
    genelist <- FALSE
}
print(genelist)

# 1) Input files

# 1.1) Phenotype
# 1.1) Phenotype - NMD efficiency

if (NMD_method == "ASE") {
    phenotype <- paste0(NMD_method,"_stopgain_",VAF)
    NMD_method_VAF <- paste0(NMD_method,"_",VAF)
} else if (NMD_method == "endogenous") {
    phenotype <- paste0(NMD_method,"_NMD_global_2_shared")
    NMD_method_VAF <- NMD_method
} else if (NMD_method == "PTC") {
    phenotype <- paste0(NMD_method,"_stopgain_NMD_triggering")
}

# GTEx NMDeff
sample_NMD_efficiencies_GTEx_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/GTEx/pantissue/NMD_efficiencies_GTEx.txt"
sample_NMD_efficiencies_GTEx <- read.table(file = sample_NMD_efficiencies_GTEx_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# Convert some columns to factors
factor_cols <- c("tissue","sample")
sample_NMD_efficiencies_GTEx[factor_cols] <- lapply(sample_NMD_efficiencies_GTEx[factor_cols], factor) 
# Filters for GTEx
sample_NMD_efficiencies_GTEx[which(sample_NMD_efficiencies_GTEx$ASE_num_PTCs_0.2 < 3),c("ASE_stopgain_0.2")] <- NA
sample_NMD_efficiencies_GTEx[which(sample_NMD_efficiencies_GTEx$ASE_num_PTCs_0.01 < 3),c("ASE_stopgain_0.01")] <- NA
input_phenotypes <- sample_NMD_efficiencies_GTEx
# Change some columns
filter <- colnames(input_phenotypes) %in% c("sample","sex")
colnames(input_phenotypes)[filter] <- c("sample_short","sex")

# 1.2) Rare germline variants
input_path <- paste0("/g/strcombio/fsupek_cancer1/gpalou/GTEx/WES_VCF/rare_germline_variants/GTEx_germline_input_variants_",rare_germline_variants_name,".txt")
variants_gnomad_allinfo <- read.csv(file = input_path, head=T,sep ="\t",stringsAsFactors = F)

# Change columns
variants_gnomad_allinfo$Gene.refGene <- variants_gnomad_allinfo$gene_name
filter <- colnames(variants_gnomad_allinfo) %in% c("GTEx_sample")
colnames(variants_gnomad_allinfo)[filter] <- c("sample_short")

# 1.4) PCA from common variants and other covariates
GTEx_tissues <- as.character(read.table(file = "/g/strcombio/fsupek_cancer1/gpalou/GTEx/RNAseq/GTEx_tissues_clean.txt")$V1)
PCA_commonVariants <- c()
for (GTEx_tissue in GTEx_tissues) {
    if (GTEx_tissue == "Brain_Spinal_cord_cervical_c1") {
        GTEx_tissue <- "Brain_Spinal_cord_cervical_c-1"
    } else if (GTEx_tissue == "Cells_EBVtransformed_lymphocytes") {
        GTEx_tissue <- "Cells_EBV-transformed_lymphocytes"
    }
    print(GTEx_tissue)
    input_data <- paste0("/g/strcombio/fsupek_cancer1/gpalou/GTEx/samples_metadata/GTEx_Analysis_v8_eQTL_covariates/",GTEx_tissue,".v8.covariates.txt")
    tryCatch({
        error <- FALSE
        tissue_covs <- read.table( file = input_data, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
        rownames(tissue_covs) <- tissue_covs$ID 
        tissue_covs$ID <- NULL
        tissue_covs <- data.frame(t(tissue_covs))
        tissue_covs$sample_short <- gsub("\\.","-",rownames(tissue_covs))
        tissue_covs$tissue <- GTEx_tissue
        tissue_covs <- tissue_covs[,-grep("InferredCov.*",colnames(tissue_covs))]
    },error = function(e) {
        print("Tissue with no covariates file...")
        error <- TRUE
    })
    if (length(PCA_commonVariants) == 0) {
        PCA_commonVariants <- tissue_covs
    } else {
        PCA_commonVariants <- rbind(PCA_commonVariants, tissue_covs)
    }
}
# Only 6  tissues do not have covariates (we will not use them)
colnames(PCA_commonVariants)[1:5] <- paste0("Dim.",1:5)
rownames(PCA_commonVariants) <- NULL
PCA_commonVariants <- PCA_commonVariants[!duplicated(PCA_commonVariants),]
print(head(PCA_commonVariants))

# 1.5) Selected genes to test
# SelectedGenes <- read.csv(file = '/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/NMD_genes.txt',head=TRUE,sep ="\t")
if (isTRUE(genelist) & randomGenes == "yes") {
    SelectedGenes <- read.table(file ="/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/random_genes_not_cancer_not_NMD_related_expanded.txt",
                                header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    selected_genes_info <- "random"
} else if (isTRUE(genelist) & randomGenes == "no") {
    SelectedGenes <- read.table(file = "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/NMD_genes_STRING_expanded.txt",
                            header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    selected_genes_info <- "selected"
} else if (!isTRUE(genelist)) {
    selected_genes_info <- "all"
}

# 1.6) Admixture samples (R Gay et al Science 2020)
# admixture_samples_path <- "/g/strcombio/fsupek_cancer1/gpalou/GTEx/samples_metadata/admixtured_samples_R_Gay_et_al_Science.txt"
# admixture_samples <- read.table(file = admixture_samples_path, header = FALSE, stringsAsFactors = FALSE, sep = "\t")$V1

# 1.7) Filter non-European samples and missing PCs samples
# Remove samples with no PCs (we only have 838 samples with PCs, those with WGS data)
# Remove non-European samples
# 1.7.1) PCA
# PCA_commonVariants_filt <- PCA_commonVariants[!duplicated(PCA_commonVariants$sample_short),]
# rownames(PCA_commonVariants_filt) <- paste0(PCA_commonVariants_filt$sample_short)
# PCA_commonVariants_filt$color <- "rest"
# PCA_commonVariants_filt[PCA_commonVariants_filt$sample_short %in% admixture_samples,"color"] <- "admixture"
# res_pca <- PCA(PCA_commonVariants_filt[,1:5], 
#             scale.unit = FALSE, ncp = 1000, graph = FALSE)
# path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/rare_germline_variants/GTEx/QC"
# png(paste0(path,"/PCA_genotype_PCs_ancestry.png"), width = 2500, height = 2000, res = 300)
# p <- fviz_pca_ind(res_pca, repel = FALSE, geom.ind = "point",
#                     axes = c(1,2),
#                     alpha.ind = 1,
#                     pointshape = 21, 
#                     fill = PCA_commonVariants_filt$color,
#                     legend.title = list(fill = "PCA")) + theme_minimal()
# print(p)
# dev.off()
# pca_df <- data.frame(res_pca$ind$coord)
# # African-american + admixture (103)
# outlier_rows_1 <- which(pca_df$Dim.1 < -0.03)
# # Asian American (12)
# outlier_rows_2 <- which(pca_df$Dim.2 < -0.1)
# # Admixture (117)
# #admixture_samples
# # All
# nonEuropean_samples <- unique(PCA_commonVariants_filt[unique(c(outlier_rows_1,outlier_rows_2)),"sample_short"])
# table(admixture_samples %in% nonEuropean_samples)
# # Samples with no PCs
# samples_with_missing_PCs <- as.character(unique(input_phenotypes[!input_phenotypes$sample_short %in% PCA_commonVariants$sample_short,"sample_short"]))
# all_outlier_samples <- unique(c(admixture_samples,nonEuropean_samples,samples_with_missing_PCs))

# table(PCA_commonVariants$sample_short %in% input_phenotypes$sample_short)

outlier_path <- paste0("/g/strcombio/fsupek_cancer1/gpalou/GTEx/samples_metadata/ancestry_and_admixture_and_PCA_outliers.txt")
# write.table(all_outlier_samples, file = outlier_path,
#         quote=FALSE, sep='\t',row.names=FALSE,col.names = F)
all_outlier_samples <- read.table(file = outlier_path)$V1

# 1.7.2) # Remove samples
input_phenotypes <- input_phenotypes %>%
                    filter(!sample_short %in% all_outlier_samples)
PCA_commonVariants <- PCA_commonVariants %>%
                    filter(!sample_short %in% all_outlier_samples)

# 2) Create dataset
variant_threshold <- 0.1
print(tissue_type)

# Take randomized phenotype?
if ( (randomization == "yes") ) {
    phenotype <- paste0(phenotype,"_randomized")
} 

####create template table with samples which should have a value
if (tissue_type == 'COAD_READ') {
    table_pheno <-
        input_phenotypes %>%
        dplyr::filter(tissue == 'COAD' | tissue == 'READ') %>%
        dplyr::select(sample_short,phenotype,tissue) %>%
        dplyr::filter(!is.na(!!rlang::sym(phenotype)) ) 
} else if (tissue_type == 'ESCA_EAC') {
    table_pheno <-
        input_phenotypes %>%
        dplyr::filter(tissue == 'ESCA') %>%
        mutate(sample_check=if_else(sample_short %in% ESCA_EAC$barcode, 'yes', 'no'
        )
        ) %>%
        dplyr::filter(sample_check == 'yes') %>%
        dplyr::select(-sample_check) %>%
        dplyr::select(sample_short,phenotype,tissue) %>%
        dplyr::filter(!is.na(!!rlang::sym(phenotype)) ) 
} else if(tissue_type == 'STAD_ESCA_EAC') {
    table_pheno <-
        input_phenotypes %>%
        dplyr::filter(tissue == 'STAD' | tissue== 'ESCA')%>%
        mutate(sample_check=if_else(sample_short %in% ESCA_ESCC$barcode, 'yes', 'no' ##take the samples which are not ESCC
        )
        ) %>%
        dplyr::filter(sample_check == 'no') %>%
        dplyr::select(-sample_check)  %>%
        dplyr::select(sample_short,phenotype,tissue) %>%
        dplyr::filter(!is.na(!!rlang::sym(phenotype)) )
} else if (tissue_type == 'pantissue') {
    table_pheno <-
        input_phenotypes %>%
        dplyr::select(sample_short,phenotype,tissue, age) %>%
        dplyr::filter(!is.na(!!rlang::sym(phenotype)) )
} else if (tissue_type == 'KIRC_KIRP') {
    table_pheno <-
        input_phenotypes %>%
        dplyr::filter(tissue == 'KIRC'|tissue=='KIRP') %>%
        dplyr::select(sample_short,phenotype,tissue) %>%
        dplyr::filter(!is.na(!!rlang::sym(phenotype)) )
} else {
    table_pheno <-
        input_phenotypes %>%
        dplyr::filter(tissue %in% tissue_type) %>%
        dplyr::select(sample_short,phenotype,tissue,age) %>%
        dplyr::filter(!is.na(!!rlang::sym(phenotype)) )
}

if ( !isTRUE(genelist) ) {
    genelist_info <- 'nogenelist'
} else {
    genelist_info <- 'withgenelist'
}

####make sure that at least 5% of the samples have a value for the phenotype!
# pheno_activity <- nrow(table_pheno[table_pheno[[phenotype]] != 0,])/nrow(table_pheno)
# if(pheno_activity > 0.0 & nrow(table_pheno) >= 50){

# Add PCA covs
matrix_analysis <- table_pheno
matrix_analysis <- merge(matrix_analysis, PCA_commonVariants, by = c("sample_short","tissue"), all.x = TRUE)
matrix_analysis <- matrix_analysis[!duplicated(matrix_analysis),]
matrix_analysis <- na.omit(matrix_analysis)
# Change covs
matrix_analysis$sex <- ifelse(matrix_analysis$sex == 1, "male", "female")
print(paste(nrow(matrix_analysis), ' samples',sep=''))

#create new output directory if not already existing
output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/rare_germline_variants/GTEx/",genelist_info,"/",rare_germline_variants_name)
if(! file.exists(output_path)){ #create directory if it does not exist
    dir.create(output_path, recursive = TRUE)
}

# 2.2) Start with Variants
# Filter Variants and add some info
variants_gnomad_samples <- variants_gnomad_allinfo %>%
    dplyr::filter(sample_short %in% matrix_analysis$sample_short)  %>%
    mutate(gene_variant_id= paste(Gene.refGene,Chr,Start,End,Ref,Alt,sep='__')) ##combine gene name with unique thingy
# At least 2 variants-samples in a gene is found

if ( !isTRUE(genelist) ) {
    variants_genes_higher2 <- variants_gnomad_samples %>%
        mutate(samplecheck=if_else(sample_short %in% matrix_analysis$sample_short, 'yes','no')) %>%
        dplyr::filter(samplecheck=='yes') %>%
        dplyr::select(-samplecheck)  %>%
        group_by(Gene.refGene,sample_short) %>%
        summarise(RDGV_perSample=n()) %>%
        ungroup() %>%
        group_by(Gene.refGene) %>%
        summarise(numberSamples=n()) %>%
        ungroup() %>%
        dplyr::filter(numberSamples >= numberGenesThres) 
} else if (isTRUE(genelist) ) {
    variants_genes_higher2 <- variants_gnomad_samples %>%
        mutate(samplecheck=if_else(sample_short %in% matrix_analysis$sample_short, 'yes','no')) %>%
        dplyr::filter(samplecheck=='yes') %>%
        dplyr::select(-samplecheck)  %>%
        mutate(samplecheck=if_else(Gene.refGene %in% SelectedGenes$Genes, 'yes','no')) %>%
        dplyr::filter(samplecheck=='yes') %>%
        dplyr::select(-samplecheck)  %>%
        group_by(Gene.refGene,sample_short) %>%
        summarise(RDGV_perSample=n()) %>%
        ungroup() %>%
        group_by(Gene.refGene) %>%
        summarise(numberSamples=n()) %>%
        ungroup() %>%
        dplyr::filter(numberSamples >= numberGenesThres) 
}

print(paste0("Number of genes to test --> ",nrow(variants_genes_higher2)))

# At least one gene to test in the tissue
if ( nrow(variants_genes_higher2) > 1) {
      
    germline_variants_matrix <- variants_gnomad_samples %>%
        dplyr::filter(Gene.refGene %in% variants_genes_higher2$Gene.refGene)
    germline_variants_matrix <- germline_variants_matrix %>%
        mutate(geno = if_else(gt == '1/1' | gt == '1|1',2,1)
        ) %>% 
        dplyr::select(sample_short,geno,gene_variant_id) %>%
        distinct() %>% # remove any duplicates, should not exist though (checked)
        spread(gene_variant_id,geno,is.na(0)) %>%
        mutate(samplecheck = if_else(sample_short %in% matrix_analysis$sample_short, 'yes','no' ##only samples with pheno
        )
        ) %>%
        dplyr::filter(samplecheck=='yes') %>%
        dplyr::select(-samplecheck)

    # Create the final matrix for testing
    combined_variants_phenotype <- germline_variants_matrix %>%
        left_join(matrix_analysis,by=c('sample_short'))
    # Put all samples without any rare variants to 0 --> NAs to zero
    for ( columns in grep("\\_\\_",colnames(combined_variants_phenotype)) ) {
        combined_variants_phenotype[,columns][is.na(combined_variants_phenotype[,columns])] <- 0
    }
    NumberVariantColumns <- nrow(variants_genes_higher2)
    samples_n <- nrow(combined_variants_phenotype)

    # Create GTEx barcode variable
    combined_variants_phenotype$GTEx_barcode <- combined_variants_phenotype$sample_short
    combined_variants_phenotype[which(table(combined_variants_phenotype$sample_short) < 5),"GTEx_barcode"] <- "other"
    
    # 2.3) Rare variants association testing
    ############## start with testing ################################################################################################################
    # Super important to initialize this matrix every time from scratch
    ResVariantsPValues <- data.frame(variant=rep('sample',NumberVariantColumns),
                                    pValue=rep(NA,NumberVariantColumns),
                                    pheno=rep(phenotype,NumberVariantColumns),
                                    n_carriers=rep(NA,NumberVariantColumns),
                                    n_variants=rep(NA,NumberVariantColumns),
                                    n_samples=rep(samples_n,NumberVariantColumns),
                                    tissue=rep(tissue_type,NumberVariantColumns),
                                    snps=rep(rare_germline_variants_name,NumberVariantColumns),
                                    model=rep('additive',NumberVariantColumns),
                                    test=rep('SKAT-O',NumberVariantColumns),
                                    rho_SKAT_O=rep(NA,NumberVariantColumns),
                                    double_del_carriers=rep(NA,NumberVariantColumns),
                                    stringsAsFactors = FALSE)
    
    if ( samples_n > 2 ) { #only if at least 50 samples go on 
        # Loop through all genes
        for(i in 1:NumberVariantColumns){
            print(i)           
            gene_to_test <- variants_genes_higher2$Gene.refGene[i] #gene to test
            # Grep the column names of the variants inside that gene
            gene_variants_to_test <- colnames(combined_variants_phenotype)[grepl(paste('^',gene_to_test,'__',sep=''),colnames(combined_variants_phenotype))] 
            number_variants <- length(gene_variants_to_test)
            # Select columns needed
            if (tissue_type == "pantissue") {
                combined_variants_phenotype_bla <- combined_variants_phenotype %>%
                            dplyr::select(all_of(gene_variants_to_test),phenotype,sex,age,Dim.1,Dim.2,Dim.3,Dim.4,Dim.5,tissue)
            } else {
                combined_variants_phenotype_bla <- combined_variants_phenotype %>%
                            dplyr::select(all_of(gene_variants_to_test),phenotype,sex,age,Dim.1,Dim.2,Dim.3,Dim.4,Dim.5,tissue,GTEx_barcode)
            }
            
            # formula <- paste0("ASE_stopgain_0.2 ~ A2M__chr12__9079692__9079692__G__A + as.factor(sex) + as.factor(age) + Dim.1 + Dim.2 + Dim.3 + Dim.4 + Dim.5 + as.factor(tissue)
            #             + (1 | GTEx_barcode)")
             
            # # stick to a generalized linear mixed-effects model for the negative binomial family
            # y <- suppressWarnings(glmer(formula = ASE_stopgain_0.2 ~ A2M__chr12__9079692__9079692__G__A + as.factor(sex) + (1 | GTEx_barcode), 
            #                             data = combined_variants_phenotype_bla))


            # Test tissue type if more than 1 tissue type in it, same for sex
            if ( length(unique(combined_variants_phenotype_bla$tissue)) > 1 & nrow(combined_variants_phenotype_bla[combined_variants_phenotype_bla$sex=='female',]) > 2 & nrow(combined_variants_phenotype_bla[combined_variants_phenotype_bla$sex=='male',]) > 2 ) {
                # Need to make dummy variable for age, sex and tissue type (discrete variables)
                combined_variants_phenotype_bla_dummy <- fastDummies::dummy_cols(combined_variants_phenotype_bla, select_columns = c("tissue","sex", "age"), remove_first_dummy = TRUE) %>%
                    dplyr::select(-tissue,-sex,-age);
                # Make the null model first, regressing the covariates on the phenotype
                SKAT_result_null <- SKAT_Null_Model(as.matrix(combined_variants_phenotype_bla_dummy[,phenotype]) ~ as.matrix(combined_variants_phenotype_bla_dummy[,(number_variants+2):ncol(combined_variants_phenotype_bla_dummy)]), out_type="C")
                # Perform the SKATO test
                SKAT_result <- SKAT(as.matrix(combined_variants_phenotype_bla_dummy[,1:number_variants]), SKAT_result_null, method="SKATO")
            } else if ( length(unique(combined_variants_phenotype_bla$tissue)) > 1 ) {
                # Need to make dummy variable for tissue
                combined_variants_phenotype_bla_dummy <- fastDummies::dummy_cols(combined_variants_phenotype_bla, select_columns = c("tissue","age","GTEx_barcode"), remove_first_dummy = TRUE) %>%
                    dplyr::select(-tissue,-sex,-age,-GTEx_barcode);
                # Make the null model first, regressing the covariates on the phenotype
                SKAT_result_null <- SKAT_Null_Model(as.matrix(combined_variants_phenotype_bla_dummy[,phenotype]) ~ as.matrix(combined_variants_phenotype_bla_dummy[,(number_variants+2):ncol(combined_variants_phenotype_bla_dummy)]), out_type="C")
                # Perform the SKATO test
                SKAT_result <- SKAT(as.matrix(combined_variants_phenotype_bla_dummy[,1:number_variants]), SKAT_result_null, method="SKATO")
            } else if ( nrow(combined_variants_phenotype_bla[combined_variants_phenotype_bla$sex=='female',]) > 2 & nrow(combined_variants_phenotype_bla[combined_variants_phenotype_bla$sex=='male',]) > 2 ) {
                # Need to make dummy variable for sex and tissue type
                combined_variants_phenotype_bla_dummy <- fastDummies::dummy_cols(combined_variants_phenotype_bla, select_columns = c("sex","age"), remove_first_dummy = TRUE) %>%
                    dplyr::select(-tissue,-sex,-age,-GTEx_barcode);
                # Make the null model first, regressing the covariates on the phenotype
                SKAT_result_null <- SKAT_Null_Model(as.matrix(combined_variants_phenotype_bla_dummy[,phenotype]) ~ as.matrix(combined_variants_phenotype_bla_dummy[,(number_variants+2):ncol(combined_variants_phenotype_bla_dummy)]), out_type="C")
                # Perform the SKATO test
                SKAT_result <- SKAT(as.matrix(combined_variants_phenotype_bla_dummy[,1:number_variants]), SKAT_result_null, method="SKATO")
            } else {
                # Need to make dummy variable for sex and tissue type
                combined_variants_phenotype_bla_dummy <- combined_variants_phenotype_bla %>%
                    dplyr::select(-tissue,-sex,-age,-GTEx_barcode)
                # Make the null model first, regressing the covariates on the phenotype
                SKAT_result_null <- SKAT_Null_Model(as.matrix(combined_variants_phenotype_bla_dummy[,phenotype]) ~ as.matrix(combined_variants_phenotype_bla_dummy[,(number_variants+2):ncol(combined_variants_phenotype_bla_dummy)]), out_type="C")
                # Perform the SKATO test
                SKAT_result <- SKAT(as.matrix(combined_variants_phenotype_bla_dummy[,1:number_variants]), SKAT_result_null, method="SKATO")                  
            } # Finished SKAT  
            
            # Result from regression WITHOUT any interaction term
            # Estimate median of logs for plotting purposes later
            ResVariantsPValues[i,1] <- gene_to_test
            ResVariantsPValues$pValue[i] <- SKAT_result$p.value
            ResVariantsPValues$n_carriers[i] <- variants_genes_higher2$numberSamples[i] #count how many samples are carriers
            ResVariantsPValues$n_variants[i] <- number_variants #count how many different variants exist
            if ( number_variants > 1 ) {
                ResVariantsPValues$double_del_carriers[i] <- length(which(apply(combined_variants_phenotype_bla_dummy[,1:number_variants]==2, 1, any))) #count how many doubles exist  
            } else {
                ResVariantsPValues$double_del_carriers[i] <- length(which(combined_variants_phenotype_bla_dummy[,1]==2)) #count how many doubles exist
            }
            # Rho value can be 0 if only one variant
            if ( is.null(SKAT_result$param$rho_est) ) {
                ResVariantsPValues$rho_SKAT_O[i] <- NA
            } else {
                ResVariantsPValues$rho_SKAT_O[i] <- SKAT_result$param$rho_est #rho value, if close to 0 it's SKAT, if close to 1 it's burden test
            }
        } # Finished testing
        ResVariantsPValues$pValue <- as.numeric(ResVariantsPValues$pValue)
        
        ######### output #########
        results_regression <- ResVariantsPValues %>%
            dplyr::filter(!is.na(pValue))
        genes_n <- nrow(results_regression)
        results_regression <- results_regression %>%
            mutate(n_genes= genes_n)
        
        #create new output directory if not already existing
        setwd(output_path)
        if (!file.exists('results')) {
            dir.create('results')
        }
        write.table(results_regression[order(results_regression$pValue),], file= paste(output_path,'/results/',NMD_method_VAF,"_",tissue_type,'_',genes_n,'genes_',samples_n,'samples_',numberGenesThres,"_genes_",selected_genes_info,"_randomization_",randomization,'.txt',sep=''),
                    quote=FALSE, sep='\t',row.names=FALSE,col.names = T)
        
        ######### testing summary #########
        testing_summary <- data.frame(tissue = tissue_type, phenotype = phenotype, n_sample = samples_n, n_genes = genes_n,
                                        cohort = 'TCGA', param = rare_germline_variants_name, model = 'additive', test = 'SKAT-O', stringsAsFactors = FALSE)
        
        #create new output directory if not already existing
        if (!file.exists('test_stats')) {
            dir.create('test_stats')
        }
        write.table(testing_summary,file= paste(output_path,'/test_stats/',
                                                '/',NMD_method_VAF,"_",tissue_type,'_tissue_',genes_n,'genes_',samples_n,'samples_',numberGenesThres,"_genes_",selected_genes_info,"_randomization_",randomization,'.txt',sep=''),
                    quote=FALSE, sep='\t',row.names=FALSE,col.names = T)
    } # Check that there are enough samples left to perform test and enough samples with activity
} # Test if enough genes to test

################################################################################################

########################################## FINISH ##############################################

################################################################################################