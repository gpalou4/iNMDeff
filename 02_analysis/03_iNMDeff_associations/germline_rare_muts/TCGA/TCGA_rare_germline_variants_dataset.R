
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

################################################################################################

########################################## SCRIPT ##############################################

################################################################################################

args <- commandArgs(trailingOnly=TRUE)
tumor <- args[1]
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
# tumor <- "pancancer"
# rare_germline_variants_name <- "PTV_Missense_CADD25_0.1perc"
# NMD_method <- "ASE"
# VAF <- "0.2"
# geneList <- "yes"
# numberGenesThres <- 2
# randomGenes <- "no"
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
# input_phenotypes <- read.csv(file= "/g/strcombio/fsupek_cancer1/gpalou/Mischan/phenotypes_crg75_least10SNVs_noHPV.txt",head=T,sep ="\t",stringsAsFactors = F)
# ##which phenotype to test
# phenotype_index=as.numeric(args[1]);
# phenotype_index <- 145
# phenotype <- colnames(input_phenotypes)[phenotype_index]
# print(phenotype)

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

# TCGA NMDeff
sample_NMD_efficiencies_TCGA_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/cancers/pancancer/NMD_efficiencies_TCGA.txt"
sample_NMD_efficiencies_TCGA <- read.table(file = sample_NMD_efficiencies_TCGA_path, header = TRUE, sep = "\t")
# Convert some columns to factors
factor_cols <- c("cancer_type","cancer_type_strat","cancer_subtype","LF_remove","purity_remove", "MSI_status",
                "batch_portion","batch_plate","batch_center","batch_vial","TCGA_full_barcode")
sample_NMD_efficiencies_TCGA[factor_cols] <- lapply(sample_NMD_efficiencies_TCGA[factor_cols], factor) 
# Remove samples with <3 PTCs
if (NMD_method %in% c("ASE")) {
    cols_rm <- paste0(NMD_method,"_stopgain_",VAF)
    sample_NMD_efficiencies_TCGA[which(sample_NMD_efficiencies_TCGA[,paste0(NMD_method,"_num_PTCs_",VAF)] < 3),cols_rm] <- NA
}
input_phenotypes <- sample_NMD_efficiencies_TCGA
# Change some columns
colnames(input_phenotypes)[1] <- "sample_short"
input_phenotypes$cancer_type <- gsub("TCGA-","",input_phenotypes$cancer_type)

# 1.2) Rare germline variants
variants_gnomad_allinfo <- read.csv(file= paste('/g/strcombio/fsupek_cancer1/gpalou/Mischan/TCGA_germline_input_variants_',rare_germline_variants_name,'.txt',sep=''),
                                    head=T,sep ="\t",stringsAsFactors = F)
                       
# 1.3) LOH
LOH_cnv_FACET_final <- read.csv(file= paste('/g/strcombio/fsupek_cancer1/gpalou/Mischan/TCGA_germline_input_variants_',rare_germline_variants_name,'_LOH.txt',sep=''),
                                head=T,sep ="\t",stringsAsFactors = F)

# 1.4) PCA from common variants and other covariates
PCA_commonVariants <- read.csv(file= "/g/strcombio/fsupek_cancer1/gpalou/Mischan/TCGA_5percent_common_ingnomAD_variants_PCA_samples_5perc_europeans.txt",
                            head=T,sep ="\t",stringsAsFactors = F) %>%
  mutate(cancer_type=sub('TCGA-','',project_id)) %>%
  mutate(sample_short=sub('_','-',sample_short)) %>%
  mutate(sample_short=sub('_','-',sample_short)) %>%
  dplyr::filter(sample_short %in% input_phenotypes$sample_short) %>%
  dplyr::select(sample_short,gender,age,Dim.1,Dim.2,Dim.3,Dim.4,Dim.5,Dim.6,Dim.7,Dim.8,Dim.9,Dim.10) 
PCA_commonVariants$age[is.na(PCA_commonVariants$age)] <- median(PCA_commonVariants$age,na.rm=T) #replace na age with median (102 NA's)

# 1.5) Selected genes to test
# SelectedGenes <- read.csv(file = '/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/NMD_genes.txt',head=TRUE,sep ="\t")
if (isTRUE(genelist) & randomGenes == "yes") {
    SelectedGenes <- read.table(file = "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/random_genes_not_cancer_not_NMD_related_expanded.txt",
                                header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    selected_genes_info <- "random"
} else if (isTRUE(genelist) & randomGenes == "no") {
    SelectedGenes <- read.table(file = "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/NMD_genes_STRING_expanded.txt",
                            header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    selected_genes_info <- "selected"
} else if (!isTRUE(genelist)) {
    selected_genes_info <- "all"
}

# 1.6) Esophagus info
TCGA_eso_info <- read.csv(file = "/g/strcombio/fsupek_cancer1/gpalou/Mischan/TCGA_oesophageal_2017.csv",head=T,sep =";")

ESCA_EAC <-
  TCGA_eso_info %>%
  dplyr::filter(Disease.code == 'ESCA') %>%
  dplyr::filter(Histological.Type == 'AC');

ESCA_ESCC <-
  TCGA_eso_info %>%
  dplyr::filter(Disease.code == 'ESCA') %>%
  dplyr::filter(Histological.Type == 'ESCC');

# 2) Create dataset
#tumor <- "pancancer"
variant_threshold <- 0.1
print(tumor)
#sampleinfo <- paste('TCGA_cadd15_',as.character(variant_threshold),'perc_',sep='')
#sampleinfo <- paste(sampleinfo,phenotype,sep='')

# c('Pancan','OV','GBM','LUAD','LUSC','UCEC','BLCA','PAAD','KIRP','STAD','LIHC','CESC','SARC',
#                'BRCA','COAD','SKCM','KIRC','READ','LGG','HNSC','ESCA','ESCA_EAC',
#                'STAD_ESCA_EAC','KIRC_KIRP','COAD_READ','PRAD')

# Take randomized phenotype?
if ( (randomization == "yes") ) {
    phenotype <- paste0(phenotype,"_randomized")
} 

####create template table with samples which should have a value
if (tumor == 'COAD_READ') {
    table_pheno <-
        input_phenotypes %>%
        dplyr::filter(cancer_type == 'COAD' | cancer_type == 'READ') %>%
        dplyr::select(sample_short,phenotype,cancer_type) %>%
        dplyr::filter(!is.na(!!rlang::sym(phenotype)) ) 
} else if (tumor == 'ESCA_EAC') {
    table_pheno <-
        input_phenotypes %>%
        dplyr::filter(cancer_type == 'ESCA') %>%
        mutate(sample_check=if_else(sample_short %in% ESCA_EAC$barcode, 'yes', 'no'
        )
        ) %>%
        dplyr::filter(sample_check == 'yes') %>%
        dplyr::select(-sample_check) %>%
        dplyr::select(sample_short,phenotype,cancer_type) %>%
        dplyr::filter(!is.na(!!rlang::sym(phenotype)) ) 
} else if(tumor == 'STAD_ESCA_EAC') {
    table_pheno <-
        input_phenotypes %>%
        dplyr::filter(cancer_type == 'STAD' | cancer_type== 'ESCA')%>%
        mutate(sample_check=if_else(sample_short %in% ESCA_ESCC$barcode, 'yes', 'no' ##take the samples which are not ESCC
        )
        ) %>%
        dplyr::filter(sample_check == 'no') %>%
        dplyr::select(-sample_check)  %>%
        dplyr::select(sample_short,phenotype,cancer_type) %>%
        dplyr::filter(!is.na(!!rlang::sym(phenotype)) )
} else if(tumor == 'pancancer') {
    table_pheno <-
        input_phenotypes %>%
        dplyr::select(sample_short,phenotype,cancer_type) %>%
        dplyr::filter(!is.na(!!rlang::sym(phenotype)) )
} else if(tumor == 'KIRC_KIRP') {
    table_pheno <-
        input_phenotypes %>%
        dplyr::filter(cancer_type == 'KIRC'|cancer_type=='KIRP') %>%
        dplyr::select(sample_short,phenotype,cancer_type) %>%
        dplyr::filter(!is.na(!!rlang::sym(phenotype)) )
} else {
    table_pheno <-
        input_phenotypes %>%
        dplyr::filter(cancer_type == tumor) %>%
        dplyr::select(sample_short,phenotype,cancer_type) %>%
        dplyr::filter(!is.na(!!rlang::sym(phenotype)) )
}

if ( !isTRUE(genelist) ) {
    genelist_info <- 'nogenelist'
    #sampleinfo <- paste(sampleinfo,'nogenelist',sep='_')
} else {
    #sampleinfo <- paste(sampleinfo,'withgenelist',sep='_')
    genelist_info <- 'withgenelist'
}
####make sure that at least 5% of the samples have a value for the phenotype!
# pheno_activity <- nrow(table_pheno[table_pheno[[phenotype]] != 0,])/nrow(table_pheno)
# if(pheno_activity > 0.0 & nrow(table_pheno) >= 50){

#saving ratios
matrix_analysis <- table_pheno
print(paste(nrow(matrix_analysis), ' samples',sep=''))

#create new output directory if not already existing
output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/rare_germline_variants/TCGA/",genelist_info,"/",rare_germline_variants_name)
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

# At least one gene to test in the cancer
if ( nrow(variants_genes_higher2) > 1) {
    
    ##make matrix also considering LOH
    germline_variants_matrix <- variants_gnomad_samples %>%
       dplyr::filter(Gene.refGene %in% variants_genes_higher2$Gene.refGene) %>%
       left_join(LOH_cnv_FACET_final,by=c('sample_short','Gene.refGene'))
    #If NA, then "no"
    germline_variants_matrix$LOH[is.na(germline_variants_matrix$LOH)] <- 'no'

    germline_variants_matrix <- germline_variants_matrix %>%
        mutate(geno=if_else(gt == '1/1' | gt == '1|1',2,
                            if_else(LOH == 'yes',2,1 #put in 2, else 1 for hetero and homo
                            )
        )
        ) %>% 
        dplyr::select(sample_short,geno,gene_variant_id) %>%
        distinct() %>% # remove any duplicates, should not exist though (checked)
        spread(gene_variant_id,geno,is.na(0)) %>%
        mutate(samplecheck=if_else(sample_short %in% matrix_analysis$sample_short, 'yes','no' ##only samples with pheno
        )
        ) %>%
        dplyr::filter(samplecheck=='yes') %>%
        dplyr::select(-samplecheck)
    
    #### TEMP ####
    # germline_variants_matrix <- variants_gnomad_samples %>%
    #     dplyr::filter(Gene.refGene %in% variants_genes_higher2$Gene.refGene)
    # germline_variants_matrix <- germline_variants_matrix %>%
    #     mutate(geno=if_else(gt == '1/1' | gt == '1|1',2,1)
    #     ) %>% 
    #     dplyr::select(sample_short,geno,gene_variant_id) %>%
    #     distinct() %>% # remove any duplicates, should not exist though (checked)
    #     spread(gene_variant_id,geno,is.na(0)) %>%
    #     mutate(samplecheck=if_else(sample_short %in% matrix_analysis$sample_short, 'yes','no' ##only samples with pheno
    #     )
    #     ) %>%
    #     dplyr::filter(samplecheck=='yes') %>%
    #     dplyr::select(-samplecheck)

    ##############

    # Use list of samples which have germline variants called and 
    matrix_variants <-  data.frame(sample_short=matrix_analysis$sample_short, stringsAsFactors=F) %>%
        left_join(germline_variants_matrix,by=c('sample_short'))
    # Put all samples without any rare variants to 0 --> NAs to zero
    for ( columns in 2:ncol(matrix_variants) ) {
        matrix_variants[,columns][is.na(matrix_variants[,columns])] <- 0
    }
    # Create the final matrix for testing
    combined_variants_phenotype <- matrix_variants %>%
        left_join(matrix_analysis,by=c('sample_short')) %>%
        left_join(PCA_commonVariants,by=c('sample_short')) 
    
    NumberVariantColumns <- nrow(variants_genes_higher2)
    samples_n <- nrow(combined_variants_phenotype)
    
    # 2.3) Rare variants association testing
    ############## start with testing ################################################################################################################
    # Super important to initialize this matrix every time from scratch
    ResVariantsPValues <- data.frame(variant=rep('sample',NumberVariantColumns),
                                    pValue=rep(NA,NumberVariantColumns),
                                    pheno=rep(phenotype,NumberVariantColumns),
                                    n_carriers=rep(NA,NumberVariantColumns),
                                    n_variants=rep(NA,NumberVariantColumns),
                                    n_samples=rep(samples_n,NumberVariantColumns),
                                    cancer_type=rep(tumor,NumberVariantColumns),
                                    snps=rep(rare_germline_variants_name,NumberVariantColumns),
                                    model=rep('additive',NumberVariantColumns),
                                    test=rep('SKAT-O',NumberVariantColumns),
                                    rho_SKAT_O=rep(NA,NumberVariantColumns),
                                    double_del_carriers=rep(NA,NumberVariantColumns),
                                    stringsAsFactors = FALSE)
    
    if ( samples_n > 2 ) { #only if at least 50 samples go on 
        # Loop through all genes
        for(i in 1:NumberVariantColumns){           
            gene_to_test <- variants_genes_higher2$Gene.refGene[i] #gene to test
            # Grep the column names of the variants inside that gene
            gene_variants_to_test <- colnames(combined_variants_phenotype)[grepl(paste('^',gene_to_test,'__',sep=''),colnames(combined_variants_phenotype))] 
            number_variants <- length(gene_variants_to_test)
            # Select columns needed
            combined_variants_phenotype_bla <- combined_variants_phenotype %>%
                dplyr::select(all_of(gene_variants_to_test),phenotype,gender,age,Dim.1,Dim.2,Dim.3,Dim.4,Dim.5,Dim.6,cancer_type)
            # Test cancer type if more than 1 cancer type in it, same for gender
            if ( length(unique(combined_variants_phenotype_bla$cancer_type)) > 1 & nrow(combined_variants_phenotype_bla[combined_variants_phenotype_bla$gender=='female',]) > 2 & nrow(combined_variants_phenotype_bla[combined_variants_phenotype_bla$gender=='male',]) > 2 ) {
                # Need to make dummy variable for gender and cancer type
                combined_variants_phenotype_bla_dummy <- fastDummies::dummy_cols(combined_variants_phenotype_bla, select_columns = c("cancer_type","gender"), remove_first_dummy = TRUE) %>%
                    dplyr::select(-cancer_type,-gender);
                # Make the null model first, regressing the covariates on the phenotype
                SKAT_result_null <- SKAT_Null_Model(as.matrix(combined_variants_phenotype_bla_dummy[,phenotype]) ~ as.matrix(combined_variants_phenotype_bla_dummy[,(number_variants+2):ncol(combined_variants_phenotype_bla_dummy)]), out_type="C")
                # Perform the SKATO test
                SKAT_result <- SKAT(as.matrix(combined_variants_phenotype_bla_dummy[,1:number_variants]), SKAT_result_null, method="SKATO")
            } else if ( length(unique(combined_variants_phenotype_bla$cancer_type)) > 1 ) {
                # Need to make dummy variable for gender and cancer type
                combined_variants_phenotype_bla_dummy <- fastDummies::dummy_cols(combined_variants_phenotype_bla, select_columns = c("cancer_type"), remove_first_dummy = TRUE) %>%
                    dplyr::select(-cancer_type,-gender);
                # Make the null model first, regressing the covariates on the phenotype
                SKAT_result_null <- SKAT_Null_Model(as.matrix(combined_variants_phenotype_bla_dummy[,phenotype]) ~ as.matrix(combined_variants_phenotype_bla_dummy[,(number_variants+2):ncol(combined_variants_phenotype_bla_dummy)]), out_type="C")
                # Perform the SKATO test
                SKAT_result <- SKAT(as.matrix(combined_variants_phenotype_bla_dummy[,1:number_variants]), SKAT_result_null, method="SKATO")
            } else if ( nrow(combined_variants_phenotype_bla[combined_variants_phenotype_bla$gender=='female',]) > 2 & nrow(combined_variants_phenotype_bla[combined_variants_phenotype_bla$gender=='male',]) > 2 ) {
                # Need to make dummy variable for gender and cancer type
                combined_variants_phenotype_bla_dummy <- fastDummies::dummy_cols(combined_variants_phenotype_bla, select_columns = c("gender"), remove_first_dummy = TRUE) %>%
                    dplyr::select(-cancer_type,-gender);
                # Make the null model first, regressing the covariates on the phenotype
                SKAT_result_null <- SKAT_Null_Model(as.matrix(combined_variants_phenotype_bla_dummy[,phenotype]) ~ as.matrix(combined_variants_phenotype_bla_dummy[,(number_variants+2):ncol(combined_variants_phenotype_bla_dummy)]), out_type="C")
                # Perform the SKATO test
                SKAT_result <- SKAT(as.matrix(combined_variants_phenotype_bla_dummy[,1:number_variants]), SKAT_result_null, method="SKATO")
            } else {
                # Need to make dummy variable for gender and cancer type
                combined_variants_phenotype_bla_dummy <- combined_variants_phenotype_bla %>%
                    dplyr::select(-cancer_type,-gender)
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
        write.table(results_regression[order(results_regression$pValue),], file= paste(output_path,'/results/',NMD_method_VAF,"_",tumor,'_',genes_n,'genes_',samples_n,'samples_',numberGenesThres,"_genes_",selected_genes_info,"_randomization_",randomization,'.txt',sep=''),
                    quote=FALSE, sep='\t',row.names=FALSE,col.names = T)
        
        ######### testing summary #########
        testing_summary <- data.frame(cancer_type= tumor, phenotype=phenotype,n_sample=samples_n,n_genes=genes_n,
                                        cohort='TCGA',param=rare_germline_variants_name, model='additive',test='SKAT-O',stringsAsFactors = FALSE)
        
        #create new output directory if not already existing
        if (!file.exists('test_stats')) {
            dir.create('test_stats')
        }
        write.table(testing_summary,file= paste(output_path,'/test_stats/',
                                                '/',NMD_method_VAF,"_",tumor,'_cancer_type_',genes_n,'genes_',samples_n,'samples_',numberGenesThres,"_genes_",selected_genes_info,"_randomization_",randomization,'.txt',sep=''),
                    quote=FALSE, sep='\t',row.names=FALSE,col.names = T)
    } # Check that there are enough samples left to perform test and enough samples with activity
} # Test if enough genes to test

################################################################################################

########################################## FINISH ##############################################

################################################################################################