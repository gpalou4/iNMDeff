# No conda
library("ggh4x")
library("GWASTools")
library("VennDiagram")
library("dplyr")
library("cowplot")
library("ggplot2")
library("RColorBrewer")

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

matched_tissue_RGWAS_replication_strat_tissue <- function(RGWAS_df_all, randomization_char, FDR, database_discovery, database_validation) {
    
    FDR <- FDR / 100
    # randomization_char <- randomization
    RGWAS_replicated_final_df <- c()
    # 1) Sig hits in discovery
    RGWAS_discovery_sig <- RGWAS_df_all %>%
                    filter(randomization == randomization_char) %>%
                    group_by(database, tissue) %>%
                    mutate(pvalue_FDR_adjusted = p.adjust(pValue, method = "fdr")) %>%
                    mutate(sig_genes_ID = ifelse( ((pvalue_FDR_adjusted < FDR) ), paste0(variant,"_",tissue), NA)) %>%
                    mutate(discovery_hits_ID = ifelse( ((pvalue_FDR_adjusted < FDR) ), paste0(database,"_",variant,"_",NMD_method,"_",dataset,"_",tissue), NA))
    sig_genes_ID_char <- as.character(na.omit(RGWAS_discovery_sig$sig_genes_ID))

    df <- RGWAS_discovery_sig[RGWAS_discovery_sig$variant == "KDM6B" & RGWAS_discovery_sig$database == "TCGA",]
    data.frame(df)
    df <- df[order(df$pValue),]
    data.frame(df)

    # 2) Validation of hits in matched discovery tissue
    for (i in 1:nrow(TCGA_GTEx_tissues)) {
        tissue <- TCGA_GTEx_tissues[i,"tissues"]
        TCGA_cancer <- TCGA_GTEx_tissues[i,"TCGA_cancers"]
        TCGA_cancer <- gsub("_MSS|_muscle","",TCGA_cancer)
        GTEx_tissue <- TCGA_GTEx_tissues[i,"GTEx_tissues"]
        if (GTEx_tissue == "pantissue") {
            GTEx_tissue_full <- "pantissue"
        } else {
            GTEx_tissue_full <- paste0(as.character(tissues_acronyms[grep(GTEx_tissue,tissues_acronyms$acronyms),"GTEx_tissues_full"]),collapse="|")
        }
        print(TCGA_GTEx_tissues[i,])

        if (database_discovery == "TCGA") {
            discovery_tissue <- TCGA_cancer
            validation_tissue <- GTEx_tissue_full
        } else {
            discovery_tissue <- GTEx_tissue_full
            validation_tissue <- TCGA_cancer
        }
        # Sig hits IDs
        sig_genes_ID_char_to_validate <- sig_genes_ID_char[grep(discovery_tissue,sig_genes_ID_char)]
        discovery_tissues <- unlist(strsplit(discovery_tissue,"\\|"))
        validation_tissues <- unlist(strsplit(validation_tissue,"\\|"))
        sig_genes_ID_char_tissue <- c()
        for (j in 1:length(discovery_tissues)) {
            for (i in 1:length(validation_tissues)) {
                tmp <- gsub(discovery_tissues[j],validation_tissues[i],sig_genes_ID_char_to_validate)
                sig_genes_ID_char_tissue <- c(sig_genes_ID_char_tissue,tmp)
            }
        }
        index <- grep(discovery_tissue,sig_genes_ID_char_tissue)
        if (length(index) != 0) {
            sig_genes_ID_char_tissue <- sig_genes_ID_char_tissue[-grep(discovery_tissue,sig_genes_ID_char_tissue)]
        }
        sig_genes_ID_char_tissue <- sig_genes_ID_char_tissue[!duplicated(sig_genes_ID_char_tissue)]
        # FDR in validation (only for the subset of genes that are significant in discovery)
        RGWAS_FDR_validation <- RGWAS_discovery_sig %>%
            filter(database == database_validation & randomization == randomization_char) %>%
            filter(paste0(variant,"_",tissue) %in% sig_genes_ID_char_tissue) %>%
            group_by(tissue) %>%
            mutate(pvalue_FDR_adjusted = p.adjust(pValue, method = "fdr")) %>%
            filter(pvalue_FDR_adjusted < FDR) %>%
            mutate(replicated_hits_ID = paste0(database,"_",variant,"_",NMD_method,"_",dataset,"_",tissue))
        RGWAS_validation_sig_replicated <- RGWAS_FDR_validation
        if (nrow(RGWAS_validation_sig_replicated) == 0 ) {next}
        print(RGWAS_validation_sig_replicated)
        # 3) Obtain the info of the replicated hits in discovery and validation
        # Dataset & NMD_method
        # Discovery hits ID
        RGWAS_discovery_sig_hits_ID <- RGWAS_discovery_sig %>%
                        filter(database == database_discovery & randomization == randomization_char) %>%
                        filter(!is.na(discovery_hits_ID)) %>%
                        mutate(tmp_ID = paste0(variant,"_",tissue))
        # Check which validated genes-tissue were in discovery
        sig_genes_ID_char_to_validate <- sig_genes_ID_char[grep(discovery_tissue,sig_genes_ID_char)]
        sig_genes_ID_char_tissue <- c()
        for (j in 1:length(discovery_tissues)) {
            tmp <- paste0(RGWAS_validation_sig_replicated$variant,"_",discovery_tissues[j])
            sig_genes_ID_char_tissue <- c(sig_genes_ID_char_tissue,tmp)
        }
        # Discovery Info
        RGWAS_discovery_sig_replicated <- RGWAS_discovery_sig_hits_ID[RGWAS_discovery_sig_hits_ID$tmp_ID %in% sig_genes_ID_char_tissue,]
        # 4) Merge Replicated hits: from Discovery & Validated
        RGWAS_final_replicated_hits <- rbind(RGWAS_discovery_sig_replicated,RGWAS_validation_sig_replicated)
        RGWAS_final_replicated_hits$matched_tissue <- tissue
        RGWAS_final_replicated_hits$database_discovery <- database_discovery
        RGWAS_final_replicated_hits$database_validation <- database_validation
        RGWAS_final_replicated_hits$disc_val <- NA
        RGWAS_final_replicated_hits <- RGWAS_final_replicated_hits %>%
                                        mutate(disc_val = ifelse( ((database %in% database_discovery) ), "discovery", "validation"))
        # Number of genes tested at discovery
        num_genes_tested <- RGWAS_df_all %>%
                            filter(randomization == randomization_char & database == database_discovery) %>%
                        summarise(genes_tested = n()) %>% data.frame()
        RGWAS_final_replicated_hits$genes_tested <- num_genes_tested$genes_tested
        RGWAS_replicated_final_df <- rbind(RGWAS_replicated_final_df,RGWAS_final_replicated_hits)
    }
    return(RGWAS_replicated_final_df)
}

cross_tissue_RGWAS_replication_strat_tissue <- function(RGWAS_df_all, randomization_char, FDR, database_discovery, database_validation) {
    FDR <- FDR / 100
    # randomization_char <- randomization
    RGWAS_replicated_final_df <- c()
    # 1) Sig hits in discovery
    RGWAS_discovery_sig <- RGWAS_df_all %>%
                    filter(randomization == randomization_char) %>%
                    group_by(database, tissue) %>%
                    mutate(pvalue_FDR_adjusted = p.adjust(pValue, method = "fdr")) %>%
                    mutate(sig_genes_ID = ifelse( ((pvalue_FDR_adjusted < FDR) ), paste0(variant,"_",tissue), NA))
    sig_genes_ID_char <- as.character(na.omit(RGWAS_discovery_sig$sig_genes_ID))

    # 2) Validation of hits in any discovery tissue
    for (i in 1:nrow(TCGA_GTEx_tissues)) {
        tissue <- TCGA_GTEx_tissues[i,"tissues"]
        TCGA_cancer <- TCGA_GTEx_tissues[i,"TCGA_cancers"]
        TCGA_cancer <- gsub("_MSS|_muscle","",TCGA_cancer)
        GTEx_tissue <- TCGA_GTEx_tissues[i,"GTEx_tissues"]
        if (GTEx_tissue == "pantissue") {
            GTEx_tissue_full <- "pantissue"
        } else {
            GTEx_tissue_full <- paste0(as.character(tissues_acronyms[grep(GTEx_tissue,tissues_acronyms$acronyms),"GTEx_tissues_full"]),collapse="|")
        }
        print(TCGA_GTEx_tissues[i,])

        if (database_discovery == "TCGA") {
            discovery_tissue <- TCGA_cancer
            validation_tissue <- GTEx_tissue_full
        } else {
            discovery_tissue <- GTEx_tissue_full
            validation_tissue <- TCGA_cancer
        }
        # Sig hits IDs
        sig_genes_ID_char_to_validate <- sig_genes_ID_char[grep(discovery_tissue,sig_genes_ID_char)]
        discovery_tissues <- unlist(strsplit(discovery_tissue,"\\|"))
        validation_tissues <- unlist(strsplit(validation_tissue,"\\|"))
        sig_genes_ID_char_tissue <- c()
        for (j in 1:length(discovery_tissues)) {
            for (i in 1:length(validation_tissues)) {
                tmp <- gsub(discovery_tissues[j],validation_tissues[i],sig_genes_ID_char_to_validate)
                sig_genes_ID_char_tissue <- c(sig_genes_ID_char_tissue,tmp)
            }
        }
        index <- grep(discovery_tissue,sig_genes_ID_char_tissue)
        if (length(index) != 0) {
            sig_genes_ID_char_tissue <- sig_genes_ID_char_tissue[-grep(discovery_tissue,sig_genes_ID_char_tissue)]
        }
        sig_genes_ID_char_tissue <- sig_genes_ID_char_tissue[!duplicated(sig_genes_ID_char_tissue)]
        # FDR in validation (only for the subset of genes that are significant in discovery)
        RGWAS_FDR_validation <- RGWAS_discovery_sig %>%
            filter(database == database_validation & randomization == randomization_char) %>%
            group_by(tissue) %>%
            filter(paste0(variant,"_",tissue) %in% sig_genes_ID_char_tissue) %>%
            mutate(pvalue_FDR_adjusted = p.adjust(pValue, method = "fdr")) %>%
            filter(pvalue_FDR_adjusted < FDR)
        #RGWAS_validation_sig <- data.frame(RGWAS_empirical_FDR_validation[RGWAS_empirical_FDR_validation$empirical_FDR < empirical_FDR,])
        RGWAS_validation_sig <- RGWAS_FDR_validation
        if (nrow(RGWAS_validation_sig) == 0 ) {next}
        print(RGWAS_validation_sig)
        # Save
        #colnames(RGWAS_validation_sig)[colnames(RGWAS_validation_sig) == "tissue"] <- "TCGA_cancer"
        RGWAS_validation_sig$matched_tissue <- tissue
        RGWAS_validation_sig$database_discovery <- database_discovery
        RGWAS_validation_sig$database_validation <- database_validation
        # Number of genes tested at discovery
        num_genes_tested <- RGWAS_df_all %>%
                            filter(randomization == randomization_char & database == database_discovery) %>%
                        summarise(genes_tested = n()) %>% data.frame()
        RGWAS_validation_sig$genes_tested <- num_genes_tested$genes_tested
        RGWAS_replicated_final_df <- rbind(RGWAS_replicated_final_df,RGWAS_validation_sig)
    }
    return(RGWAS_replicated_final_df)
}

venn_diagram_within_databases <- function(list, title) {
    vd <- venn.diagram(
            x = list,
            category.names = c("ASE" , "Endogenous"),
            filename = NULL,
            output=F,
            main = title,
            # Output features
            imagetype="png" ,
            height = 320 , 
            width = 320 , 
            resolution = 300,
            compression = "lzw",
            # Circles
            lwd = 1,
            lty = 'blank',
            fill = myCol,
            # Numbers
            cex = 1,
            fontface = "bold",
            #fontfamily = "sans",
            # Set names
            cat.cex = 0.7,
            cat.fontface = "bold",
            cat.default.pos = "outer",
            cat.pos = c(-27, 27),
            cat.dist = c(0.055, 0.055),
            cat.fontfamily = "sans"
            #rotation = 1
    )
    return(vd)
}

calculate_OR <- function(vd_list_genes) {
    # OR
    num_ASE_hits <- length(vd_list_genes$ASE_hits)
    num_END_hits <- length(vd_list_genes$END_hits)
    num_ASE_END_hits <- length(vd_list_genes$ASE_END_overlapping_hits)
    num_total <- length(vd_list_genes$testable_genes)
    # Create a 2x2 contingency table
    observed <- matrix(c(num_ASE_END_hits, c(num_ASE_hits-num_ASE_END_hits), num_END_hits-num_ASE_END_hits, 
                    c(num_total-num_ASE_hits-num_END_hits+num_ASE_END_hits)), nrow = 2, byrow = TRUE)
    # Perform Fisher's exact test and store the result
    result <- fisher.test(observed)
    # Print the odds ratio
    odds_ratio <- result$estimate
    cat("Odds Ratio:", odds_ratio, "\n")
    # Print the p-value
    p_value <- result$p.value
    cat("p-value:", p_value, "\n")
    return(round(odds_ratio,2))
}

sig_hits_within_databases <- function(RGWAS_df_all, ddbb, FDR, randomization_char) {
    FDR_num <- as.numeric(FDR/100)
    # Sig hits in ASE
    RGWAS_ASE_sig <- RGWAS_df_all %>%
                    filter(randomization == randomization_char & database == ddbb & NMD_method == "ASE") %>%
                    group_by(tissue,dataset,NMD_method) %>%
                    mutate(pvalue_FDR_adjusted = p.adjust(pValue, method = "bonferroni")) %>%
                    filter( pvalue_FDR_adjusted < FDR_num )
    # Sig hits in Endogenous
    RGWAS_END_sig <- RGWAS_df_all %>%
                    filter(randomization == randomization_char & database == ddbb & NMD_method == "Endogenous") %>%
                    group_by(tissue,dataset,NMD_method) %>%
                    mutate(pvalue_FDR_adjusted = p.adjust(pValue, method = "bonferroni")) %>%
                    filter( pvalue_FDR_adjusted < FDR_num )
    # Overlapping hits
    ASE_hits <- as.character(unique(RGWAS_ASE_sig$variant))
    END_hits <- as.character(unique(RGWAS_END_sig$variant))
    ASE_END_overlapping_hits <- intersect(ASE_hits,END_hits)

    RGWAS_END_sig %>% filter(variant %in% ASE_END_overlapping_hits) %>%
                    group_by(tissue,dataset,NMD_method) %>%
                    summarise(n())
    RGWAS_ASE_sig %>% filter(variant %in% ASE_END_overlapping_hits) %>%
                    group_by(tissue,dataset,NMD_method) %>%
                    summarise(n())
    # Total testable genes
    ASE_genes <- RGWAS_df_all %>%
                    filter(randomization == randomization_char & database == ddbb & NMD_method == "ASE") %>%
                    pull(variant)
    END_genes <- RGWAS_df_all %>%
                    filter(randomization == randomization_char & database == ddbb & NMD_method == "Endogenous") %>%
                    pull(variant)
    testable_genes <- ASE_genes
    vd_genes <- list(ASE_hits = ASE_hits, END_hits = END_hits, ASE_END_overlapping_hits = ASE_END_overlapping_hits, testable_genes = testable_genes )
    return(vd_genes)
}

GTEx_tissue_names_path <- "/g/strcombio/fsupek_cancer1/gpalou/GTEx/RNAseq/GTEx_tissues_clean.txt"
GTEx_tissues <- read.table(file = GTEx_tissue_names_path, stringsAsFactors = FALSE)$V1
TCGA_cancer_names_path <- "/g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/TCGA_projects_names.txt"
TCGA_cancers <- read.table(file = TCGA_cancer_names_path, stringsAsFactors = FALSE)$V1
NMD_genes <- read.table( file = "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/NMD_genes.txt",
                        sep = "\t", header = TRUE)
output_path <- paste0("/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/TCGA_GTEx_match.txt")
TCGA_GTEx_tissues <- read.table(file = output_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
acronyms <- c("ADPSBQ","ADPVSC","ADRNLG","ARTAORT","ARTCRN","ARTTBL","BLDDER","BRNAMY","BRNACC","BRNCDT","BRNCHB","BRNCHA","BRNCTXA","BRNCTXB",
              "BRNHPP","BRNHPT","BRNNCC","BRNPTM","BRNSPC","BRNSNG","BREAST","FIBRBLS","LCL","CML","CVXECT","CVSEND","CLNSGM","CLNTRN","ESPGEJ","ESPMCS","ESPMSL",
              "FLLPNT","HRTAA","HRTLV","KDNCTX","KDNMDL","LIVER","LUNG","SLVRYG","MSCLSK","NERVET","OVARY","PNCREAS","PTTARY","PRSTTE","SKINNS","SKINS",
              "SNTTRM","SPLEEN","STMACH","TESTIS","THYROID","UTERUS","VAGINA","WHLBLD")
tissues_acronyms <- data.frame(GTEx_tissues_full = GTEx_tissues, acronyms = acronyms )
TCGA_GTEx_tissues[19,] <- c("pantissue","pantissue","pancancer")
GTEx_tissues <- c(GTEx_tissues,"pantissue")
TCGA_cancers <- c(TCGA_cancers,"pancancer")

# 1) Obtain datasets
geneList <- "no"
if (geneList == "yes") {
    genelist <- "withgenelist"
} else if (geneList == "no") {
    genelist <- "nogenelist"
}
databases <- c("GTEx","TCGA")

# GTEx
GTEx_RGWAS_df_all <- c()
for (randomization in c("yes","no")) {
    for (dataset in c("PTV_Missense_CADD15_0.1perc","PTV_Missense_CADD25_0.1perc","PTV_0.1perc")) {
        print(dataset)
        RWGAS_ASE <- RGWAS_df(dataset_name = dataset, NMD_method = "ASE", geneList = geneList, randomGenes = "all",
                    randomization = randomization, geneVariantsThres = 2, VAF = 0.2, tissues = GTEx_tissues, database = "GTEx")
        RWGAS_ASE$NMD_method <- "ASE"
        RWGAS_end <- RGWAS_df(dataset_name = dataset, NMD_method = "endogenous", geneList = geneList, randomGenes = "all",
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

# TCGA
TCGA_RGWAS_df_all <- c()
for (randomization in c("yes","no")) {
    for (dataset in c("PTV_Missense_CADD15_0.1perc","PTV_Missense_CADD25_0.1perc","PTV_0.1perc","MRT_25thCentile_0.1perc","CCR_90orhigher_0.1perc")) {
        print(dataset)
        RWGAS_ASE <- RGWAS_df(dataset_name = dataset, NMD_method = "ASE", geneList = geneList, randomGenes = "all",
                    randomization = randomization, geneVariantsThres = 2, VAF = 0.2, tissues = TCGA_cancers, database = "TCGA")
        RWGAS_ASE$NMD_method <- "ASE"
        RWGAS_end <- RGWAS_df(dataset_name = dataset, NMD_method = "endogenous", geneList = geneList, randomGenes = "all",
                    randomization = randomization, geneVariantsThres = 2, VAF = 0.2, tissues = TCGA_cancers, database = "TCGA")
        RWGAS_end$NMD_method <- "Endogenous"
        if (is.null(nrow(RWGAS_end))) {
            RWGAS_dataset <- RWGAS_ASE
        } else {
            RWGAS_dataset <- rbind(RWGAS_ASE,RWGAS_end)
        }  
        RWGAS_dataset$dataset <- dataset
        RWGAS_dataset$randomization <- randomization
        RWGAS_dataset$random_geneset <- "all"
        if (length(TCGA_RGWAS_df_all) == 0) {
            TCGA_RGWAS_df_all <- RWGAS_dataset
        } else {
            TCGA_RGWAS_df_all <- rbind(TCGA_RGWAS_df_all,RWGAS_dataset)
        }
    }
    
}
print(dim(TCGA_RGWAS_df_all))
colnames(TCGA_RGWAS_df_all)[colnames(TCGA_RGWAS_df_all) %in% "cancer_type"] <- "tissue"

# Merge
TCGA_RGWAS_df_all$database <- "TCGA"
GTEx_RGWAS_df_all$database <- "GTEx"
RGWAS_df_all <- rbind(GTEx_RGWAS_df_all,TCGA_RGWAS_df_all)

RGWAS_df_all <- RGWAS_df_all %>%
                        filter(!dataset %in% c("CCR_90orhigher_0.1perc","MRT_25thCentile_0.1perc"))

# 1.2) Filter models with inflation
# TCGA
input_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/rare_germline_variants/TCGA/",genelist,"/TCGA_",genelist,"_all_lambdas.txt")
TCGA_lambdas <- read.table( file = input_path, sep = "\t", header = TRUE)
# Plot of Lambdas across tissues
p <- ggplot(TCGA_lambdas, aes(x = tissue, y = lambda_non_randomization, color = factor(dataset))) + 
        geom_point(size = 5) +
        ylim(c(0,4)) +
        facet_grid(. ~ NMD_method) +
        theme_bw(base_size = 25) +
        theme(axis.text.x = element_text(size = 14, angle = 45, hjust=1),
            legend.position="top",
            legend.text = element_text(size=12, face="bold"),
            legend.title = element_blank()) +
        guides(colour = guide_legend(override.aes = list(size=10))) 
output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/rare_germline_variants/replication/",genelist,"/TCGA_lambdas.png")
png(output_path, width = 5000, height = 3500, res = 300)
print(p)
dev.off()

# GTEx
input_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/rare_germline_variants/GTEx/",genelist,"/GTEx_",genelist,"_all_lambdas.txt")
GTEx_lambdas <- read.table( file = input_path, sep = "\t", header = TRUE)
# Plot of Lambdas across tissues
p <- ggplot(GTEx_lambdas, aes(x = tissue, y = lambda_non_randomization, color = factor(dataset))) + 
        #geom_bar(stat='identity',position=position_dodge()) +
        geom_point(size = 5) +
        ylim(c(0,4)) +
        facet_grid(. ~ NMD_method) +
        theme_bw(base_size = 25) +
        theme(axis.text.x = element_text(size = 10, angle = 45, hjust=1),
            legend.position="top",
            legend.text = element_text(size=14, face="bold"),
            legend.title = element_blank()) +
        guides(colour = guide_legend(override.aes = list(size=10))) 
output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/rare_germline_variants/replication/",genelist,"/GTEx_lambdas.png")
png(output_path, width = 5500, height = 3500, res = 300)
print(p)
dev.off()

# Filter
TCGA_models_to_remove_non_rand <- TCGA_lambdas[which(TCGA_lambdas$lambda_non_randomization >= 1.5),]
TCGA_models_to_remove_non_rand <- TCGA_models_to_remove_non_rand %>%
        mutate(ID = paste(database,NMD_method,dataset,tissue,"no",sep = "_"))
TCGA_models_to_remove_rand <- TCGA_lambdas[which(TCGA_lambdas$lambda_randomization >= 1.5),]
TCGA_models_to_remove_rand <- TCGA_models_to_remove_rand %>%
        mutate(ID = paste(database,NMD_method,dataset,tissue,"yes",sep = "_"))
GTEx_models_to_remove_non_rand <- GTEx_lambdas[which(GTEx_lambdas$lambda_non_randomization >= 1.5),]
GTEx_models_to_remove_non_rand <- GTEx_models_to_remove_non_rand %>%
        mutate(ID = paste(database,NMD_method,dataset,tissue,"no",sep = "_"))
GTEx_models_to_remove_rand <- GTEx_lambdas[which(GTEx_lambdas$lambda_randomization >= 1.5),]
GTEx_models_to_remove_rand <- GTEx_models_to_remove_rand %>%
        mutate(ID = paste(database,NMD_method,dataset,tissue,"yes",sep = "_"))
RGWAS_res_filt <- RGWAS_df_all %>%
    mutate(ID = paste(database,NMD_method,dataset,tissue,randomization,sep = "_")) %>%
    filter(! (ID %in% TCGA_models_to_remove_non_rand$ID ) ) %>%
    filter(! (ID %in% TCGA_models_to_remove_rand$ID ) ) %>%
    filter(! (ID %in% GTEx_models_to_remove_non_rand$ID ) ) %>%
    filter(! (ID %in% GTEx_models_to_remove_rand$ID ) )
RGWAS_res_filt_copy <- RGWAS_res_filt

# 1.3) Save
output_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/rare_germline_variants/replication/RGWAS_res_filtered.txt"
RGWAS_res_filt <- read.table(file = output_path, header = TRUE, sep = "\t")
# write.table(data.frame(RGWAS_res_filt), file = output_path, sep = "\t", quote = FALSE,
#                       col.names = TRUE, row.names = FALSE)

# 2) Results

# Descriptive analysis
# Number of tests for each database
RGWAS_res_tmp <- RGWAS_res_filt %>%
                    filter(!tissue %in% c("pancancer","pantissue")) %>%
                    group_by(database,randomization) %>%
                    summarise(tests = n()) %>%
                    arrange(desc(tests)) %>% as.data.frame()
RGWAS_res_tmp

# Are significant hits (48) significantly higher than randomized (24) ?

# Given data
N_obs <- 850721 +221613 # Total number of genes tested in the observed dataset
N_rand <- 813312 + 221842  # Total number of genes tested in the randomized dataset

observed_significant <- 48  # Observed number of significant candidates
random_significant <- 24  # Randomized number of significant candidates

# Proportions
p_obs <- observed_significant / N_obs  # Proportion of observed significant genes
p_rand <- random_significant / N_rand  # Proportion of randomized significant genes

# Standard error for the difference in proportions
SE <- sqrt((p_rand * (1 - p_rand)) / N_rand + (p_obs * (1 - p_obs)) / N_obs)

# Z-score calculation
z_score <- (p_obs - p_rand) / SE

# Two-tailed p-value
p_value <- 2 * (1 - pnorm(abs(z_score)))

# Output results
z_score
p_value

# Alternatively use a binomial test OR PERMUTATION test?

# 2.1) Sig hits --> WITHIN databases and BETWEEN methods (any direction)

# 2.2) Sig hits --> BETWEEN databases
# Split by tissue only and matched turmor-normal

databases <- c("GTEx","TCGA")
 
# 2.2.1) Proportion of replicated hits with stratification: only by tissue
RGWAS_replication_hits_final <- c()
for (randomization in c("yes","no")) {
    for (FDR in c(1,2,3,4,5,10)) {
        for (database_discovery in databases) {
            for (database_validation in databases) {
                if (database_discovery == database_validation) {next}
                    RGWAS_replication_hits <- matched_tissue_RGWAS_replication_strat_tissue(RGWAS_df_all = RGWAS_res_filt, randomization_char = randomization,
                                                    FDR = FDR, database_discovery = database_discovery, database_validation = database_validation)
                    if (!is.null(RGWAS_replication_hits)) {
                        RGWAS_replication_hits$FDR_threshold_used <- FDR 
                        RGWAS_replication_hits_final <- rbind(RGWAS_replication_hits_final,RGWAS_replication_hits)
                    }
            }
        }
    }
}
# Save
genelist <- "nogenelist"
output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/rare_germline_variants/replication/",genelist,"/between/RGWAS_replication_hits.txt")
RGWAS_replication_hits_final <- read.table(file = output_path, header = TRUE, sep = "\t")
# write.table(RGWAS_replication_hits_final, file = output_path , sep = "\t", quote = FALSE,
#             col.names = TRUE, row.names = FALSE)
write.table(RGWAS_replication_hits_final, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig4/fig4B_C.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
saveRDS(RGWAS_replication_hits_final, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/Fig4/fig4B_C.RData")

df <- RGWAS_replication_hits_final %>% 
            filter(randomization == "no") %>%
            filter(FDR_threshold_used == 2)
dim(df)
df <- RGWAS_replication_hits_final %>% 
            filter(randomization == "yes") %>%
            filter(FDR_threshold_used == 2)
dim(df)
# The number of hits should not count the number of datasets!
tmp_df <- RGWAS_replication_hits_final[,c("variant","FDR_threshold_used","database_discovery","database_validation","randomization","matched_tissue")]
tmp_df <- tmp_df[!duplicated(tmp_df),]

p <- tmp_df %>%
        filter(!matched_tissue %in% c("pantissue","pancancer")) %>%
        group_by(FDR_threshold_used, database_discovery, database_validation, randomization) %>%
        # group_by(FDR_threshold_used,randomization,database, disc_val) %>%
        summarise(num_replicated_hits = n()) %>%
        mutate(database_validation = paste0(database_validation," - validation")) %>%
        mutate(database_discovery = paste0(database_discovery," - discovery")) %>% data.frame() %>%
            ggplot(aes(x = factor(FDR_threshold_used), y = num_replicated_hits, fill = randomization)) + 
                geom_bar(stat='identity',position=position_dodge()) +
                facet_grid(. ~ database_discovery + database_validation) +
                scale_fill_brewer(palette = "Set2", direction = -1, labels = c("Observed", "Randomization")) +
                labs(title = "", x = "% FDR threshold", y = "Replicated hits", fill = "") + 
                theme_bw(base_size = 35) + scale_x_discrete(labels=c("1","2","3","4","5","10")) +
                theme(axis.text.x = element_text(size = 32),
                    strip.text = element_text(size = 35),
                    axis.text.y = element_text(size = 35),
                    legend.position = "top",
                    legend.text = element_text(size = 35)) +
                guides(fill = guide_legend(override.aes = list(size = 10)))
output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/rare_germline_variants/replication/",genelist,"/between/num_replicated_hits_strat_tissue.png")
png(output_path, width = 5000, height = 3500, res = 300)
print(p)
dev.off()

write.table(tmp_df, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig16/SuppFig16A.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
saveRDS(tmp_df, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig16/SuppFig16A.RData")

# Actual calculated FDR at a 2% empirical FDR 
# Ratio of random/observed hits in one direction (GTEx discovery and TCGA as validation) 
# --> 24/(48) = 0.5. 
FDR_observed_calculated <- tmp_df %>%
            group_by(FDR_threshold_used) %>%
            filter(database_discovery == "GTEx" & database_validation == "TCGA") %>%
            summarise(num_hits_rand = sum(randomization == "yes"), 
                        num_hits_obs = sum(randomization == "no"),
                    FDR_observed = num_hits_rand / num_hits_obs)
FDR_observed_calculated

# 2.2.2) Proportion of replicated hits with stratification: tissue, and validated in any other tissue

# 2.3) Barplots of numbers of replicated hits
# Descriptive
# Number of hits by FDR
# num_hits_FDR <- RGWAS_replication_hits_final %>%
#                 filter(randomization == "no") %>%
#                 group_by(FDR_threshold_used) %>%
#                 summarise(replicated_hits = n()) %>%
#                 arrange(replicated_hits)
# num_hits_FDR
# Number of genes == 19
num_genes <- RGWAS_replication_hits_final %>%
                filter(randomization == "no" & FDR_threshold_used == 2) %>%
                pull(variant)
unique(num_genes)
# Number of gene-tissue pairs == 19
num_genes_tissue <- RGWAS_replication_hits_final %>%
                filter(randomization == "no" & FDR_threshold_used == 2) %>%
                select(variant,matched_tissue)
num_genes_tissue[,c("database","tissue")] <- NULL
num_genes_tissue <- num_genes_tissue[!duplicated(num_genes_tissue),]
unique(num_genes_tissue$variant)

# By all
RGWAS_replication_hits_final_filt <- RGWAS_replication_hits_final %>%
    mutate(database_validation = paste0(database_validation,"_val")) %>%
    mutate(database_discovery = paste0(database_discovery,"_disc"))

p <- RGWAS_replication_hits_final_filt %>%
        group_by(matched_tissue, dataset,NMD_method,disc_val,database) %>%
        filter(FDR_threshold_used == 2 & randomization == "no") %>%
        summarise(replicated_hits = n()) %>%
        arrange(replicated_hits) %>%
        #mutate(matched_tissue = factor(matched_tissue, levels = rev(matched_tissue))) %>%
        ggplot(aes(x = dataset, y = replicated_hits, fill = matched_tissue)) + 
            geom_bar(stat='identity',position=position_dodge()) +
            #ylim(c(0,0.25)) +
            #facet_grid(. ~ NMD_method_discovery + NMD_method_validation) +
            #facet_grid(. ~ database_discovery + database_validation) +
            facet_grid(disc_val ~ NMD_method + database) +
            theme_bw(base_size = 25) +
            theme(axis.text.x = element_text(size = 18, angle = 45, hjust=1))
output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/rare_germline_variants/replication/",genelist,"/between/Number_replicated_hits_by_all.png")
png(output_path, width = 4000, height = 4500, res = 300)
print(p)
dev.off()

# By Dataset
p <- RGWAS_replication_hits_final %>%
        filter(FDR_threshold_used == 2 & randomization == "no") %>%
        group_by(dataset) %>%
        summarise(replicated_hits = n()) %>%
        ggplot(aes(x = dataset, y = replicated_hits)) + 
            geom_bar(stat='identity',position=position_dodge()) +
            theme_bw(base_size = 25) +
            theme(axis.text.x = element_text(size = 18, angle = 45, hjust=1))
output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/rare_germline_variants/replication/",genelist,"/between/Number_replicated_hits_by_dataset.png")
png(output_path, width = 4000, height = 4500, res = 300)
print(p)
dev.off()

# By NMD Method
p <- RGWAS_replication_hits_final %>%
        filter(FDR_threshold_used == 2 & randomization == "no") %>%
        group_by(NMD_method) %>%
        summarise(replicated_hits = n()) %>%
        ggplot(aes(x = NMD_method, y = replicated_hits)) + 
            geom_bar(stat='identity',position=position_dodge()) +
            theme_bw(base_size = 25) +
            theme(axis.text.x = element_text(size = 18, angle = 45, hjust=1))
output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/rare_germline_variants/replication/",genelist,"/between/Number_replicated_hits_by_NMD_method.png")
png(output_path, width = 4000, height = 4500, res = 300)
print(p)
dev.off()

# By tissue
p <- RGWAS_replication_hits_final %>%
        filter(randomization == "no") %>%
        group_by(matched_tissue,FDR_threshold_used,NMD_method) %>%
        summarise(replicated_hits = n()) %>%
        arrange(replicated_hits) %>%
        #mutate(matched_tissue = factor(matched_tissue, levels = rev(matched_tissue))) %>%
            ggplot(aes(x = matched_tissue, y = replicated_hits, fill = factor(FDR_threshold_used))) + 
                theme_bw() + scale_fill_brewer(palette = "OrRd", direction = -1) +
                labs(title = "", x = "", y = "Replicated hits", fill = "% FDR threshold") + 
                geom_bar(stat='identity',position=position_dodge()) +
                theme_bw(base_size = 45) +
                facet_grid(. ~ NMD_method) +
                theme(axis.text.x = element_text(size = 38, angle = 45, hjust=1),
                axis.text.y = element_text(size = 35),
                legend.position = "top", legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
                legend.title = element_text(size = 34)) +
                guides(fill = guide_legend(override.aes = list(size = 8)))

output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/rare_germline_variants/replication/",genelist,"/between/Number_replicated_hits_by_tissue.png")
png(output_path, width = 5000, height = 4500, res = 300)
print(p)
dev.off()

# By gene
# By tissue
p <- RGWAS_replication_hits_final %>%
        filter(FDR_threshold_used == 2 & randomization == "no") %>%
        group_by(variant) %>%
        summarise(replicated_hits = n()) %>%
        arrange(replicated_hits) %>%
        mutate(variant = factor(variant, levels = rev(variant))) %>%
        ggplot(aes(x = variant, y = replicated_hits)) + 
            geom_bar(stat='identity',position=position_dodge()) +
            theme_bw(base_size = 25) +
            #facet_grid(. ~ NMD_method_validation) +
            theme(axis.text.x = element_text(size = 18, angle = 45, hjust=1))
output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/rare_germline_variants/replication/",genelist,"/between/Number_replicated_hits_by_gene.png")
png(output_path, width = 4000, height = 4500, res = 300)
print(p)
dev.off()

# 2.4) Genes table
RGWAS_replication_hits_final$gene_tissue <- paste0(RGWAS_replication_hits_final$variant," - ",RGWAS_replication_hits_final$matched_tissue)
RGWAS_replication_hits_final$database_translation <- paste0(RGWAS_replication_hits_final$database_discovery,"/",RGWAS_replication_hits_final$database_validation)

RGWAS_replication_hits_table <- RGWAS_replication_hits_final %>%
            filter(FDR_threshold_used == 2 & randomization == "no") %>%
            mutate(facet_order = factor(paste0(database, " - ",disc_val),
                            levels = c("GTEx - discovery", "TCGA - validation",
                                        "TCGA - discovery", "GTEx - validation")))

# Manual change
RGWAS_replication_hits_table <- RGWAS_replication_hits_table %>%
                mutate(n_carriers = ifelse(variant == "LAMC1" & NMD_method == "Endogenous" & dataset == "PTV_Missense_CADD15_0.1perc", 4, n_carriers)) %>%
                mutate(n_carriers = ifelse(variant == "CRTC1" & NMD_method == "Endogenous" & dataset == "PTV_Missense_CADD15_0.1perc", 3, n_carriers))
# Keep only GTEx as discovery and TCGA as validation

p <- RGWAS_replication_hits_table %>%
        filter(database_discovery == "GTEx" & database_validation == "TCGA") %>%
        group_by(disc_val) %>%
            ggplot(aes(x = dataset,y = gene_tissue)) +
            geom_tile(aes(fill = n_carriers), size = 2) +
            geom_text(aes(label = n_carriers),color = "black", size=10) +
            theme_bw(base_size = 35) +
            facet_nested(. ~ NMD_method + facet_order) +
            theme(axis.text.x = element_text(size = 35, angle = 45, hjust=1),
                    axis.text.y = element_text(size = 35),
                    legend.title=element_text(size=35)) +
            scale_fill_gradientn(colours = c('grey','#F21A00')) +
            labs(fill = '# of individuals') + xlab("") + ylab("") +
            guides(fill = guide_legend(override.aes = list(size = 20)))
output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/rare_germline_variants/replication/",genelist,"/between/table_genes_replicated.png")
png(output_path, width = 11000, height = 6500, res = 300)
print(p)
dev.off()

### Supplementary
RGWAS_replication_hits_final$gene_tissue <- paste0(RGWAS_replication_hits_final$variant," - ",RGWAS_replication_hits_final$matched_tissue)
RGWAS_replication_hits_final$database_translation <- paste0(RGWAS_replication_hits_final$database_discovery,"/",RGWAS_replication_hits_final$database_validation)

RGWAS_replication_hits_table <- RGWAS_replication_hits_final %>%
            filter(FDR_threshold_used == 5 & randomization == "no") %>%
            mutate(facet_order = factor(paste0(database, " - ",disc_val),
                            levels = c("GTEx - discovery", "TCGA - validation",
                                        "TCGA - discovery", "GTEx - validation")))
# Manual change
RGWAS_replication_hits_table <- RGWAS_replication_hits_table %>%
                mutate(n_carriers = ifelse(variant == "LAMC1" & NMD_method == "Endogenous" & dataset == "PTV_Missense_CADD15_0.1perc", 4, n_carriers)) %>%
                mutate(n_carriers = ifelse(variant == "CRTC1" & NMD_method == "Endogenous" & dataset == "PTV_Missense_CADD15_0.1perc", 3, n_carriers))

RGWAS_replication_hits_table$dataset <- gsub("PTV","pLoF",RGWAS_replication_hits_table$dataset)
RGWAS_replication_hits_table$dataset <- gsub("_0.1perc","",RGWAS_replication_hits_table$dataset)

write.table(RGWAS_replication_hits_table, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig16/SuppFig16B.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
saveRDS(RGWAS_replication_hits_table, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig16/SuppFig16B.RData")


