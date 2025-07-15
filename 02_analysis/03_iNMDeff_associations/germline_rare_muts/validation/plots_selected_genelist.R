# No conda
library("GWASTools")
library("VennDiagram")
#library("ggVennDiagram")
library("dplyr")
library("cowplot")
library("ggplot2")
library("RColorBrewer")

sig_hits_within_databases_old <- function(RGWAS_df_all, discovery, validation, empirical_FDR) {
    empirical_FDR_char <- paste0(empirical_FDR,"%")
    empirical_FDR <- as.numeric(empirical_FDR/100)
    # Empirical FDR in discovery
    RGWAS_empirical_FDR <- RGWAS_df_all %>%
            filter(randomization == "yes") %>%
            group_by(database, tissue, NMD_method, random_geneset) %>%
            summarise(empirical_FDR = as.numeric(quantile(pValue, seq(0,1,0.01))[empirical_FDR_char]))
    # Sig hits in discovery
    RGWAS_discovery_sig <- RGWAS_df_all %>%
                    filter(randomization == "no" & random_geneset == "no") %>%
                    group_by(database, tissue, random_geneset, NMD_method) %>%
                    mutate(empirical_FDR = as.numeric(RGWAS_empirical_FDR[RGWAS_empirical_FDR$NMD_method %in% NMD_method & 
                    RGWAS_empirical_FDR$tissue %in% tissue &
                    RGWAS_empirical_FDR$database %in% database &
                    RGWAS_empirical_FDR$random_geneset %in% random_geneset, "empirical_FDR"])) %>%
                    mutate(sig_genes_ID = ifelse( ((pValue < empirical_FDR) & (NMD_method == discovery)), paste0(database,"_",variant,"_",tissue,"_",random_geneset), NA))
    sig_genes_ID_char <- as.character(na.omit(RGWAS_discovery_sig$sig_genes_ID))
    # Empirical FDR in validation (only for the subset of genes that are significant in discovery)
    RGWAS_empirical_FDR_validation <- RGWAS_discovery_sig %>%
        group_by(database,tissue,random_geneset) %>%
        filter(NMD_method == validation & paste0(database,"_",variant,"_",tissue,"_",random_geneset) %in% sig_genes_ID_char) %>%
        mutate(empirical_FDR = quantile(pValue, seq(0,1,0.01))[empirical_FDR_char]) %>%
        filter(pValue < empirical_FDR)
    RGWAS_validation_sig <- data.frame(RGWAS_empirical_FDR_validation[RGWAS_empirical_FDR_validation$empirical_FDR < empirical_FDR,])
    return(RGWAS_validation_sig)
}


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

replicated_hits <- function(discovery, validation, FDR, geneVariantsThres) {
    # Variant-gene threshold
    discovery <- discovery[discovery$ASE_n_carriers >= geneVariantsThres,]
    print(paste0("N from discovery --> ",nrow(discovery)))
    # Variant-gene threshold
    validation <- validation[validation$END_n_carriers >= geneVariantsThres,]
    print(paste0("N from validation --> ",nrow(validation)))
    # Significant hits in discovery
    discovery$ASE_pValueAdjust <- p.adjust(discovery$ASE_pValue, method = "fdr")
    discovery_filt <- discovery[discovery$ASE_pValueAdjust < FDR,]
    genes_overlap <- intersect(as.character(discovery_filt$ASE_variant), validation$END_variant) 
    # Replicated hits in validation (correct FDR only in the significants from discovery)
    validation_filt <- validation[validation$END_variant %in% genes_overlap,]
    validation_filt$END_pValueAdjust <- p.adjust(validation_filt$END_pValue, method = "fdr")
    validation_filt <- validation_filt[validation_filt$END_pValueAdjust < FDR,]
    replicated_genes <- as.character(validation_filt$END_variant)
    # Merge
    discovery_replicated <- discovery[discovery$ASE_variant %in% replicated_genes,]
    validation_replicated <- validation[validation$END_variant %in% replicated_genes,]
    final_df <- merge(discovery_replicated,validation_replicated, by.x = "ASE_variant", by.y = "END_variant", all.x = TRUE)
    return(final_df)
}

matched_tissue_RGWAS_replication <- function(NMD_method_discovery, NMD_method_validation, random_geneset_char,
                                                empirical_FDR, database_discovery, database_validation) {
    empirical_FDR_char <- paste0(empirical_FDR,"%")
    #empirical_FDR_num <- as.numeric(empirical_FDR/100)
    RGWAS_replicated_final_df <- c()
    # Empirical FDR in discovery
    RGWAS_empirical_FDR <- RGWAS_df_all %>%
            filter(randomization == "yes") %>%
            group_by(database, tissue, NMD_method, dataset, random_geneset) %>%
            summarise(empirical_FDR = as.numeric(quantile(pValue, seq(0,1,0.01))[empirical_FDR_char]))
    # 1) Sig hits in discovery
    RGWAS_discovery_sig <- RGWAS_df_all %>%
                    filter(randomization == "no" & random_geneset == random_geneset_char) %>%
                    group_by(database, tissue, dataset, random_geneset, NMD_method) %>%
                    mutate(empirical_FDR = as.numeric(RGWAS_empirical_FDR[RGWAS_empirical_FDR$NMD_method %in% NMD_method & 
                    RGWAS_empirical_FDR$tissue %in% tissue & RGWAS_empirical_FDR$dataset %in% dataset & 
                    RGWAS_empirical_FDR$database %in% database &
                    RGWAS_empirical_FDR$random_geneset %in% random_geneset, "empirical_FDR"])) %>%
                    mutate(sig_genes_ID = ifelse( ((pValue < empirical_FDR) & (NMD_method == NMD_method_discovery)), paste0(variant,"_",tissue,"_",dataset,"_",random_geneset), NA))
    sig_genes_ID_char <- as.character(na.omit(RGWAS_discovery_sig$sig_genes_ID))
    # Number of genes tested at discovery
    num_genes_tested <- RGWAS_df_all %>%
                        filter(randomization == "no" & random_geneset == random_geneset_char & database == database_discovery) %>%
                        summarise(genes_tested = n()) %>% data.frame()

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
        # Empirical FDR in validation (only for the subset of genes that are significant in discovery)
        RGWAS_empirical_FDR_validation <- RGWAS_discovery_sig %>%
            filter(database == database_validation) %>%
            group_by(tissue,dataset) %>%
            filter(NMD_method == NMD_method_validation & 
                paste0(variant,"_",tissue,"_",dataset,"_",random_geneset) %in% sig_genes_ID_char_tissue) %>%
            mutate(empirical_FDR_validation = quantile(pValue, seq(0,1,0.01))[empirical_FDR_char]) %>%
            filter(pValue < empirical_FDR_validation)
        #RGWAS_validation_sig <- data.frame(RGWAS_empirical_FDR_validation[RGWAS_empirical_FDR_validation$empirical_FDR < empirical_FDR,])
        RGWAS_validation_sig <- RGWAS_empirical_FDR_validation
        if (nrow(RGWAS_validation_sig) == 0 ) {next}
        print(RGWAS_validation_sig)
        # Save
        #colnames(RGWAS_validation_sig)[colnames(RGWAS_validation_sig) == "tissue"] <- "TCGA_cancer"
        RGWAS_validation_sig$matched_tissue <- tissue
        RGWAS_validation_sig$NMD_method_discovery <- NMD_method_discovery
        RGWAS_validation_sig$NMD_method_validation <- NMD_method_validation
        RGWAS_validation_sig$database_discovery <- database_discovery
        RGWAS_validation_sig$database_validation <- database_validation
        RGWAS_validation_sig$genes_tested <- num_genes_tested$genes_tested
        RGWAS_replicated_final_df <- rbind(RGWAS_replicated_final_df,RGWAS_validation_sig)
    }
    return(RGWAS_replicated_final_df)
}

matched_tissue_RGWAS_replication_no_strat <- function(random_geneset_char, empirical_FDR, database_discovery, database_validation) {
    empirical_FDR_char <- paste0(empirical_FDR,"%")
    #empirical_FDR_num <- as.numeric(empirical_FDR/100)
    RGWAS_replicated_final_df <- c()
    # Empirical FDR in discovery
    RGWAS_empirical_FDR <- RGWAS_df_all %>%
            filter(randomization == "yes") %>%
            group_by(database, tissue, random_geneset) %>%
            summarise(empirical_FDR = as.numeric(quantile(pValue, seq(0,1,0.01))[empirical_FDR_char]))
    # 1) Sig hits in discovery
    RGWAS_discovery_sig <- RGWAS_df_all %>%
                    filter(randomization == "no" & random_geneset == random_geneset_char) %>%
                    group_by(database, tissue, random_geneset) %>%
                    mutate(empirical_FDR = as.numeric(RGWAS_empirical_FDR[RGWAS_empirical_FDR$tissue %in% tissue & RGWAS_empirical_FDR$database %in% database &
                    RGWAS_empirical_FDR$random_geneset %in% random_geneset, "empirical_FDR"])) %>%
                    mutate(sig_genes_ID = ifelse( ((pValue < empirical_FDR)), paste0(variant,"_",tissue,"_",random_geneset), NA))
    sig_genes_ID_char <- as.character(na.omit(RGWAS_discovery_sig$sig_genes_ID))
    # Number of genes tested at discovery
    num_genes_tested <- RGWAS_df_all %>%
                        filter(randomization == "no" & random_geneset == random_geneset_char & database == database_discovery) %>%
                        summarise(genes_tested = n()) %>% data.frame()

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
        # Empirical FDR in validation (only for the subset of genes that are significant in discovery)
        RGWAS_empirical_FDR_validation <- RGWAS_discovery_sig %>%
            filter(database == database_validation) %>%
            group_by(tissue) %>%
            filter(paste0(variant,"_",tissue,"_",random_geneset) %in% sig_genes_ID_char_tissue) %>%
            mutate(empirical_FDR_validation = quantile(pValue, seq(0,1,0.01))[empirical_FDR_char]) %>%
            filter(pValue < empirical_FDR_validation)
        #RGWAS_validation_sig <- data.frame(RGWAS_empirical_FDR_validation[RGWAS_empirical_FDR_validation$empirical_FDR < empirical_FDR,])
        RGWAS_validation_sig <- RGWAS_empirical_FDR_validation
        if (nrow(RGWAS_validation_sig) == 0 ) {next}
        print(RGWAS_validation_sig)
        # Save
        #colnames(RGWAS_validation_sig)[colnames(RGWAS_validation_sig) == "tissue"] <- "TCGA_cancer"
        RGWAS_validation_sig$matched_tissue <- tissue
        RGWAS_validation_sig$database_discovery <- database_discovery
        RGWAS_validation_sig$database_validation <- database_validation
        RGWAS_validation_sig$genes_tested <- num_genes_tested$genes_tested
        RGWAS_replicated_final_df <- rbind(RGWAS_replicated_final_df,RGWAS_validation_sig)
    }
    return(RGWAS_replicated_final_df)
}

matched_tissue_RGWAS_replication_newEmpiricalFDR <- function(NMD_method_discovery, NMD_method_validation, random_geneset_char,
                                                empirical_FDR, database_discovery, database_validation) {
    empirical_FDR_char <- paste0(empirical_FDR,"%")
    RGWAS_replicated_final_df <- c()

    # hits in random
    RGWAS_FDR_hits_randomization <- RGWAS_df_all %>%
            filter(randomization == "yes") %>%
            group_by(database, NMD_method, dataset, tissue, random_geneset) %>%
            mutate(empirical_FDR = as.numeric(quantile(pValue, seq(0,1,0.01))[empirical_FDR_char])) %>%
            summarise(num_hits_randomization = sum(pValue < empirical_FDR))
    # Number of genes tested
    num_genes_tested <- RGWAS_df_all %>%
                        filter(randomization == "yes") %>%
                        group_by(database,NMD_method, dataset, tissue, random_geneset) %>%
                        summarise(genes_tested = n()) %>% data.frame()
    RGWAS_FDR_hits_randomization <- left_join(RGWAS_FDR_hits_randomization,num_genes_tested)
    RGWAS_FDR_hits_randomization$prop_num_hits_randomization <- RGWAS_FDR_hits_randomization$num_hits_randomization / RGWAS_FDR_hits_randomization$genes_tested
    RGWAS_FDR_hits_randomization$genes_tested <- NULL
    # hits in non random
    RGWAS_hits_non_randomization <- RGWAS_df_all %>%
            filter(randomization == "no") %>%
            group_by(database, dataset, tissue, random_geneset) %>%
            mutate(empirical_FDR = as.numeric(quantile(pValue, seq(0,1,0.01))[empirical_FDR_char])) %>%
            summarise(num_hits_non_randomization = sum(pValue < empirical_FDR))
    # Number of genes tested
    num_genes_tested <- RGWAS_df_all %>%
                        filter(randomization == "no") %>%
                        group_by(database,tissue, dataset, NMD_method,random_geneset) %>%
                        summarise(genes_tested = n()) %>% data.frame()
    RGWAS_hits_non_randomization <- left_join(RGWAS_hits_non_randomization,num_genes_tested)
    RGWAS_hits_non_randomization$prop_num_hits_non_randomization <- RGWAS_hits_non_randomization$num_hits_non_randomization / RGWAS_hits_non_randomization$genes_tested
    RGWAS_hits_non_randomization$genes_tested <- NULL
    # FDR empirical
    RGWAS_empirical_FDR <- left_join(RGWAS_hits_non_randomization,RGWAS_FDR_hits_randomization)
    RGWAS_empirical_FDR$empirical_FDR <- RGWAS_empirical_FDR$prop_num_hits_randomization / RGWAS_empirical_FDR$prop_num_hits_non_randomization

    # 1) Sig hits in discovery
    RGWAS_discovery_sig <- RGWAS_df_all %>%
                    filter(randomization == "no" & random_geneset == random_geneset_char) %>%
                    group_by(database, tissue, dataset, random_geneset, NMD_method) %>%
                    mutate(empirical_FDR = as.numeric(RGWAS_empirical_FDR[RGWAS_empirical_FDR$NMD_method %in% NMD_method & 
                    RGWAS_empirical_FDR$tissue %in% tissue & RGWAS_empirical_FDR$dataset %in% dataset & 
                    RGWAS_empirical_FDR$database %in% database &
                    RGWAS_empirical_FDR$random_geneset %in% random_geneset, "empirical_FDR"])) %>%
                    mutate(sig_genes_ID = ifelse( ((pValue < empirical_FDR) & (NMD_method == NMD_method_discovery)), paste0(variant,"_",tissue,"_",dataset,"_",random_geneset), NA))
    sig_genes_ID_char <- as.character(na.omit(RGWAS_discovery_sig$sig_genes_ID))

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
        # Empirical FDR in validation (only for the subset of genes that are significant in discovery)
        RGWAS_empirical_FDR_validation <- RGWAS_discovery_sig %>%
            filter(database == database_validation) %>%
            group_by(tissue,dataset) %>%
            filter(NMD_method == NMD_method_validation & 
                paste0(variant,"_",tissue,"_",dataset,"_",random_geneset) %in% sig_genes_ID_char_tissue) %>%
            mutate(empirical_FDR_validation = quantile(pValue, seq(0,1,0.01))[empirical_FDR_char]) %>%
            filter(pValue < empirical_FDR_validation)
        #RGWAS_validation_sig <- data.frame(RGWAS_empirical_FDR_validation[RGWAS_empirical_FDR_validation$empirical_FDR < empirical_FDR,])
        RGWAS_validation_sig <- RGWAS_empirical_FDR_validation
        if (nrow(RGWAS_validation_sig) == 0 ) {next}
        print(RGWAS_validation_sig)
        # Save
        #colnames(RGWAS_validation_sig)[colnames(RGWAS_validation_sig) == "tissue"] <- "TCGA_cancer"
        RGWAS_validation_sig$matched_tissue <- tissue
        RGWAS_validation_sig$NMD_method_discovery <- NMD_method_discovery
        RGWAS_validation_sig$NMD_method_validation <- NMD_method_validation
        RGWAS_validation_sig$database_discovery <- database_discovery
        RGWAS_validation_sig$database_validation <- database_validation
        # Number of genes tested at discovery
        num_genes_tested <- RGWAS_df_all %>%
                            filter(randomization == "no" & random_geneset == random_geneset_char & database == database_discovery) %>%
                        summarise(genes_tested = n()) %>% data.frame()
        RGWAS_validation_sig$genes_tested <- num_genes_tested$genes_tested
        RGWAS_replicated_final_df <- rbind(RGWAS_replicated_final_df,RGWAS_validation_sig)
    }
    return(RGWAS_replicated_final_df)
}

matched_tissue_RGWAS_replication_newEmpiricalFDR_no_strat <- function(random_geneset_char, empirical_FDR, database_discovery, database_validation) {
    empirical_FDR_char <- paste0(empirical_FDR,"%")
    RGWAS_replicated_final_df <- c()
    # hits in random
    RGWAS_FDR_hits_randomization <- RGWAS_df_all %>%
            filter(randomization == "yes") %>%
            group_by(database, tissue, NMD_method, random_geneset) %>%
            mutate(empirical_FDR = as.numeric(quantile(pValue, seq(0,1,0.01))[empirical_FDR_char])) %>%
            summarise(num_hits_randomization = sum(pValue < empirical_FDR))
    # Number of genes tested
    num_genes_tested <- RGWAS_df_all %>%
                        filter(randomization == "yes") %>%
                        group_by(database, tissue, NMD_method,random_geneset) %>%
                        summarise(genes_tested = n()) %>% data.frame()
    RGWAS_FDR_hits_randomization <- left_join(RGWAS_FDR_hits_randomization,num_genes_tested)
    RGWAS_FDR_hits_randomization$prop_num_hits_randomization <- RGWAS_FDR_hits_randomization$num_hits_randomization / RGWAS_FDR_hits_randomization$genes_tested
    RGWAS_FDR_hits_randomization$genes_tested <- NULL
    # hits in non random
    RGWAS_hits_non_randomization <- RGWAS_df_all %>%
            filter(randomization == "no") %>%
            group_by(database, tissue, NMD_method, random_geneset) %>%
            mutate(empirical_FDR = as.numeric(quantile(pValue, seq(0,1,0.01))[empirical_FDR_char])) %>%
            summarise(num_hits_non_randomization = sum(pValue < empirical_FDR))
    # Number of genes tested
    num_genes_tested <- RGWAS_df_all %>%
                        filter(randomization == "no") %>%
                        group_by(database,tissue,NMD_method,random_geneset) %>%
                        summarise(genes_tested = n()) %>% data.frame()
    RGWAS_hits_non_randomization <- left_join(RGWAS_hits_non_randomization,num_genes_tested)
    RGWAS_hits_non_randomization$prop_num_hits_non_randomization <- RGWAS_hits_non_randomization$num_hits_non_randomization / RGWAS_hits_non_randomization$genes_tested
    RGWAS_hits_non_randomization$genes_tested <- NULL
    # FDR empirical
    RGWAS_empirical_FDR <- left_join(RGWAS_hits_non_randomization,RGWAS_FDR_hits_randomization)
    RGWAS_empirical_FDR$empirical_FDR <- RGWAS_empirical_FDR$prop_num_hits_randomization / RGWAS_empirical_FDR$prop_num_hits_non_randomization


    # df <- RGWAS_df_all[RGWAS_df_all$tissue == "Artery_Tibial" & RGWAS_df_all$random_geneset == "no" & RGWAS_df_all$randomization == "no" &
    #                     RGWAS_df_all$NMD_method == "Endogenous",]
    # df$empirical_FDR <- as.numeric(quantile(df$pValue, seq(0,1,0.01))[empirical_FDR_char])
    # df <- df[df$pValue < df$empirical_FDR,]
    # dim(df)

    # 1) Sig hits in discovery
    RGWAS_discovery_sig <- RGWAS_df_all %>%
                    filter(randomization == "no" & random_geneset == random_geneset_char) %>%
                    group_by(database, tissue,NMD_method, random_geneset) %>%
                    mutate(empirical_FDR = as.numeric(RGWAS_empirical_FDR[RGWAS_empirical_FDR$tissue %in% tissue & RGWAS_empirical_FDR$database %in% database & 
                                RGWAS_empirical_FDR$random_geneset %in% random_geneset &
                                RGWAS_empirical_FDR$NMD_method %in% NMD_method, "empirical_FDR"])) %>%
                    mutate(sig_genes_ID = ifelse( ((pValue < empirical_FDR) ), paste0(variant,"_",tissue,"_",NMD_method,"_",random_geneset), NA))
    sig_genes_ID_char <- as.character(na.omit(RGWAS_discovery_sig$sig_genes_ID))

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
        # Empirical FDR in validation (only for the subset of genes that are significant in discovery)
        RGWAS_empirical_FDR_validation <- RGWAS_discovery_sig %>%
            filter(database == database_validation) %>%
            group_by(tissue,NMD_method) %>%
            filter(paste0(variant,"_",tissue,"_",NMD_method,"_",random_geneset) %in% sig_genes_ID_char_tissue) %>%
            mutate(empirical_FDR_validation = quantile(pValue, seq(0,1,0.01))[empirical_FDR_char]) %>%
            filter(pValue < empirical_FDR_validation)
        #RGWAS_validation_sig <- data.frame(RGWAS_empirical_FDR_validation[RGWAS_empirical_FDR_validation$empirical_FDR < empirical_FDR,])
        RGWAS_validation_sig <- RGWAS_empirical_FDR_validation
        if (nrow(RGWAS_validation_sig) == 0 ) {next}
        print(RGWAS_validation_sig)
        # Save
        #colnames(RGWAS_validation_sig)[colnames(RGWAS_validation_sig) == "tissue"] <- "TCGA_cancer"
        RGWAS_validation_sig$matched_tissue <- tissue
        RGWAS_validation_sig$database_discovery <- database_discovery
        RGWAS_validation_sig$database_validation <- database_validation
        # Number of genes tested at discovery
        num_genes_tested <- RGWAS_df_all %>%
                            filter(randomization == "no" & random_geneset == random_geneset_char & database == database_discovery) %>%
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


matched_tissue_RGWAS_replication_normal_FDR_no_strat <- function(random_geneset_char, FDR, database_discovery, database_validation) {
    FDR <- FDR / 100
    # randomization_char <- randomization
    RGWAS_replicated_final_df <- c()
    # 1) Sig hits in discovery
    RGWAS_discovery_sig <- RGWAS_df_all %>%
                    filter(random_geneset == random_geneset_char) %>%
                    group_by(database, tissue, NMD_method) %>%
                    mutate(pvalue_FDR_adjusted = p.adjust(pValue, method = "fdr")) %>%
                    mutate(sig_genes_ID = ifelse( ((pvalue_FDR_adjusted < FDR) ), paste0(variant,"_",tissue,"_",NMD_method), NA))
    sig_genes_ID_char <- as.character(na.omit(RGWAS_discovery_sig$sig_genes_ID))

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
            filter(database == database_validation & random_geneset == random_geneset_char) %>%
            group_by(tissue,NMD_method,random_geneset) %>%
            filter(paste0(variant,"_",tissue,"_",NMD_method) %in% sig_genes_ID_char_tissue) %>%
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
                            filter(random_geneset == random_geneset_char & database == database_discovery) %>%
                        summarise(genes_tested = n()) %>% data.frame()
        RGWAS_validation_sig$genes_tested <- num_genes_tested$genes_tested
        RGWAS_replicated_final_df <- rbind(RGWAS_replicated_final_df,RGWAS_validation_sig)
    }
    return(RGWAS_replicated_final_df)
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

GTEx_tissue_names_path <- "/g/strcombio/fsupek_cancer1/gpalou/GTEx/RNAseq/GTEx_tissues_clean.txt"
GTEx_tissues <- read.table(file = GTEx_tissue_names_path, stringsAsFactors = FALSE)$V1
TCGA_cancer_names_path <- "/g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/TCGA_projects_names.txt"
TCGA_cancers <- read.table(file = TCGA_cancer_names_path, stringsAsFactors = FALSE)$V1
TCGA_cancers <- c(TCGA_cancers,"pancancer")
NMD_genes <- read.table( file = "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/NMD_genes.txt",
                        sep = "\t", header = TRUE)
output_path <- paste0("/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/TCGA_GTEx_match.txt")
TCGA_GTEx_tissues <- read.table(file = output_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
acronyms <- c("ADPSBQ","ADPVSC","ADRNLG","ARTAORT","ARTCRN","ARTTBL","BLDDER","BRNAMY","BRNACC","BRNCDT","BRNCHB","BRNCHA","BRNCTXA","BRNCTXB",
              "BRNHPP","BRNHPT","BRNNCC","BRNPTM","BRNSPC","BRNSNG","BREAST","FIBRBLS","LCL","CML","CVXECT","CVSEND","CLNSGM","CLNTRN","ESPGEJ","ESPMCS","ESPMSL",
              "FLLPNT","HRTAA","HRTLV","KDNCTX","KDNMDL","LIVER","LUNG","SLVRYG","MSCLSK","NERVET","OVARY","PNCREAS","PTTARY","PRSTTE","SKINNS","SKINS",
              "SNTTRM","SPLEEN","STMACH","TESTIS","THYROID","UTERUS","VAGINA","WHLBLD")
tissues_acronyms <- data.frame(GTEx_tissues_full = GTEx_tissues, acronyms = acronyms )
GTEx_tissues <- c(GTEx_tissues,"pantissue")
TCGA_GTEx_tissues[19,] <- c("pantissue","pantissue","pancancer")

# 1) Obtain datasets
geneList <- "yes"
if (geneList == "yes") {
    genelist <- "withgenelist"
} else if (geneList == "no") {
    genelist <- "nogenelist"
}
databases <- c("GTEx","TCGA")

# GTEx
GTEx_RGWAS_df_all <- c()
for (randomization in c("yes","no")) {
    for (randomGenes in c("yes","no")) {
        for (dataset in c("PTV_Missense_CADD15_0.1perc","PTV_Missense_CADD25_0.1perc","PTV_0.1perc")) {
            print(dataset)
            RWGAS_ASE <- RGWAS_df(dataset_name = dataset, NMD_method = "ASE", geneList = geneList, randomGenes = randomGenes,
                        randomization = randomization, geneVariantsThres = 2, VAF = 0.2, tissues = GTEx_tissues, database = "GTEx")
            RWGAS_ASE$NMD_method <- "ASE"
            RWGAS_end <- RGWAS_df(dataset_name = dataset, NMD_method = "endogenous", geneList = geneList, randomGenes = randomGenes,
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
#table(duplicated(GTEx_RGWAS_df_all))

# TCGA
TCGA_RGWAS_df_all <- c()
for (randomization in c("yes","no")) {
    for (randomGenes in c("yes","no")) {
        for (dataset in c("PTV_Missense_CADD15_0.1perc","PTV_Missense_CADD25_0.1perc","PTV_0.1perc","MRT_25thCentile_0.1perc","CCR_90orhigher_0.1perc")) {
            print(dataset)
            RWGAS_ASE <- RGWAS_df(dataset_name = dataset, NMD_method = "ASE", geneList = geneList, randomGenes = randomGenes,
                        randomization = randomization, geneVariantsThres = 2, VAF = 0.2, tissues = TCGA_cancers, database = "TCGA")
            RWGAS_ASE$NMD_method <- "ASE"
            RWGAS_end <- RGWAS_df(dataset_name = dataset, NMD_method = "endogenous", geneList = geneList, randomGenes = randomGenes,
                        randomization = randomization, geneVariantsThres = 2, VAF = 0.2, tissues = TCGA_cancers, database = "TCGA")
            RWGAS_end$NMD_method <- "Endogenous"
            if (is.null(nrow(RWGAS_end))) {
                RWGAS_dataset <- RWGAS_ASE
            } else {
                RWGAS_dataset <- rbind(RWGAS_ASE,RWGAS_end)
            }  
            RWGAS_dataset$dataset <- dataset
            RWGAS_dataset$random_geneset <- randomGenes
            RWGAS_dataset$randomization <- randomization
            if (length(TCGA_RGWAS_df_all) == 0) {
                TCGA_RGWAS_df_all <- RWGAS_dataset
            } else {
                TCGA_RGWAS_df_all <- rbind(TCGA_RGWAS_df_all,RWGAS_dataset)
            }
        }
    }
}
print(dim(TCGA_RGWAS_df_all))
#table(duplicated(TCGA_RGWAS_df_all))
colnames(TCGA_RGWAS_df_all)[colnames(TCGA_RGWAS_df_all) %in% "cancer_type"] <- "tissue"
# Merge
TCGA_RGWAS_df_all$database <- "TCGA"
GTEx_RGWAS_df_all$database <- "GTEx"
RGWAS_df_all <- rbind(GTEx_RGWAS_df_all,TCGA_RGWAS_df_all)
# output_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/rare_germline_variants/replication/genes_list.txt"
# write.table(data.frame(a), file = output_path, sep = "\t", quote = FALSE,
#                       col.names = FALSE, row.names = FALSE)

############ FDR correction for each analysis ############
# From Mischan paper -->
# We calculated empirical FDR thresholds for each cancer type (or pan-cancer) separately. For instance, 
# the p-value at which 1% of the associations from the randomized run would have been called as a hit 
# (false discovery) corresponds to a FDR of 1%.

# Based on the conservative hypothesis that there would be no real associations from the random list of genes, 
# we calculated FDRs at different empirical FDR thresholds by dividing the number of hits, which were detected 
# via the random list of genes by the number of genes detected at the same empirical FDR with our pre-selected list of genes. 
# For instance, at an empirical FDR of 1% we identified 44 hits with our random list of genes and 207 
# hits with our pre-selected list of genes. Thus, we estimated a FDR of 44/207 ≈ 21% at our empirical FDR of 1%.

# 2) Results

# 2.1) Number of significant hits
# FDR empirical correction splitting by NMD_method/dataset/database/geneset


#############

# hits in random
RGWAS_FDR_hits_randomization <- RGWAS_df_all %>%
        filter(randomization == "yes") %>%
        group_by(database, tissue, random_geneset) %>%
        mutate(empirical_FDR = as.numeric(quantile(pValue, seq(0,1,0.01))["5%"])) %>%
        summarise(num_hits_randomization = sum(pValue < empirical_FDR))
# Number of genes tested
num_genes_tested <- RGWAS_df_all %>%
                    filter(randomization == "yes") %>%
                    group_by(database, tissue, random_geneset) %>%
                    summarise(genes_tested = n()) %>% data.frame()
RGWAS_FDR_hits_randomization <- left_join(RGWAS_FDR_hits_randomization,num_genes_tested)
RGWAS_FDR_hits_randomization$prop_num_hits_randomization <- RGWAS_FDR_hits_randomization$num_hits_randomization / RGWAS_FDR_hits_randomization$genes_tested
RGWAS_FDR_hits_randomization$genes_tested <- NULL
# hits in non random
RGWAS_hits_non_randomization <- RGWAS_df_all %>%
        filter(randomization == "no") %>%
        group_by(database, tissue, random_geneset) %>%
        mutate(empirical_FDR = as.numeric(quantile(pValue, seq(0,1,0.01))["5%"])) %>%
        summarise(num_hits_non_randomization = sum(pValue < empirical_FDR))
# Number of genes tested
num_genes_tested <- RGWAS_df_all %>%
                    filter(randomization == "no") %>%
                    group_by(database,tissue,random_geneset) %>%
                    summarise(genes_tested = n()) %>% data.frame()
RGWAS_hits_non_randomization <- left_join(RGWAS_hits_non_randomization,num_genes_tested)
RGWAS_hits_non_randomization$prop_num_hits_non_randomization <- RGWAS_hits_non_randomization$num_hits_non_randomization / RGWAS_hits_non_randomization$genes_tested
RGWAS_hits_non_randomization$genes_tested <- NULL
# FDR empirical
RGWAS_empirical_FDR <- left_join(RGWAS_hits_non_randomization,RGWAS_FDR_hits_randomization)

RGWAS_empirical_FDR$empirical_FDR <- RGWAS_empirical_FDR$num_hits_randomization / RGWAS_empirical_FDR$num_hits_non_randomization
RGWAS_empirical_FDR$empirical_FDR_prop <- RGWAS_empirical_FDR$prop_num_hits_randomization / RGWAS_empirical_FDR$prop_num_hits_non_randomization


RGWAS_empirical_FDR$empirical_FDR_prop
RGWAS_empirical_FDR <- RGWAS_empirical_FDR[order(RGWAS_empirical_FDR$empirical_FDR_prop),]

df <- RGWAS_empirical_FDR %>%
            filter(random_geneset == "no")
df <- df[order(df$empirical_FDR_prop),]

head(data.frame(df))      

# ##############

RGWAS_empirical_FDR <- RGWAS_df_all %>%
        filter(randomization == "yes") %>%
        group_by(database, tissue, NMD_method, dataset, random_geneset) %>%
        summarise(empirical_FDR = as.numeric(quantile(pValue, seq(0,1,0.01))["5%"]))

df <- RGWAS_df_all %>%
        filter(randomization == "no") %>%
        group_by(database,tissue, NMD_method, dataset, random_geneset) %>%
        mutate(empirical_FDR = as.numeric(RGWAS_empirical_FDR[RGWAS_empirical_FDR$tissue %in% tissue & 
            RGWAS_empirical_FDR$NMD_method %in% NMD_method & RGWAS_empirical_FDR$dataset %in% dataset & 
            RGWAS_empirical_FDR$random_geneset %in% random_geneset & 
            RGWAS_empirical_FDR$database %in% database, "empirical_FDR"])) %>% 
        summarise(num_hits = sum(pValue < empirical_FDR)) %>% data.frame()

num_hits <- na.omit(df) %>%
    group_by(database,NMD_method,dataset,random_geneset) %>%
    summarise(sig_hits = sum(num_hits)) %>% data.frame()

# Number of genes tested
num_genes_tested <- RGWAS_df_all %>%
                    filter(randomization == "no") %>%
                    group_by(database, NMD_method,dataset,random_geneset) %>%
                    summarise(genes_tested = n()) %>% data.frame()
num_hits_total <- merge(num_hits, num_genes_tested)
num_hits_total$prop_sig_hits <- round(num_hits_total$sig_hits / num_hits_total$genes_tested,3)

p <- ggplot(num_hits_total, aes(x= dataset, y = prop_sig_hits, fill = factor(random_geneset))) + 
        geom_bar(stat='identity',position=position_dodge()) +
        ylim(c(0,0.25)) +
        facet_grid(database ~ NMD_method) +
        theme_classic(base_size = 25) +
        theme(axis.text.x = element_text(size = 14, angle = 45, hjust=1))
output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/rare_germline_variants/replication/",genelist,"/Number_sig_hits_full_split.png")
png(output_path, width = 4500, height = 3500, res = 300)
print(p)
dev.off()

# FDR empirical correction without splitting by NMD_method/dataset
RGWAS_empirical_FDR <- RGWAS_df_all %>%
        filter(randomization == "yes" & !dataset %in% c("CCR_90orhigher_0.1perc","MRT_25thCentile_0.1perc")) %>%
        group_by(database,tissue, random_geneset) %>%
        summarise(empirical_FDR = as.numeric(quantile(pValue, seq(0,1,0.01))["5%"]))

df <- RGWAS_df_all %>%
        filter(randomization == "no" & !dataset %in% c("CCR_90orhigher_0.1perc","MRT_25thCentile_0.1perc")) %>%
        group_by(database, tissue, random_geneset) %>%
        mutate(empirical_FDR = as.numeric(RGWAS_empirical_FDR[RGWAS_empirical_FDR$tissue %in% tissue &  
            RGWAS_empirical_FDR$random_geneset %in% random_geneset & RGWAS_empirical_FDR$database %in% database, "empirical_FDR"])) %>% 
        summarise(num_hits = sum(pValue < empirical_FDR)) %>% data.frame()
num_hits_no_split <- df %>%
            group_by(database,random_geneset) %>%
            summarise(sig_hits = sum(num_hits)) %>% data.frame()
num_hits_split <- num_hits %>%
            filter(!dataset %in% c("CCR_90orhigher_0.1perc","MRT_25thCentile_0.1perc")) %>%
            group_by(database,random_geneset) %>%
            summarise(sig_hits = sum(sig_hits)) %>% data.frame()
num_hits_split$split <- "split"
num_hits_no_split$split <- "no_split"
num_hits_total <- rbind(num_hits_no_split,num_hits_split)
# Number of genes tested
num_genes_tested <- RGWAS_df_all %>%
                    filter(randomization == "no") %>%
                    group_by(database, random_geneset) %>%
                    summarise(genes_tested = n()) %>% data.frame()
num_hits_total <- merge(num_hits_total, num_genes_tested)
num_hits_total$prop_sig_hits <- round(num_hits_total$sig_hits / num_hits_total$genes_tested,3)

p <- ggplot(num_hits_total, aes(x= split, y = prop_sig_hits, fill = factor(random_geneset))) + 
        geom_bar(stat='identity',position=position_dodge()) +
        facet_wrap(~database) +
        theme_classic(base_size = 25) +
        theme(axis.text.x = element_text(size = 14, angle = 45, hjust=1))

output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/rare_germline_variants/replication/",genelist,"/Number_sig_hits_total_split_vs_no_split.png")
png(output_path, width = 4500, height = 3500, res = 300)
print(p)
dev.off()


# 2.2) Sig hits --> WITHIN databases and BETWEEN methods (any direction)
# Split by tissue

sig_hits_within_databases <- function(RGWAS_df_all, ddbb, empirical_FDR, random) {
    if (ddbb == "GTEx") {
        RGWAS_df_all <- RGWAS_df_all %>%
                    filter(tissue != "pantissue")
    }
    empirical_FDR_char <- paste0(empirical_FDR,"%")
    empirical_FDR_num <- as.numeric(empirical_FDR/100)
    # Empirical FDR in discovery
    RGWAS_empirical_FDR <- RGWAS_df_all %>%
            filter(randomization == "yes") %>%
            group_by(database, tissue, random_geneset) %>%
            summarise(empirical_FDR = as.numeric(quantile(pValue, seq(0,1,0.01))[empirical_FDR_char]))
    # Sig hits in ASE
    RGWAS_ASE_sig <- RGWAS_df_all %>%
                    filter(randomization == "no" & random_geneset == random & database == ddbb & NMD_method == "ASE") %>%
                    group_by(tissue) %>%
                    mutate(empirical_FDR = as.numeric(RGWAS_empirical_FDR[RGWAS_empirical_FDR$tissue %in% tissue &
                    RGWAS_empirical_FDR$random_geneset %in% random_geneset, "empirical_FDR"])) %>%
                    filter( pValue < empirical_FDR_num )
    # Sig hits in Endogenous
    RGWAS_END_sig <- RGWAS_df_all %>%
                    filter(randomization == "no" & random_geneset == random & database == ddbb & NMD_method == "Endogenous") %>%
                    group_by(tissue) %>%
                    mutate(empirical_FDR = as.numeric(RGWAS_empirical_FDR[RGWAS_empirical_FDR$tissue %in% tissue &
                    RGWAS_empirical_FDR$random_geneset %in% random_geneset, "empirical_FDR"])) %>%
                    filter( pValue < empirical_FDR_num )
    # Overlapping hits
    ASE_hits <- as.character(unique(RGWAS_ASE_sig$variant))
    END_hits <- as.character(unique(RGWAS_END_sig$variant))
    ASE_END_overlapping_hits <- intersect(ASE_hits,END_hits)
    # Total testable genes
    ASE_genes <- RGWAS_df_all %>%
                    filter(randomization == "no" & random_geneset == random & database == ddbb & NMD_method == "ASE") %>%
                    pull(variant)
    END_genes <- RGWAS_df_all %>%
                    filter(randomization == "no" & random_geneset == random & database == ddbb & NMD_method == "Endogenous") %>%
                    pull(variant)
    testable_genes <- ASE_genes
    vd_genes <- list(ASE_hits = ASE_hits, END_hits = END_hits, ASE_END_overlapping_hits = ASE_END_overlapping_hits, testable_genes = testable_genes )
    return(vd_genes)
}

# Venn Diagrams and Odds Ratio
myCol <- brewer.pal(2, "Pastel2")
myCol <- c("#B3E2CD","#CBD5E8")

for ( random_genes in c("yes","no") ) {
    if (random_genes == "yes") {
        genes <- "Random"
    } else{
        genes <- "Selected"
    }
    # GTEx
    vd_list_genes_GTEx <- sig_hits_within_databases(RGWAS_df_all, ddbb = "GTEx", empirical_FDR = 1, random = random_genes)
    OR_GTEx <- calculate_OR(vd_list_genes_GTEx)
    vd1 <- venn_diagram_within_databases(list = vd_list_genes_GTEx[1:2], title = paste0("GTEx \n OR: ",OR_GTEx))
    # TCGA
    vd_list_genes_TCGA <- sig_hits_within_databases(RGWAS_df_all, ddbb = "TCGA", empirical_FDR = 1, random = random_genes)
    OR_TCGA <- calculate_OR(vd_list_genes_TCGA)
    vd2 <- venn_diagram_within_databases(list = vd_list_genes_TCGA[1:2], title = paste0("TCGA \n OR: ",OR_TCGA))

    # Plot
    vd_plots <- list(vd1,vd2)
    output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/rare_germline_variants/replication/",genelist,"/WITHIN_databases_BETWEEN_NMD_methods_venn_diagram_random_genes_",random_genes,".png")
    png(output_path, width = 1750, height = 1000, res = 300)
    p <- cowplot::plot_grid( plotlist = vd_plots, labels = "AUTO", align = "v") #ncol = 4, nrow = 4)
    title_name <- paste0(genes," genes")
    title <- ggdraw() + draw_label(title_name, size = 18, fontface='bold')
    final_p <- plot_grid(title, p, ncol=1, rel_heights=c(0.1, 1)) # rel_heights values control title margins
    grid.newpage()
    grid.draw(final_p)
    dev.off()

}

# 2.3) Sig hits --> BETWEEN databases
# Split by tissue only and matched turmor-normal

# NMD_methods <- c("ASE","Endogenous")
# for (empirical_FDR in c(5,10,15,20,25)) {
#     for (database_discovery in databases) {
#         for (database_validation in databases) {
#             if (database_discovery == database_validation) {next}
#             for (NMD_method_discovery in NMD_methods) {
#                 for (NMD_method_validation in NMD_methods) {
#                     vd_plots <- matched_tissue_RGWAS_replication(NMD_method_discovery = NMD_method_discovery, NMD_method_validation = NMD_method_validation,
#                                                     empirical_FDR = empirical_FDR, database_discovery = database_discovery, database_validation = database_validation)
#                     output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/rare_germline_variants/replication/",geneList,"/between/tissues_matched_vdg_discovery_",database_discovery,"_discovery_",database_validation,"_validation_",NMD_method_discovery,"_discovery_",NMD_method_validation,"_validation_FDR_discovery_",empirical_FDR,"_genes_",randomGenes,".png")
#                     png(output_path, width = 5500, height = 4500, res = 300)
#                     p <- cowplot::plot_grid( plotlist = vd_plots, labels = "AUTO", align = "v") #ncol = 4, nrow = 4)
#                     title_name <- paste0("RGWAS\nDatabase Discovery: ",database_discovery," - Database Validation: ",database_validation,"\n",
#                                         "NMD method Discovery: ",NMD_method_discovery," - NMD method Validation: ",NMD_method_validation,"\n",
#                                         "Empirical FDR: ",empirical_FDR)
#                     title <- ggdraw() + draw_label(title_name, size = 18, fontface='bold')
#                     final_p <- plot_grid(title, p, ncol=1, rel_heights=c(0.1, 1)) # rel_heights values control title margins
#                     grid.newpage()
#                     grid.draw(final_p)
#                     dev.off()
#                 }
#             }
#         }
#     }
# }

NMD_methods <- c("ASE","Endogenous")
RGWAS_replication_hits_final <- c()
for (random_geneset_char in c("yes","no")) {
    for (empirical_FDR in c(1,2,3,4,5)) {
        for (database_discovery in databases) {
            for (database_validation in databases) {
                if (database_discovery == database_validation) {next}
                RGWAS_replication_hits <- matched_tissue_RGWAS_replication_no_strat( random_geneset_char = random_geneset_char, empirical_FDR = empirical_FDR,
                                                 database_discovery = database_discovery, database_validation = database_validation)
                RGWAS_replication_hits$empirical_FDR_threshold_used <- empirical_FDR 
                RGWAS_replication_hits_final <- rbind(RGWAS_replication_hits_final,RGWAS_replication_hits)
            }
        }
    }
}

# Proportion of replicated hits without stratification by dataset/NMD_method
# Number of replicated hits are normalized by the number of total tests in the discovery database

p <- RGWAS_replication_hits_final %>%
        group_by(empirical_FDR_threshold_used, database_discovery, database_validation, random_geneset, genes_tested) %>%
        summarise(num_replicated_hits = n()) %>%
        mutate(prop_replicated_hits = round(num_replicated_hits /  genes_tested,4)*100) %>% 
        mutate(database_validation = paste0(database_validation,"_val")) %>%
        mutate(database_discovery = paste0(database_discovery,"_disc")) %>% data.frame() %>%
                ggplot(aes(x = empirical_FDR_threshold_used, y = prop_replicated_hits, fill = random_geneset)) + 
                geom_bar(stat='identity',position=position_dodge()) +
                facet_grid(. ~ database_discovery + database_validation) +
                theme_bw(base_size = 25) +
                theme(axis.text.x = element_text(size = 18, angle = 45, hjust=1))
output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/rare_germline_variants/replication/",genelist,"/between/proportion_replicated_hits_without_stratification.png")
png(output_path, width = 4500, height = 4000, res = 300)
print(p)
dev.off()

# Split by NMD method

p <- RGWAS_replication_hits_final %>%
        group_by(empirical_FDR_threshold_used, database_discovery, database_validation, NMD_method,random_geneset, genes_tested) %>%
        summarise(num_replicated_hits = n()) %>%
        mutate(prop_replicated_hits = round(num_replicated_hits /  genes_tested,4)*100) %>% 
        mutate(database_validation = paste0(database_validation,"_val")) %>%
        mutate(database_discovery = paste0(database_discovery,"_disc")) %>% data.frame() %>%
                ggplot(aes(x = empirical_FDR_threshold_used, y = prop_replicated_hits, fill = random_geneset)) + 
                geom_bar(stat='identity',position=position_dodge()) +
                facet_grid(NMD_method ~ database_discovery + database_validation) +
                theme_bw(base_size = 25) +
                theme(axis.text.x = element_text(size = 18, angle = 45, hjust=1))
output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/rare_germline_variants/replication/",genelist,"/between/proportion_replicated_hits_without_stratification_NMD_method.png")
png(output_path, width = 4500, height = 4000, res = 300)
print(p)
dev.off()
        
##########################################################################################
##########################################################################################
##########################################################################################

# Proportion of replicated hits WITH stratification

RGWAS_replication_hits_final <- c()
NMD_methods <- c("ASE","Endogenous")
for (random_geneset_char in c("yes","no")) {
    for (empirical_FDR in c(1,2,3,4,5)) {
        for (database_discovery in databases) {
            for (database_validation in databases) {
                if (database_discovery == database_validation) {next}
                for (NMD_method_discovery in NMD_methods) {
                    for (NMD_method_validation in NMD_methods) {
                        RGWAS_replication_hits <- matched_tissue_RGWAS_replication(NMD_method_discovery = NMD_method_discovery, NMD_method_validation = NMD_method_validation, random_geneset_char = random_geneset_char,
                                                        empirical_FDR = empirical_FDR, database_discovery = database_discovery, database_validation = database_validation)
                        RGWAS_replication_hits$empirical_FDR_threshold_used <- empirical_FDR 
                        RGWAS_replication_hits_final <- rbind(RGWAS_replication_hits_final,RGWAS_replication_hits)
                    }
                }
            }
        }
    }
}

p <- RGWAS_replication_hits_final %>%
        group_by(empirical_FDR_threshold_used,database_discovery, database_validation, random_geneset, genes_tested) %>%
        summarise(num_replicated_hits = n()) %>%
        mutate(prop_replicated_hits = round(num_replicated_hits /  genes_tested,4)*100) %>% 
        mutate(database_validation = paste0(database_validation,"_val")) %>%
        mutate(database_discovery = paste0(database_discovery,"_disc")) %>% data.frame() %>%
            ggplot(aes(x = empirical_FDR_threshold_used, y = prop_replicated_hits, fill = random_geneset)) + 
                geom_bar(stat='identity',position=position_dodge()) +
                facet_grid(. ~ database_discovery + database_validation) +
                theme_bw(base_size = 25) +
                theme(axis.text.x = element_text(size = 18, angle = 45, hjust=1))
output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/rare_germline_variants/replication/",genelist,"/between/proportion_replicated_hits_stratification.png")
png(output_path, width = 4500, height = 4000, res = 300)
print(p)
dev.off()

# Proportion of replicated hits without stratification

NMD_methods <- c("ASE","Endogenous")
RGWAS_replication_hits_final <- c()
for (random_geneset_char in c("yes","no")) {
    for (empirical_FDR in c(1,2,3,4,5)) {
        for (database_discovery in databases) {
            for (database_validation in databases) {
                if (database_discovery == database_validation) {next}
                RGWAS_replication_hits <- matched_tissue_RGWAS_replication_no_strat( random_geneset_char = random_geneset_char,
                                                empirical_FDR = empirical_FDR, database_discovery = database_discovery, database_validation = database_validation)
                RGWAS_replication_hits$empirical_FDR_threshold_used <- empirical_FDR 
                RGWAS_replication_hits_final <- rbind(RGWAS_replication_hits_final,RGWAS_replication_hits)
            }
        }
    }
}

p <- RGWAS_replication_hits_final %>%
        group_by(empirical_FDR_threshold_used, database_discovery, database_validation, random_geneset, genes_tested) %>%
        summarise(num_replicated_hits = n()) %>%
        mutate(prop_replicated_hits = round(num_replicated_hits /  genes_tested,4)*100) %>% 
        mutate(database_validation = paste0(database_validation,"_val")) %>%
        mutate(database_discovery = paste0(database_discovery,"_disc")) %>% data.frame() %>%
                ggplot(aes(x = empirical_FDR_threshold_used, y = prop_replicated_hits, fill = random_geneset)) + 
                geom_bar(stat='identity',position=position_dodge()) +
                facet_grid(. ~ database_discovery + database_validation) +
                theme_bw(base_size = 25) +
                theme(axis.text.x = element_text(size = 18, angle = 45, hjust=1))
output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/rare_germline_variants/replication/",genelist,"/between/proportion_replicated_hits_without_stratification.png")
png(output_path, width = 4500, height = 4000, res = 300)
print(p)
dev.off()
        

# Proportion of replicated hits WITH stratification and using NEW EMPIRICAL FDR CORRECTION

RGWAS_replication_hits_final <- c()
NMD_methods <- c("ASE","Endogenous")
for (random_geneset_char in c("yes","no")) {
    for (empirical_FDR in c(1,2,3,4,5)) {
        for (database_discovery in databases) {
            for (database_validation in databases) {
                if (database_discovery == database_validation) {next}
                for (NMD_method_discovery in NMD_methods) {
                    for (NMD_method_validation in NMD_methods) {
                        RGWAS_replication_hits <- matched_tissue_RGWAS_replication_newEmpiricalFDR(NMD_method_discovery = NMD_method_discovery, NMD_method_validation = NMD_method_validation, random_geneset_char = random_geneset_char,
                                                        empirical_FDR = empirical_FDR, database_discovery = database_discovery, database_validation = database_validation)
                        RGWAS_replication_hits$empirical_FDR_threshold_used <- empirical_FDR 
                        RGWAS_replication_hits_final <- rbind(RGWAS_replication_hits_final,RGWAS_replication_hits)
                    }
                }
            }
        }
    }
}

p <- RGWAS_replication_hits_final %>%
        group_by(empirical_FDR_threshold_used,database_discovery, database_validation, random_geneset, genes_tested) %>%
        summarise(num_replicated_hits = n()) %>%
        mutate(prop_replicated_hits = round(num_replicated_hits /  genes_tested,4)*100) %>% 
        mutate(database_validation = paste0(database_validation,"_val")) %>%
        mutate(database_discovery = paste0(database_discovery,"_disc")) %>% data.frame() %>%
            ggplot(aes(x = empirical_FDR_threshold_used, y = prop_replicated_hits, fill = random_geneset)) + 
                geom_bar(stat='identity',position=position_dodge()) +
                facet_grid(. ~ database_discovery + database_validation) +
                theme_bw(base_size = 25) +
                theme(axis.text.x = element_text(size = 18, angle = 45, hjust=1))
output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/rare_germline_variants/replication/",genelist,"/between/proportion_replicated_hits_stratification_NEW_FDR.png")
png(output_path, width = 4500, height = 4000, res = 300)
print(p)
dev.off()

# Proportion of replicated hits WITHOUT stratification and using NEW EMPIRICAL FDR CORRECTION

RGWAS_replication_hits_final <- c()
NMD_methods <- c("ASE","Endogenous")
for (random_geneset_char in c("yes","no")) {
    for (empirical_FDR in c(1,2,3,4,5,10)) {
        for (database_discovery in databases) {
            for (database_validation in databases) {
                if (database_discovery == database_validation) {next}
                        RGWAS_replication_hits <- matched_tissue_RGWAS_replication_newEmpiricalFDR_no_strat(random_geneset_char = random_geneset_char,
                                                        empirical_FDR = empirical_FDR, database_discovery = database_discovery, database_validation = database_validation)
                        RGWAS_replication_hits$empirical_FDR_threshold_used <- empirical_FDR 
                        RGWAS_replication_hits_final <- rbind(RGWAS_replication_hits_final,RGWAS_replication_hits)
            }
        }
    }
}

p <- RGWAS_replication_hits_final %>% 
        filter(matched_tissue!="pantissue") %>%
        group_by(empirical_FDR_threshold_used,database_discovery, database_validation, random_geneset, genes_tested) %>%
        summarise(num_replicated_hits = n()) %>%
        mutate(prop_replicated_hits = round(num_replicated_hits /  genes_tested,4)*100) %>% 
        mutate(database_validation = paste0(database_validation,"_val")) %>%
        mutate(database_discovery = paste0(database_discovery,"_disc")) %>% data.frame() %>%
            ggplot(aes(x = empirical_FDR_threshold_used, y = prop_replicated_hits, fill = random_geneset)) + 
                geom_bar(stat='identity',position=position_dodge()) +
                facet_grid(. ~ database_discovery + database_validation) +
                theme_bw(base_size = 25) +
                theme(axis.text.x = element_text(size = 18, angle = 45, hjust=1))
output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/rare_germline_variants/replication/",genelist,"/between/proportion_replicated_hits_without_stratification_NEW_FDR.png")
png(output_path, width = 4500, height = 4000, res = 300)
print(p)
dev.off()


# Proportion of replicated hits WITHOUT stratification and using NORMAL FDR CORRECTION

RGWAS_replication_hits_final <- c()
for (random_geneset in c("yes","no")) {
    for (FDR in c(1,2,3,4,5,10,15,20,25,30,35,40,45,50)) {
        for (database_discovery in databases) {
            for (database_validation in databases) {
                if (database_discovery == database_validation) {next}
                    RGWAS_replication_hits <- matched_tissue_RGWAS_replication_normal_FDR_no_strat(random_geneset_char = random_geneset,
                                                    FDR = FDR, database_discovery = database_discovery, database_validation = database_validation)
                    if (!is.null(RGWAS_replication_hits)) {
                        RGWAS_replication_hits$FDR_threshold_used <- FDR 
                        RGWAS_replication_hits_final <- rbind(RGWAS_replication_hits_final,RGWAS_replication_hits)
                    }
            }
        }
    }
}

p <- RGWAS_replication_hits_final %>%
        filter(matched_tissue != "pantissue" & randomization == "no") %>%
        group_by(FDR_threshold_used,random_geneset,database_discovery, database_validation) %>%
        summarise(num_replicated_hits = n()) %>%
        mutate(database_validation = paste0(database_validation,"_val")) %>%
        mutate(database_discovery = paste0(database_discovery,"_disc")) %>% data.frame() %>%
            ggplot(aes(x = factor(FDR_threshold_used), y = num_replicated_hits, fill = random_geneset)) + 
                geom_bar(stat='identity',position=position_dodge()) +
                facet_grid(. ~ database_discovery + database_validation) +
                theme_bw(base_size = 25) + #scale_x_discrete(labels=c("1","2","3","4","5","10","25")) +
                theme(axis.text.x = element_text(size = 18, angle = 45, hjust=1))
output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/rare_germline_variants/replication/",genelist,"/between/num_replicated_hits_without_stratification_normal_FDR.png")
png(output_path, width = 4500, height = 4000, res = 300)
print(p)
dev.off()

##########################################################################################
##########################################################################################
##########################################################################################

# 2.3.1) Barplots of numbers of replicated hits
# By all
RGWAS_replication_hits_final <- RGWAS_replication_hits_final %>%
    # mutate(NMD_method_validation = paste0(NMD_method_validation,"_val")) %>%
    # mutate(NMD_method_discovery = paste0(NMD_method_discovery,"_disc")) %>%
    mutate(database_validation = paste0(database_validation,"_val")) %>%
    mutate(database_discovery = paste0(database_discovery,"_disc"))

p <- RGWAS_replication_hits_final %>%
        group_by(matched_tissue, dataset,NMD_method,database_discovery,database_validation) %>%
        filter(empirical_FDR_threshold_used == 5 & random_geneset == "no") %>%
        filter(!dataset %in% c("CCR_90orhigher_0.1perc","MRT_25thCentile_0.1perc")) %>%
        summarise(replicated_hits = n()) %>%
        arrange(replicated_hits) %>%
        #mutate(matched_tissue = factor(matched_tissue, levels = rev(matched_tissue))) %>%
        ggplot(aes(x = dataset, y = replicated_hits, fill = matched_tissue)) + 
            geom_bar(stat='identity',position=position_dodge()) +
            #ylim(c(0,0.25)) +
            #facet_grid(. ~ NMD_method_discovery + NMD_method_validation) +
            #facet_grid(. ~ database_discovery + database_validation) +
            facet_grid(NMD_method ~ database_discovery + database_validation) +
            theme_bw(base_size = 25) +
            theme(axis.text.x = element_text(size = 18, angle = 45, hjust=1))
output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/rare_germline_variants/replication/",genelist,"/between/Number_replicated_hits_by_all.png")
png(output_path, width = 4000, height = 4500, res = 300)
print(p)
dev.off()

# By Dataset
p <- RGWAS_replication_hits_final %>%
        filter(empirical_FDR_threshold_used == 5 & random_geneset == "no") %>%
        filter(!dataset %in% c("CCR_90orhigher_0.1perc","MRT_25thCentile_0.1perc")) %>%
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
        filter(empirical_FDR_threshold_used == 5 & random_geneset == "no") %>%
        filter(!dataset %in% c("CCR_90orhigher_0.1perc","MRT_25thCentile_0.1perc")) %>%
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
        filter(empirical_FDR_threshold_used == 5 & random_geneset == "no") %>%
        filter(!dataset %in% c("CCR_90orhigher_0.1perc","MRT_25thCentile_0.1perc")) %>%
        group_by(matched_tissue) %>%
        summarise(replicated_hits = n()) %>%
        arrange(replicated_hits) %>%
        mutate(matched_tissue = factor(matched_tissue, levels = rev(matched_tissue))) %>%
        ggplot(aes(x = matched_tissue, y = replicated_hits)) + 
            geom_bar(stat='identity',position=position_dodge()) +
            theme_bw(base_size = 25) +
            #facet_grid(. ~ NMD_method_validation) +
            theme(axis.text.x = element_text(size = 18, angle = 45, hjust=1))
output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/rare_germline_variants/replication/",genelist,"/between/Number_replicated_hits_by_tissue.png")
png(output_path, width = 4000, height = 4500, res = 300)
print(p)
dev.off()
# By gene
# By tissue
p <- RGWAS_replication_hits_final %>%
        filter(empirical_FDR_threshold_used == 5 & random_geneset == "no") %>%
        filter(!dataset %in% c("CCR_90orhigher_0.1perc","MRT_25thCentile_0.1perc")) %>%
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

# 2.3.2) Genes table
RGWAS_replication_hits_final$gene_tissue <- paste0(RGWAS_replication_hits_final$variant," - ",RGWAS_replication_hits_final$matched_tissue)
#RGWAS_replication_hits_final$NMD_method_translation <- gsub("Endogenous","END",paste0(RGWAS_replication_hits_final$NMD_method_discovery,"/",RGWAS_replication_hits_final$NMD_method_validation))
RGWAS_replication_hits_final$database_translation <- paste0(RGWAS_replication_hits_final$database_discovery,"/",RGWAS_replication_hits_final$database_validation)
# RGWAS_replication_hits_table <- RGWAS_replication_hits_final %>%
#             filter(empirical_FDR_threshold_used == 5 & random_geneset == "no") %>%
#             filter(!dataset %in% c("CCR_90orhigher_0.1perc","MRT_25thCentile_0.1perc")) %>%
#             group_by(gene_tissue,dataset, database_translation) %>%
#             summarise(count = n()) %>%
#             ungroup()
RGWAS_replication_hits_table <- RGWAS_replication_hits_final %>%
            filter(empirical_FDR_threshold_used == 2 & random_geneset == "yes") %>%
            filter(!dataset %in% c("CCR_90orhigher_0.1perc","MRT_25thCentile_0.1perc")) %>%
            group_by(gene_tissue,dataset, database_translation)
#RGWAS_replication_hits_table <- RGWAS_replication_hits_table[order(RGWAS_replication_hits_table$gene_tissue, decreasing = FALSE),]
p <- RGWAS_replication_hits_table %>%
            ggplot(aes(x = dataset,y = gene_tissue)) +
            geom_tile(aes(fill = n_carriers), size = 1) +
            geom_text(aes(label = n_carriers),color = "black", size=6) +
            theme_bw(base_size = 24) +
            facet_grid(. ~ database_translation + NMD_method) +
            theme(axis.text.x = element_text(size = 18, angle = 45, hjust=1),
                    legend.title=element_text(size=16)) +
            scale_fill_gradientn(colours = c('grey','#F21A00')) +
            labs(fill = '# of individuals')
output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/rare_germline_variants/replication/",genelist,"/between/table_genes_replicated.png")
png(output_path, width = 5000, height = 5500, res = 300)
print(p)
dev.off()

# 2.4) Confident sig and replicated hits
# 2.4.1) Overlap between 2.2 (WITHIN) and 2.3 (BETWEEN)
# Thus, overlap between WITHIN and BETWEEN
# sig_hits_total
# RGWAS_replication_hits_final

# Hacer un solo venn diagram de replicated hits across the combinations (confident sets?)

# 2.5) Frequency of significant sig hits across tissues + its coefficients. Also in End and GTEx if possible









