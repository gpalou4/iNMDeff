GWAS_df <- function(NMD_method, VAF, TCGA_cancers, PCs) {

    if (NMD_method == "ASE") {
        NMD_method_VAF <- paste0(NMD_method,"_",VAF)
    } else {
        NMD_method_VAF <- NMD_method
    }

    dir_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/somatic_associations/")
    plots_path <- paste0(dir_path,"/analysis_plots")
    GWAS_res_path <-paste0(dir_path,"/final_associations")

    if (! dir.exists(plots_path) ) {
        dir.create(plots_path)
    }

    GWAS_df_cancer_all <- c()
    for (TCGA_cancer in TCGA_cancers) {
        print(TCGA_cancer)
        # GWAS results
        GWAS_res_name <- paste0(TCGA_cancer,"_",NMD_method_VAF,"_CGC_somatic_mut_CNV_PCs_",PCs,".txt")
        tryCatch({
            error <<- FALSE
            GWAS_res_cancer <- read.table( file = paste0(GWAS_res_path,"/",GWAS_res_name),
                            sep = "\t", header = TRUE, stringsAsFactors = FALSE)
            GWAS_res_cancer$cancer_type <- TCGA_cancer
        },error = function(e) {
            print("Cancer GWAS results not found...")
            error <<- TRUE
        })
        if (isTRUE(error)) {next}
        if (nrow(GWAS_res_cancer) == 0 ) {next}
        if (length(GWAS_df_cancer_all) == 0) {
            GWAS_df_cancer_all <- GWAS_res_cancer
        } else {
            GWAS_df_cancer_all <- rbind(GWAS_df_cancer_all,GWAS_res_cancer)
        }
    }
    return(GWAS_df_cancer_all)
}

calculate_lambda <- function(p_values) {
    # Lambda calculation
    chisq <- qchisq(1-p_values,1)
    lambda = round(median(chisq,na.rm=TRUE)/qchisq(0.5,1),2)
    return(lambda)
}

QQplot_and_lambdas <- function(data, qqplot_df = NULL) {
            
    lambdas <- data.frame(NMD_method = NA, dataset = NA, genesets = NA, cancer_type = NA)
    # Obtain name info
    cancer_type <- as.character(unique(data$cancer_type))
    genesets_name_char <- as.character(unique(data$genesets))
    dataset_name <- as.character(unique(data$dataset))
    NMD_method_name <- as.character(unique(data$NMD_method))
    # Paths
    dir_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/somatic_associations")
    plots_path <- paste0(dir_path,"/analysis_plots")
    # Sig genes of interest
    DLBC_genes <- c("SPECC1","ACVR1")
    PCPG_genes <- c("PSMA1","NUP98","MYOD1","LMO1","HRAS","FANCF","CARS","USP47","DKK3")
    CESC_genes <- c("PTPN11")

    if (cancer_type == "TCGA-DLBC") {
        genes <- DLBC_genes
    } else if (cancer_type == "TCGA-PCPG") {
        genes <- PCPG_genes
    } else if (cancer_type == "TCGA-CESC") {
        genes <- CESC_genes
    } else if (cancer_type == "TCGA-LUAD") {
        genes <- "TLX1"
    } else if (cancer_type == "TCGA-MESO") {
        genes <- "SMO"
    } else if (cancer_type == "TCGA-BRCA") {
        genes <- "CDH1"
    } else {
        genes <- "SMG5"
    }

    # 1) Lambdas
    # Calculate lambda
    # Remove p-values with NA
    data_filt <- data[!is.na(data[,"p_value"]),]
    tryCatch({
        error <<- FALSE
        p_values_random <- data_filt[data_filt$genesets == "Random","p_value"] %>% pull()
        p_values_non_random <- data_filt[data_filt$genesets == "Selected","p_value"]  %>% pull()
    },error = function(e) {
        error <<- TRUE
    })
    if (isTRUE(error)) {
        p_values_random <- data_filt[data_filt$genesets == "Random","p_value"]
        p_values_non_random <- data_filt[data_filt$genesets == "Selected","p_value"]
    }
    lambda_random <- calculate_lambda(p_values_random)
    lambda_non_random <- calculate_lambda(p_values_non_random)
    # 2) Save lambda
    lambdas$dataset <- dataset_name
    lambdas$NMD_method <- NMD_method_name
    lambdas$lambda_non_random <- as.numeric(lambda_non_random)
    lambdas$lambda_random <- as.numeric(lambda_random)
    lambdas$cancer_type <- cancer_type
    #lambdas$genesets <- genesets_name_char
    
    # 3) QQplot
    if (cancer_type %in% c("TCGA-PCPG","TCGA-DLBC","TCGA-CESC","TCGA-LUAD","TCGA-MESO","TCGA-BRCA","pancancer")) {
        # 3.1) Non Random Genes
        plot_title <- paste0("CGC + NMD Genes --> Lambda = ",lambda_non_random)
        df_non_random <- data_filt[data_filt$genesets == "Selected",]
        df_non_random <- df_non_random[order(df_non_random[,"p_value"]),]
        ylim_max <- round(max(c(-log10(p_values_random),-log10(p_values_non_random)), na.rm = TRUE),2)+0.25

        # Expected p-values
        n1 <- length(p_values_non_random)
        if (n1 != 0) {
            a <- 1:n1
            pvals_exp <- a/n1
            df_non_random[,"p_value_expected"] <- pvals_exp    
            df_non_random$lambda <- lambda_non_random               
            # QQplot       
            p1 <- ggplot(df_non_random, aes(x = -log10(eval(parse(text="p_value_expected"))), y = -log10(eval(parse(text="p_value"))) )) +
                geom_point(size = 2) + ylim(c(0,ylim_max)) +
                geom_abline(intercept = 0, slope = 1) +
                geom_label_repel(data = subset(df_non_random, Gene_symbol %in% c(genes)), max.overlaps = 100, 
                    aes(label = Gene_symbol), color = "black", size = 5, arrow = arrow(length = unit(0.01, 'npc')),
                        point.padding = 0.2, nudge_x = .08, nudge_y = .08) +
                ylab("-log10(P-values)") + xlab("-log10(expected P-values)") + ggtitle(plot_title) +
                theme_classic(base_size = 16)
        }
        # 3.2) Random Genes
        plot_title <- paste0("Random Genes --> Lambda = ",lambda_random)
        df_random <- data_filt[data_filt$genesets == "Random",]
        df_random <- df_random[order(df_random[,"p_value"]),]

        # Expected p-values
        n2 <- length(p_values_random)
        if (n2 != 0) {
            a <- 1:n2
            pvals_exp <- a/n2
            df_random[,"p_value_expected"] <- pvals_exp
            df_random$lambda <- lambda_random
            # QQplot       
            p2 <- ggplot(df_random, aes(x = -log10(eval(parse(text="p_value_expected"))), y = -log10(eval(parse(text="p_value"))) )) +
                geom_point(size = 2) + ylim(c(0,ylim_max)) +
                geom_abline(intercept = 0, slope = 1) +
                ylab("-log10(P-values)") + xlab("-log10(expected P-values)") + ggtitle(plot_title) +
                theme_classic(base_size = 16)
        }
        if ( (n1 != 0) & (n2 != 0) ) {
            # Final qqplot
            qqplot <- paste0(plots_path,"/",cancer_type,"_",NMD_method_name,"_qqplot_genes_",dataset_name,".png")
            p <- cowplot::plot_grid(plotlist=list(p1,p2), labels = "AUTO", align = "v", ncol = 2, nrow = 1)
            ggsave(filename = qqplot, device = "png", plot = p, width = 5000, height = 2500, dpi = 300, units = "px")
        }
        if (!is.null(qqplot_df)) {
            return(list(df_random = df_random, df_non_random = df_non_random))
        }
    }
    return(lambdas)
}

GWAS_replication_strat_cancer <- function(GWAS_df_all, genesets_char, FDR, NMD_method_discovery, NMD_method_validation) {
    FDR <- FDR / 100
    # genesets_char <- genesets
    GWAS_replicated_final_df <- c()
    # 1) Sig hits in discovery
    GWAS_discovery_sig <- GWAS_df_all %>%
                    filter(genesets == genesets_char) %>%
                    group_by(cancer_type, dataset) %>%
                    mutate(pvalue_FDR_adjusted = p.adjust(p_value, method = "fdr")) %>%
                    mutate(sig_genes_ID = ifelse( ((pvalue_FDR_adjusted < FDR) ), paste0(NMD_method,"_",Gene_symbol,"_",dataset,"_",cancer_type), NA)) %>%
                    mutate(discovery_hits_ID = ifelse( ((pvalue_FDR_adjusted < FDR) ), paste0(NMD_method,"_",Gene_symbol,"_",dataset,"_",cancer_type), NA))
    sig_genes_ID_char <- as.character(na.omit(GWAS_discovery_sig$sig_genes_ID))

    # 2) Validation of hits in NMD_method discovery within the same cancer type
    for (i in 1:length(TCGA_cancers)) {
        cancer <- TCGA_cancers[i]   
        #print(cancer)
        # Sig hits IDs in the cancer
        sig_genes_ID_char_to_validate <- sig_genes_ID_char[grep(cancer,sig_genes_ID_char)]
        # Sig hits in the NMD_method discovery
        index <- grep(NMD_method_discovery,sig_genes_ID_char_to_validate)
        if (length(index) != 0) {
            sig_genes_ID_char_cancer <- sig_genes_ID_char_to_validate[grep(NMD_method_discovery,sig_genes_ID_char_to_validate)]
        } else {
            next
        }
        # Swap by NMD_method validation
        sig_genes_ID_char_cancer <- gsub(NMD_method_discovery,NMD_method_validation,sig_genes_ID_char_cancer)
        sig_genes_ID_char_cancer <- sig_genes_ID_char_cancer[!duplicated(sig_genes_ID_char_cancer)]

        # FDR in NMD_method validation (only for the subset of genes that are significant in discovery)
        GWAS_FDR_validation <- GWAS_discovery_sig %>%
            filter(NMD_method == NMD_method_validation & genesets == genesets_char) %>%
            filter(paste0(NMD_method,"_",Gene_symbol,"_",dataset,"_",cancer_type) %in% sig_genes_ID_char_cancer) %>%
            group_by(cancer_type,dataset) %>%
            mutate(pvalue_FDR_adjusted = p.adjust(p_value, method = "fdr")) %>%
            filter(pvalue_FDR_adjusted < FDR) %>%
            mutate(replicated_hits_ID = paste0(NMD_method,"_",Gene_symbol,"_",dataset,"_",cancer_type))
        GWAS_validation_sig_replicated <- GWAS_FDR_validation
        if (nrow(GWAS_validation_sig_replicated) == 0 ) {next}
        print(GWAS_validation_sig_replicated)
        
        # # 3) Obtain the info of the replicated hits in discovery and validation
        # Dataset & NMD_method
        # Discovery hits ID

        # GWAS_discovery_sig_hits_ID <- GWAS_discovery_sig %>%
        #                 filter(NMD_method == NMD_method_discovery & genesets == genesets_char) %>%
        #                 filter(!is.na(discovery_hits_ID)) %>%
        #                 mutate(tmp_ID = paste0(NMD_method,"_",Gene_symbol,"_",dataset,"_",cancer_type))
        # # Check which validated genes-cancer were in discovery
        # sig_genes_ID_char_to_validate <- sig_genes_ID_char[grep(cancer,sig_genes_ID_char)]
        # # Sig hits in the NMD_method discovery
        # index <- grep(NMD_method_discovery,sig_genes_ID_char_to_validate)
        # if (length(index) != 0) {
        #     sig_genes_ID_char_cancer <- sig_genes_ID_char_to_validate[grep(NMD_method_discovery,sig_genes_ID_char_to_validate)]
        # }

        GWAS_discovery_sig_hits_ID <- GWAS_discovery_sig %>%
                        filter(NMD_method == NMD_method_discovery & genesets == genesets_char) %>%
                        filter(!is.na(discovery_hits_ID)) %>%
                        mutate(replicated_hits_ID = paste0(NMD_method,"_",Gene_symbol,"_",dataset,"_",cancer_type))
        # Discovery Info
        tmp_ID <- gsub(NMD_method_validation,NMD_method_discovery,GWAS_validation_sig_replicated$replicated_hits_ID)
        GWAS_discovery_sig_replicated <- GWAS_discovery_sig_hits_ID[GWAS_discovery_sig_hits_ID$replicated_hits_ID %in% tmp_ID,]

        # 4) Merge Replicated hits: from Discovery & Validated
        GWAS_final_replicated_hits <- rbind(GWAS_discovery_sig_replicated,GWAS_validation_sig_replicated)
        GWAS_final_replicated_hits$matched_cancer <- cancer
        GWAS_final_replicated_hits$NMD_method_discovery <- NMD_method_discovery
        GWAS_final_replicated_hits$NMD_method_validation <- NMD_method_validation
        GWAS_final_replicated_hits$disc_val <- NA
        GWAS_final_replicated_hits <- GWAS_final_replicated_hits %>%
                                        mutate(disc_val = ifelse( ((NMD_method %in% NMD_method_discovery) ), "discovery", "validation"))
        # Number of genes tested at discovery
        num_genes_tested <- GWAS_df_all %>%
                            filter(genesets == genesets_char & NMD_method == NMD_method_discovery) %>%
                            summarise(genes_tested = n()) %>% data.frame()
        GWAS_final_replicated_hits$genes_tested <- num_genes_tested$genes_tested
        GWAS_replicated_final_df <- rbind(GWAS_replicated_final_df,GWAS_final_replicated_hits)
    }
    return(GWAS_replicated_final_df)
}

library(dplyr)
library(ggplot2)
# facet_nested
library("ggh4x")
library(RColorBrewer)
library(cowplot)
library(ggrepel)
library(scales)

# 1) Data
# 1.1) TCGA cancer type
TCGA_cancers <- read.table(file = "/g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/TCGA_projects_names.txt", 
                    header = FALSE, stringsAsFactors = FALSE, sep ="\t")$V1
TCGA_cancers <- c(TCGA_cancers,"pancancer")
TCGA_cancers <- gsub("TCGA-","", TCGA_cancers)
# 1.2) CGC gene info
CGC <- read.csv(file = "/g/strcombio/fsupek_cancer1/gpalou/COSMIC/cancer_gene_census_updated.tsv", 
                    header = TRUE, stringsAsFactors = FALSE)
# 1.3) NMD factors
NMD_genes <- read.table(file = "/g/strcombio/fsupek_home/gpalou/data/NMD_project/conversor_tables/ENSEMBL/NMD_genes.txt",
                header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# 1.4) ENSEMBL transcripts IDs hg38 GTF
ensembl_v88_gtf <- rtracklayer::import("/g/strcombio/fsupek_home/gpalou/data/NMD_project/conversor_tables/ENSEMBL/gencode.v26.annotation.gtf")
ensembl_v88_gtf <- as.data.frame(ensembl_v88_gtf)
ensembl_v88_gtf_filt <- ensembl_v88_gtf[ensembl_v88_gtf$type == "gene",]
ensembl_v88_gtf_filt$genome_location <- gsub("chr","",paste0(ensembl_v88_gtf_filt$seqnames,":",ensembl_v88_gtf_filt$start,"-",ensembl_v88_gtf_filt$end))
ensembl_v88_gtf_filt <- ensembl_v88_gtf_filt[,c("gene_name","genome_location")]

# 2) Obtain GWAS results
GWAS_end <- GWAS_df(NMD_method = "endogenous", VAF = 0.2, TCGA_cancers = TCGA_cancers, PCs = 1)
GWAS_end$NMD_method <- "Endogenous"
GWAS_ASE <- GWAS_df(NMD_method = "ASE", VAF = 0.2, TCGA_cancers = TCGA_cancers, PCs = 1)
GWAS_ASE$NMD_method <- "ASE"
GWAS_all <- rbind(GWAS_ASE,GWAS_end)

# 3) Add info
# CGC gene
GWAS_all <- merge(GWAS_all,CGC, by.x = c("Gene_symbol","ENSEMBL_gene"), by.y = c("Gene.Symbol","Synonyms"), all.x = TRUE)
# ENSEMBL chr & genome location
GWAS_all <- merge(GWAS_all, ensembl_v88_gtf_filt, by.x = "Gene_symbol", by.y ="gene_name", all.x = TRUE)
GWAS_all$chr <- gsub("(.*)\\:.*","\\1",GWAS_all$genome_location)
# Fix TSG/OG label
GWAS_all$Role.in.Cancer <- ifelse(GWAS_all$Role.in.Cancer == "oncogene, fusion","oncogene",GWAS_all$Role.in.Cancer)
GWAS_all$Role.in.Cancer <- ifelse(GWAS_all$Role.in.Cancer == "TSG, fusion","TSG",GWAS_all$Role.in.Cancer)
GWAS_all$Role.in.Cancer <- ifelse(GWAS_all$Role.in.Cancer == "oncogene, TSG, fusion","mix",GWAS_all$Role.in.Cancer)
GWAS_all$Role.in.Cancer <- ifelse(GWAS_all$Role.in.Cancer == "oncogene, TSG","mix",GWAS_all$Role.in.Cancer)
GWAS_all$Role.in.Cancer <- ifelse(GWAS_all$Role.in.Cancer == "",NA,GWAS_all$Role.in.Cancer)
GWAS_all$Role.in.Cancer <- ifelse(GWAS_all$Role.in.Cancer == "fusion","fusion",GWAS_all$Role.in.Cancer)
table(GWAS_all$Role.in.Cancer)
# Obtain gene classification (NMD, cancer, random)
CGC_genes <- CGC[which(CGC$Role.in.Cancer != "non_CGC_NMD"),"Gene.Symbol"]
NMD_genes <- CGC[which(!is.na(CGC$NMD_type)),"Gene.Symbol"]
random_genes <- CGC[which(CGC$Role.in.Cancer == "non_CGC_NMD"),"Gene.Symbol"]
GWAS_all$genesets <- ""
GWAS_all[GWAS_all$Gene_symbol %in% c(CGC_genes,NMD_genes),"genesets"] <- "Selected"
GWAS_all[!GWAS_all$Gene_symbol %in% c(CGC_genes,NMD_genes),"genesets"] <- "Random"
# Split by dataset
df_tmp <- GWAS_all[,grep("coeff",colnames(GWAS_all))]
df_tmp_stacked <- stack(df_tmp)
colnames(df_tmp_stacked) <- c("beta_coefficient","dataset")
df_tmp_2 <- GWAS_all[,grep("pval",colnames(GWAS_all))]
df_tmp_stacked_2 <- stack(df_tmp_2)
colnames(df_tmp_stacked_2) <- c("p_value","dataset")
# Merge
df_tmp_stacked$dataset <- gsub("_coeff","",df_tmp_stacked$dataset)
df_final_tmp <- cbind(df_tmp_stacked,df_tmp_stacked_2[,1,drop=FALSE])
# Add info
df_tmp_3 <- GWAS_all[,-grep("som",colnames(GWAS_all))]
df_final <- c()
for (i in 1:5) {
    print(i)
    df_final <- rbind(df_final,df_tmp_3)
}
GWAS_final_df <- cbind(df_final_tmp,df_final)

output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/somatic_associations/final_associations/GWAS_somatic_mut_CNV.txt")
GWAS_final_df <- read.table(file = output_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# write.table(data.frame(GWAS_final_df), file = output_path, sep = "\t", quote = FALSE,
#                        col.names = TRUE, row.names = FALSE)

# 4) Plots
# 4.1) pot of Lambdas across cancers
lambdas_res <- GWAS_final_df %>%
                    group_by(NMD_method,cancer_type,dataset) %>%
                    do( QQplot_and_lambdas(.) )
df <- GWAS_final_df %>% filter(NMD_method == "Endogenous" & cancer_type == "pancancer" & dataset == "som_CNV_amp") 
pancancer_df_list <- QQplot_and_lambdas(data = df, qqplot_df = "yes")
# Save
# write.table(pancancer_df_list, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig8/SuppFig8A.txt", 
#                 sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
saveRDS(pancancer_df_list, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig8/SuppFig8A.RData")
df <- GWAS_final_df %>% filter(NMD_method == "ASE" & cancer_type == "pancancer" & dataset == "som_CNV_amp") 
pancancer_df_list <- QQplot_and_lambdas(data = df, qqplot_df = "yes")
# Save
# write.table(pancancer_df_list, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig8/SuppFig8B.txt", 
#                 sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
saveRDS(pancancer_df_list, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig8/SuppFig8B.RData")

lambdas_res <- lambdas_res %>%
  mutate(dataset = case_when(
    dataset == "som_mut_truncating" ~ "Truncating",
    dataset == "som_mut_missense" ~ "Missense",
    dataset == "som_mut_synonymous" ~ "Synonymous",
    dataset == "som_CNV_amp" ~ "CNA amp",
    dataset == "som_CNV_del" ~ "CNA del",
    TRUE ~ dataset
  ))
lambdas_res$cancer_type <- gsub("TCGA-","",lambdas_res$cancer_type)

# Non Random genes
p <- ggplot(lambdas_res, aes(x = cancer_type, y = lambda_non_random, color = factor(dataset))) + 
        geom_point(size = 5) +
        ylim(c(0,4)) + xlab("") + ylab("Lambda") + ggtitle("Cancer and NMD genes") +
        facet_grid(. ~ NMD_method) +
        theme_bw(base_size = 30) +
        theme(axis.text.x = element_text(size = 25, angle = 45, hjust=1),
            plot.title = element_text(size = 35, hjust = 0.5),
            legend.position="top",
            strip.text = element_text(size = 32),
            legend.text = element_text(size=25, face="bold"),
            legend.title = element_blank()) +
        guides(colour = guide_legend(override.aes = list(size=10))) 

output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/somatic_associations/analysis_plots/QC/by_cancers/TCGA_lambdas_non_random_genes.png")
png(output_path, width = 8500, height = 3500, res = 300)
print(p)
dev.off()

# Random genes
p <- ggplot(lambdas_res, aes(x = cancer_type, y = lambda_random, color = factor(dataset))) + 
        geom_point(size = 5) +
        ylim(c(0,4)) + xlab("") + ylab("Lambda") + ggtitle("Random genes") +
        facet_grid(. ~ NMD_method) +
        theme_bw(base_size = 30) +
        theme(axis.text.x = element_text(size = 25, angle = 45, hjust=1),
            plot.title = element_text(size = 35, hjust = 0.5),
            legend.position="top",
            strip.text = element_text(size = 32),
            legend.text = element_text(size=25, face="bold"),
            legend.title = element_blank()) +
        guides(colour = guide_legend(override.aes = list(size=10))) 

output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/somatic_associations/analysis_plots/QC/by_cancers/TCGA_lambdas_random_genes.png")
png(output_path, width = 8500, height = 3500, res = 300)
print(p)
dev.off()

lambdas_res$NMD_method <- ifelse(lambdas_res$NMD_method == "Endogenous","ETG","ASE")
# Save
write.table(lambdas_res, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig7/SuppFig7A_B.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
saveRDS(lambdas_res, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig7/SuppFig7A_B.RData")

# 4.2) Remove models with high lambda values
TCGA_models_to_remove_non_rand <- lambdas_res[which(lambdas_res$lambda_non_random >= 1.5),]
TCGA_models_to_remove_non_rand <- TCGA_models_to_remove_non_rand %>%
        mutate(ID = paste(NMD_method,dataset,cancer_type, sep = "_"))
TCGA_models_to_remove_rand <- lambdas_res[which(lambdas_res$lambda_random >= 1.5),]
TCGA_models_to_remove_rand <- TCGA_models_to_remove_rand %>%
        mutate(ID = paste(NMD_method,dataset,cancer_type, sep = "_"))

GWAS_final_df_filt <- GWAS_final_df %>%
    mutate(ID = paste(NMD_method,dataset,cancer_type,sep = "_")) %>%
    filter(! (ID %in% TCGA_models_to_remove_non_rand$ID ) ) %>%
    filter(! (ID %in% TCGA_models_to_remove_rand$ID ) ) 

GWAS_final_df_filt$cancer_type <- gsub("TCGA-","",GWAS_final_df_filt$cancer_type)

# 4.3) Number of replicated hits random vs non random  with stratification: only by tissue

NMD_methods <- c("ASE","Endogenous")
GWAS_replication_hits <- c()
for (geneset in c("Random","Selected")) {
    for (FDR in c(1,2,3,4,5,10)) {
    # for (FDR in c(25)) {
        print(FDR)
        for (NMD_method_discovery in NMD_methods) {
            print(NMD_method_discovery)
            for (NMD_method_validation in NMD_methods) {
                print(NMD_method_validation)
                if (NMD_method_discovery == NMD_method_validation) {next}
                    GWAS_replication_hits_tmp <- GWAS_replication_strat_cancer(GWAS_df_all = GWAS_final_df_filt, genesets_char = geneset,
                                                    FDR = FDR, NMD_method_discovery = NMD_method_discovery, NMD_method_validation = NMD_method_validation)
                    if (!is.null(GWAS_replication_hits_tmp)) {
                        GWAS_replication_hits_tmp$FDR_threshold_used <- FDR 
                        GWAS_replication_hits <- rbind(GWAS_replication_hits,GWAS_replication_hits_tmp)
                    }
            }
        }
    }
}

output_file <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/somatic_associations/analysis_plots/all_replicated_hits.txt")
# write.table(data.frame(GWAS_replication_hits), file = output_file, sep = "\t", quote = FALSE,
#                        col.names = TRUE, row.names = FALSE)
GWAS_replication_hits <- read.table(file = output_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
data.frame(GWAS_replication_hits[GWAS_replication_hits$Gene_symbol == "CDH1",])

GWAS_replication_hits <- GWAS_replication_hits %>%
  mutate(dataset = case_when(
    dataset == "som_mut_truncating" ~ "Truncating",
    dataset == "som_mut_missense" ~ "Missense",
    dataset == "som_mut_synonymous" ~ "Synonymous",
    dataset == "som_CNV_amp" ~ "CNA amp",
    dataset == "som_CNV_del" ~ "CNA del",
    TRUE ~ dataset
  ))

p <- GWAS_replication_hits %>%
        filter(!cancer_type %in% c("PCPG","CESC")) %>%
        group_by(FDR_threshold_used, NMD_method_discovery, NMD_method_validation, dataset, genesets) %>%
        summarise(num_replicated_hits = ( n() / genes_tested ) * 100) %>%
        mutate(NMD_method_validation = paste0(NMD_method_validation," - validation")) %>%
        mutate(NMD_method_discovery = paste0(NMD_method_discovery," - discovery")) %>% data.frame() %>%
            ggplot(aes(x = factor(FDR_threshold_used), y = num_replicated_hits, fill = genesets)) + 
                geom_bar(stat='identity',position=position_dodge()) +
                facet_grid(NMD_method_validation + NMD_method_discovery ~ dataset ) +
                #scale_fill_brewer(palette = "Set2", direction = -1, labels = c("Non Random", "Random")) +
                labs(title = "", x = "% FDR threshold", y = "% of Replicated Hits", fill = "") + 
                theme_bw(base_size = 30) + scale_x_discrete(labels=c("1","2","3","4","5","10")) +
                theme(axis.text.x = element_text(size = 35),
                    axis.title.x = element_text(size = 40),
                    strip.text = element_text(size = 30),
                    axis.text.y = element_text(size = 32),
                    axis.title.y = element_text(size = 40),
                    legend.position = "top",
                    legend.text = element_text(size = 35)) +
                guides(fill = guide_legend(override.aes = list(size = 14), nrow = 1))
output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/somatic_associations/analysis_plots/num_replicated_hits_strat_tissue.png")
png(output_path, width = 7000, height = 4500, res = 300)
print(p)
dev.off()

GWAS_replication_hits$NMD_method <- ifelse(GWAS_replication_hits$NMD_method == "Endogenous","ETG","ASE")
GWAS_replication_hits$NMD_method_discovery <- ifelse(GWAS_replication_hits$NMD_method_discovery == "Endogenous","ETG","ASE")
GWAS_replication_hits$NMD_method_validation <- ifelse(GWAS_replication_hits$NMD_method_validation == "Endogenous","ETG","ASE")

# Save
write.table(GWAS_replication_hits, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig7/SuppFig7C.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
saveRDS(GWAS_replication_hits, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig7/SuppFig7C.RData")

# 4.4) Table of gene-cancer replicated hits at 5% FDR
GWAS_replication_hits$gene_cancer <- paste0(GWAS_replication_hits$Gene_symbol," - ",GWAS_replication_hits$cancer_type)
GWAS_replication_hits$NMD_method_translation <- paste0(GWAS_replication_hits$NMD_method_discovery,"/",GWAS_replication_hits$NMD_method_validation)

# check the colours
cols <- brewer.pal(n = 5, name = "RdBu")

# Take into account that beta coefficients are inversed (When I did the associations I did nnot convert NMDeff high to positive values)
# So we will revert the values in the plots

# Selected genes
GWAS_replication_hits$disc_val <- ifelse(GWAS_replication_hits$disc_val == "discovery","disc","val")
GWAS_replication_hits_table <- GWAS_replication_hits %>%
            filter(FDR_threshold_used == 5 & genesets == "Selected") %>%
            mutate(facet_order = factor(paste0(NMD_method, " - ",disc_val),
                            levels = c("ASE - disc", "ETG - val",
                                        "ETG - disc", "ASE - val")))
p <- GWAS_replication_hits_table %>%
        filter(!cancer_type %in% c("PCPG","CESC")) %>%
        group_by(disc_val) %>%
            ggplot(aes(x = dataset, y = gene_cancer, fill = beta_coefficient)) +
            geom_tile(size = 2) +
            geom_text(aes(label = round(-beta_coefficient,2)), color = "black", size = 10) +
            theme_bw(base_size = 35) +
            facet_nested(. ~ facet_order) +
            theme(axis.text.x = element_text(size = 40, angle = 45, hjust=1),
                    axis.text.y = element_text(size = 40),
                    strip.text = element_text(size = 33),
                    legend.title=element_text(size=35)) +
            labs(fill = 'Beta Coefficient') + xlab("") + ylab("") +
            # scale_fill_gradientn(colours = c('grey','#F21A00'))
            scale_fill_gradientn(colours = cols, 
                                    values = rescale(c(-2, -1, 0, 1, 2)),
                                    guide = guide_colorbar(barheight = 15, barwidth = 2),
                                    limits=c(-2, 2))
            #guides(fill = guide_legend(override.aes = list(size = 20)))
output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/somatic_associations/analysis_plots/table_selected_genes_replicated.png")
png(output_path, width = 9000, height = 3000, res = 300)
print(p)
dev.off()

# Save
write.table(GWAS_replication_hits_table, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig7/SuppFig7D.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
saveRDS(GWAS_replication_hits_table, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig7/SuppFig7D.RData")

# Random Genes
GWAS_replication_hits_table <- GWAS_replication_hits %>%
            filter(FDR_threshold_used == 5 & genesets == "Random") %>%
            mutate(facet_order = factor(paste0(NMD_method, " - ",disc_val),
                            levels = c("ASE - disc", "ETG - val",
                                        "ETG - disc", "ASE - val")))

p <- GWAS_replication_hits_table %>%
        #filter(!cancer_type %in% c("TCGA-DLBC","TCGA-CESC")) %>%
        group_by(disc_val) %>%
            ggplot(aes(x = dataset, y = gene_cancer, fill = beta_coefficient)) +
            geom_tile(size = 2) +
            geom_text(aes(label = round(-beta_coefficient,2)), color = "black", size = 10) +
            theme_bw(base_size = 35) +
            facet_nested(. ~ facet_order) +
            theme(axis.text.x = element_text(size = 35, angle = 45, hjust=1),
                    axis.text.y = element_text(size = 35),
                    strip.text = element_text(size = 35),
                    panel.spacing = unit(3, "lines"),
                    legend.title=element_text(size=35)) +
            labs(fill = 'Beta Coefficient') + xlab("") + ylab("") +
            scale_fill_gradientn(colours = cols, 
                                    guide = guide_colorbar(barheight = 20, barwidth = 2),
                                    values = rescale(c(-4, -1, 0, 1, 4)),
                                    limits=c(-4, 4))
output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/somatic_associations/analysis_plots/table_random_genes_replicated.png")
png(output_path, width = 10000, height = 8000, res = 300)
print(p)
dev.off()

# 4.5) Beta Coefficients for all cancers in the interested hits
selected_hits <- c("TLX1","CDH1")

GWAS_selected_hits <- GWAS_final_df_filt %>%  
            filter(genesets == "Selected") %>%
            filter( (Gene_symbol %in% "CDH1" & dataset == "som_mut_truncating") |
            (Gene_symbol %in% "TLX1" & dataset == "som_mut_missense"))# |
            #(Gene_symbol %in% "PSMA1" & dataset == "som_CNV_amp") |
            #(Gene_symbol %in% "PTPN11" & dataset == "som_CNV_del"))

p <- GWAS_selected_hits %>%
            ggplot(aes(x = Gene_symbol, y = cancer_type, fill = beta_coefficient)) +
            geom_tile(size = 2) +
            geom_text(aes(label = round(-beta_coefficient,2)), color = "black", size = 13) +
            theme_bw(base_size = 35) +
            facet_nested(. ~ NMD_method) +
            theme(axis.text.x = element_text(size = 43),
                    axis.text.y = element_text(size = 40),
                    strip.text = element_text(size = 42),
                    panel.spacing = unit(3, "lines"),
                    #legend.position="top",
                    legend.title=element_text(size=35)) +
            geom_text(aes(label = ifelse(p_value < 0.02, "*", "")), 
                               size = 25, hjust = -2, color = "black") +
            labs(fill = 'Beta Coefficient') + xlab("") + ylab("") +
            scale_fill_gradientn(colours = cols, na.value = 'white',
                                    guide = guide_colorbar(barheight = 20, barwidth = 2),
                                    values = rescale(c(-2, -1, 0, 1, 2)),
                                    limits=c(-2, 2))
output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/somatic_associations/analysis_plots/selected_replicated_hits_all_cancers.png")
png(output_path, width = 10000, height = 8000, res = 300)
print(p)
dev.off()

# Save
write.table(GWAS_selected_hits, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig7/SuppFig7E.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
saveRDS(GWAS_selected_hits, "/g/strcombio/fsupek_home/gpalou/Manuscript/figures_data/SuppFig7/SuppFig7E.RData")
