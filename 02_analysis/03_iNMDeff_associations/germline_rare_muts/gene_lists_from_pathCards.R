# conda activate /home/gpalou/anaconda3_envs/general
library("queryPathcards")
library(dplyr)
library(Seurat)
library(httr)
library(rvest)
library(stringr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(tidyverse)
library("rbioapi")

## It's a modified version of the queryPathcard function
# It extracts Gene Lists from a given number of selected pathways from PathCards database
# Pathway Names must be exactly as it appear on the web page

GeneListfromPathCards <- function (pathways)  {
    url_template_pathway <- "https://pathcards.genecards.org/card/"

    if (identical(which(pathways %in% "NA"), integer(0)) == FALSE) {
        pathways <- pathways[-which(pathways %in% "NA")]
    }
    pathway_list <- vector("list", length(pathways))
    names(pathway_list) <- pathways
    pathways <- tolower(pathways)
    for (i in 1:length(pathways)) {
        name_simplified <- fixNames(pathways[i])
        url_3 <- paste(url_template_pathway, name_simplified, sep = "")
        #print(url_3)
        gene_data <- url_3 %>% read_html() %>% html_nodes(".geneslist") %>% 
            html_text() %>% str_remove_all("\r\n") %>% str_remove_all("via the multiplicity of each gene in the constituent pathways.") %>% 
            str_remove_all("\\*Darkness represents the genes rank within the SuperPath,") %>% 
            trimws() %>% as.data.frame()
        if (nrow(gene_data) >= 2) {
            gene_data <- gene_data[-1, ]
        } else {
            gene_data <- gene_data[1, ]
        }
        pathway_list[i] <- gene_data %>% as.character() %>% 
            strsplit("\\s+")
        print(as.character(pathways[i]))
    }
    return(pathway_list)
}

# 1) Obtain genes from PathCards database
query_pathways <- c("Processing of Capped Intron-Containing Pre-mRNA",
                    "Unfolded Protein Response (UPR)",
                    "Processing of Capped Intronless Pre-mRNA",
                    "Translational Control",
                    "Transport of Mature Transcript to Cytoplasm",
                    "Gene expression (Transcription)",
                    "RNA Polymerase III Transcription Initiation",
                    "RNA Polymerase I Transcription Termination",
                    "RNA Polymerase I Promoter Opening",
                    "RNA Polymerase II Transcription Initiation And Promoter Clearance",
                    "Deadenylation-dependent mRNA decay",
                    "Chromatin Regulation / Acetylation",
                    "Chromatin organization")
pathway_genes_list <- GeneListfromPathCards(pathways = query_pathways)
lapply(pathway_genes_list, function(pathway) {length(pathway)})
# Convert the list to a dataframe
pathway_genes_df <- tibble(Pathways = names(pathway_genes_list), Genes = pathway_genes_list) %>%
                    unnest(Genes) %>% as.data.frame()
# Print the dataframe
print(head(pathway_genes_df))

# 2) NMD genes from experimental papers
NMD_genes <- read.table(file = "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/NMD_genes.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
NMD_genes_filt <- NMD_genes[,c(4,1)]
colnames(NMD_genes_filt) <- c("Pathways","Genes")
pathway_genes_all <- rbind(NMD_genes_filt,pathway_genes_df)

# 3) NMD targets
NMD_targets <- read.table( file = "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/NMD_targets/ENSEMBL_transcripts/NMD_global_ensembl.txt",
                 header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 4) Check 

# 4.1) Overlap between genes of different pathways

# Create a matrix to store the overlapped gene counts
overlap_matrix <- matrix(0, nrow = length(pathway_genes_list), ncol = length(pathway_genes_list), dimnames = list(names(pathway_genes_list), names(pathway_genes_list)))
# Iterate over each combination of gene pathways
for (i in 1:length(pathway_genes_list)) {
  for (j in 1:length(pathway_genes_list)) {
    if (i != j) {
      genes_i <- pathway_genes_list[[i]]
      genes_j <- pathway_genes_list[[j]]
      #overlap <- length(intersect(genes_i, genes_j))
      # In %
      overlap <- round(length(intersect(genes_i, genes_j)) / max(length(genes_i), length(genes_j)),2)
      overlap_matrix[i, j] <- overlap
    }
  }
}
# Print the overlap matrix
print(overlap_matrix)

# 4.2) Missing genes
# Missing genes are those not matching with our old gene set of ensembl v88, because the pathCards are current updated genes
# We will get the "old genes" from those missing ones
ensembl_v88_gene_transcript_genesymbol <- read.table(file = "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/ensembl_v88_gene_transcript_genesymbol.txt",
                                                    header = TRUE, sep = "\t", stringsAsFactors = FALSE)
ensembl_v88_gene_symbol_updated <- read.csv(file = "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/ensembl_v88_gene_symbol_updated.csv")
ensembl_v88_gene_symbol_updated <- ensembl_v88_gene_symbol_updated[,c("Input","Approved.symbol")]
table(unique(pathway_genes_all$Genes) %in% ensembl_v88_gene_transcript_genesymbol$gene_name)
# table(missing_genes %in% ensembl_v88_gene_symbol_updated$Approved.symbol)
# missing_genes[!missing_genes %in% ensembl_v88_gene_symbol_updated$Approved.symbol]
missing_genes_df <- pathway_genes_all[!pathway_genes_all$Genes %in% ensembl_v88_gene_transcript_genesymbol$gene_name,]
old_genes_df <- ensembl_v88_gene_symbol_updated[ensembl_v88_gene_symbol_updated$Approved.symbol %in% missing_genes_df$Genes,]
recovered_genes_df <- merge(missing_genes_df,old_genes_df, by.x = "Genes", by.y = "Approved.symbol", all.x = TRUE)
# Only 22 genes are not recovered. Remove them 
table(is.na(recovered_genes_df$Input))
recovered_genes_df$Genes <- recovered_genes_df$Input
recovered_genes_df$Input <- NULL
recovered_genes_df <- na.omit(recovered_genes_df)
# Merge
matching_genes_df <- pathway_genes_all[pathway_genes_all$Genes %in% ensembl_v88_gene_transcript_genesymbol$gene_name,]
pathway_updated_genes <- rbind(matching_genes_df,recovered_genes_df)
table(pathway_updated_genes$Genes %in% ensembl_v88_gene_transcript_genesymbol$gene_name)

# Save NMD-related pathway
NMD_related_pathways <- c("Processing of Capped Intron-Containing Pre-mRNA",
                        "Processing of Capped Intronless Pre-mRNA",
                        "Translational Control",
                        "Transport of Mature Transcript to Cytoplasm"
                        )
NMD_pathways <- c("DECID_complex","EJC","EJC_associated","NMD_ER","NMD_factor","NMD_related")
pathway_updated_genes_NMD_related <- pathway_updated_genes %>% 
                        filter(Pathways %in% c(NMD_related_pathways,NMD_pathways))

output_path <- "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/pathCards_pathways_NMD_related.txt"
write.table(pathway_updated_genes_NMD_related, file = output_path,
            quote = FALSE, sep='\t', row.names = FALSE, col.names = TRUE)

# Save Gene Regulation pathway
output_path <- "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/pathCards_pathways_mRNA_regulation.txt"
write.table(pathway_updated_genes, file = output_path,
            quote = FALSE, sep='\t', row.names = FALSE, col.names = TRUE)

# 4.3) STRING expand nodes from the previous NMD gene lists
#STRING_info <- read.table(file = "/g/strcombio/fsupek_cancer1/gpalou/STRING/protein.physical.links.v11.5.score_filt_500.txt", header = TRUE, sep = "\t")

pathway_updated_genes_NMD_related_expanded <- rba_string_interactions_network(
  pathway_updated_genes_NMD_related$Genes,
  species = 9606,
  required_score = 500,
  add_nodes = 10000,
  network_type = "physical"
)

# NMD confident genes only (experimental)
pathway_genes_all_NMD <- pathway_updated_genes_NMD_related[pathway_updated_genes_NMD_related$Pathways %in% NMD_pathways,]
pathway_updated_genes_NMD_confident_expanded <- rba_string_interactions_network(
  pathway_genes_all_NMD$Genes,
  species = 9606,
  required_score = 500,
  add_nodes = 100,
  network_type = "physical"
)

genes <- unique(c(pathway_updated_genes_NMD_confident_expanded$preferredName_A,pathway_updated_genes_NMD_confident_expanded$preferredName_B))
NMD_genes_tmp <- data.frame(Genes = genes)
NMD_genes_confident_expanded <- merge(NMD_genes_tmp,NMD_genes, by.x = "Genes", by.y = "gene_symbol", all.x = TRUE)

# Save
output_path <- "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/NMD_genes_STRING_expanded.txt"
write.table(NMD_genes_confident_expanded, file = output_path,
            quote = FALSE, sep='\t', row.names = FALSE, col.names = TRUE)

# Update Missing genes (I can't recover any)
# table(unique(pathway_updated_genes_NMD_related_expanded$preferredName_A) %in% ensembl_v88_gene_transcript_genesymbol$gene_name)
# table(unique(pathway_updated_genes_NMD_related_expanded$preferredName_B) %in% ensembl_v88_gene_transcript_genesymbol$gene_name)
# update_genes <- function(genes_df, column) {
#     missing_genes_df <- genes_df[!genes_df[[column]] %in% ensembl_v88_gene_transcript_genesymbol$gene_name,]
#     old_genes_df <- ensembl_v88_gene_symbol_updated[ensembl_v88_gene_symbol_updated$Approved.symbol %in% missing_genes_df[[column]],]
#     recovered_genes_df <- merge(missing_genes_df,old_genes_df, by.x = column, by.y = "Approved.symbol", all.x = TRUE)
#     table(is.na(recovered_genes_df$Input))
#     recovered_genes_df[[column]] <- recovered_genes_df$Input
#     recovered_genes_df$Input <- NULL
#     recovered_genes_df <- na.omit(recovered_genes_df)
#     # Merge
#     matching_genes_df <- genes_df[genes_df[[column]] %in% ensembl_v88_gene_transcript_genesymbol$gene_name,]
#     pathway_updated_genes <- rbind(matching_genes_df,recovered_genes_df)
#     table(pathway_updated_genes[[column]] %in% ensembl_v88_gene_transcript_genesymbol$gene_name)
#     return(pathway_updated_genes)
# }
# pathway_updated_genes_NMD_expanded_A <- update_genes(genes_df = pathway_updated_genes_NMD_related_expanded, column = "preferredName_A")
# pathway_updated_genes_NMD_expanded_B <- update_genes(genes_df = pathway_updated_genes_NMD_related_expanded, column = "preferredName_B")

# 5) Random genes
# Check random genes that contain rare germline variants
# GTEx
input_path <- paste0("/g/strcombio/fsupek_cancer1/gpalou/GTEx/WES_VCF/rare_germline_variants/GTEx_germline_input_variants_PTV_Missense_CADD15_0.1perc.txt")
GTEx_all_rare_variants <- read.csv(file = input_path, head=T, sep ="\t",stringsAsFactors = F)
samples_outlier_path <- paste0("/g/strcombio/fsupek_cancer1/gpalou/GTEx/samples_metadata/ancestry_and_admixture_and_PCA_outliers.txt")
all_outlier_samples <- read.table(file = samples_outlier_path)$V1
GTEx_all_rare_variants_filt <- GTEx_all_rare_variants %>%
                          filter(!GTEx_sample %in% all_outlier_samples)
# TCGA
input_path <- paste0("/g/strcombio/fsupek_cancer1/gpalou/Mischan/TCGA_germline_input_variants_PTV_Missense_CADD15_0.1perc.txt")
TCGA_all_rare_variants <- read.csv(file = input_path, head=T, sep ="\t",stringsAsFactors = F)
# Cancer genes
input_path <- "/g/strcombio/fsupek_cancer1/gpalou/COSMIC/cancer_gene_census_updated.tsv"
CGC <- read.csv(file = input_path, head=T,sep =",",stringsAsFactors = F)
cancer_genes <- CGC[which(CGC$Role.in.Cancer != "non_CGC_NMD"),"Gene.Symbol"]
# LOEUF score
LOEUF <- read.table(file = "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/LOEUF_score.tsv",
            header = TRUE, sep = "\t")
LOEUF_filt <- LOEUF %>% 
    group_by(gene) %>%
    summarise(median(oe_lof))
colnames(LOEUF_filt)[2] <- "oe_lof"
# Essential genes
essential <- read.table(file = "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/essential_CSEGs_CEGs.txt",
            header = TRUE, sep = "\t", skip = 11)
essential <- essential[essential$essentiality == "CEGs",]

# Filter genes from
# gene transcript regulation network // cancer genes // NMD_related // NMD_related expanded network
random_genes <- ensembl_v88_gene_transcript_genesymbol %>% 
                filter(!gene_name %in% pathway_updated_genes$Genes) %>%
                filter(!gene_name %in% cancer_genes) %>%
                filter(!gene_name %in% pathway_updated_genes_NMD_related_expanded$preferredName_A) %>%
                filter(!gene_name %in% pathway_updated_genes_NMD_related_expanded$preferredName_B) %>%
                filter(!gene_name %in% NMD_targets$gene_symbol) %>%
                pull(gene_name) %>% unique()
GTEx_freq_variants <- data.frame(table(GTEx_all_rare_variants_filt$gene_name))
GTEx_genes_with_variants <- unique(GTEx_freq_variants[GTEx_freq_variants$Freq >= 2,"Var1"])
GTEx_random_genes_with_variants <- as.character(GTEx_genes_with_variants[GTEx_genes_with_variants %in% random_genes])
TCGA_freq_variants <- data.frame(table(TCGA_all_rare_variants$Gene.refGene))
TCGA_genes_with_variants <- unique(TCGA_freq_variants[TCGA_freq_variants$Freq >= 2,"Var1"])
TCGA_random_genes_with_variants <- as.character(TCGA_genes_with_variants[TCGA_genes_with_variants %in% random_genes])
random_genes_with_variants <- intersect(GTEx_random_genes_with_variants,TCGA_random_genes_with_variants)

# X genes with variants in TCGA and, from those, X with variants in GTEx - matched sizes with the NMD list
GTEx_n <- sum(unique(NMD_genes_confident_expanded$Genes) %in% GTEx_genes_with_variants)
TCGA_n <- sum(unique(NMD_genes_confident_expanded$Genes) %in% TCGA_genes_with_variants)

############ SAMPLING MATCHING BY LOEUF SCORE & ESSENTIAL GENES PROPORTION ###########
# Sample the same number of genes that have rare variants in our NMD_related gene list
# Match so that both genesets have a similar LOEUF score
# Median of LOEUF score for NMD genes
LOEUF_NMD_median <- LOEUF_filt[LOEUF_filt$gene %in% NMD_genes_confident_expanded$Genes,"oe_lof"] %>% 
                    summarise(median(oe_lof,na.rm = TRUE)) %>% as.numeric()
# Filter gene with LOEUF score similar to NMD genes
LOEUF_filt2 <- LOEUF_filt[LOEUF_filt$oe_lof < (LOEUF_NMD_median+0.15),]
TCGA_random_genes_with_variants_filt <- TCGA_random_genes_with_variants[TCGA_random_genes_with_variants %in% LOEUF_filt2$gene]
# Proportion of essential genes of NMD genes
NMD_genes_mutated <- unique(NMD_genes_confident_expanded$Genes)[unique(NMD_genes_confident_expanded$Genes) %in% TCGA_genes_with_variants]
NMD_genes_essential_prop <- round(as.numeric(table(NMD_genes_mutated %in% essential$gene)/length(NMD_genes_mutated))[2],2)

# Create a function to check if a set of genes meets the criteria
check_criteria <- function(sampled_genes) {
  # Find the LOEUF scores of the sampled genes
  LOEUF_median <- LOEUF_filt[LOEUF_filt$gene %in% sampled_genes,"oe_lof"] %>% 
                    summarise(median(oe_lof,na.rm = TRUE)) %>% as.numeric()
  print(LOEUF_median)
  criteria_1 <- LOEUF_median > (LOEUF_NMD_median-0.005) & LOEUF_median < (LOEUF_NMD_median+0.005)
  # Find the proportion of essential genes of the samples genes
  # random_genes_essential_prop <- round(as.numeric(table(sampled_genes %in% essential$gene)/length(sampled_genes))[2],2)
  # print(random_genes_essential_prop)
  # criteria_2 <- random_genes_essential_prop > (NMD_genes_essential_prop-0.005) & random_genes_essential_prop < (NMD_genes_essential_prop+0.005) 
  # criteria_2 <- TRUE
  # Check if the criterias are tre
  if (isTRUE(criteria_1)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}
# Randomly sample genes until a set that meets the criterias is found
set.seed(333)
non_essential_random_genes <- TCGA_random_genes_with_variants_filt[!TCGA_random_genes_with_variants_filt %in% essential$gene]
essential_random_genes <- TCGA_random_genes_with_variants_filt[TCGA_random_genes_with_variants_filt %in% essential$gene]
sampled_genes <- NULL
while (is.null(sampled_genes)) {
  non_essential_random_genes_reduced <- sample(non_essential_random_genes,size=length(essential_random_genes)*15)
  non_essential_random_genes_reduced <- non_essential_random_genes_reduced[!non_essential_random_genes_reduced %in% "SLC6A1"]
  #sampled_genes <- sample(non_essential_random_genes, size = TCGA_n-length(essential_random_genes))
  sampled_genes <- sample(c(essential_random_genes,non_essential_random_genes_reduced), size = TCGA_n)
  if (!check_criteria(sampled_genes)) {
    sampled_genes <- NULL
  }
}

#final_gene_set_random <- c(sampled_genes,essential_random_genes)
final_gene_set_random <- c(sampled_genes)
GTEx_size_random <- table(final_gene_set_random %in% GTEx_random_genes_with_variants)[2]
# Match the same sample size in GTEx as our NMD_related gene list
# Remove exceed genes
GTEx_genes_to_remove <- final_gene_set_random[final_gene_set_random %in% GTEx_random_genes_with_variants]
GTEx_size_remove <- as.numeric(GTEx_size_random) - GTEx_n
final_gene_set_random <- final_gene_set_random[!final_gene_set_random %in% sample(GTEx_genes_to_remove,GTEx_size_remove)]
# Add genes
# TCGA genes not in GTEx
random_genes_with_variants_inTCGA_noGTEx <- TCGA_random_genes_with_variants[!(TCGA_random_genes_with_variants %in% GTEx_random_genes_with_variants)]
random_genes_with_variants_inTCGA_noGTEx <- random_genes_with_variants_inTCGA_noGTEx[!(random_genes_with_variants_inTCGA_noGTEx %in% final_gene_set_random)]
final_gene_set_random <- c(final_gene_set_random,sample(random_genes_with_variants_inTCGA_noGTEx,GTEx_size_remove))
# Check
table(final_gene_set_random %in% GTEx_random_genes_with_variants)
table(final_gene_set_random %in% TCGA_random_genes_with_variants)
final_gene_set_random <- data.frame(Genes = final_gene_set_random)

output_path <- "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/random_genes_not_cancer_not_NMD_related_expanded.txt"
write.table(final_gene_set_random, file = output_path,
            quote = FALSE, sep='\t', row.names = FALSE, col.names = TRUE)


