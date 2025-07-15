library("dplyr")

genes_update <- function(germline_variants_final_df, ensembl_v88_gene_transcript_genesymbol) {
    # Add gene symbols / ensembl gene id to germline_variants_final_df
    ensembl_gene_id_1 <- unlist(lapply(strsplit(germline_variants_final_df$Gene.refGene,";"), function(x){
        x[1]
    }))
    ensembl_gene_id_2 <- unlist(lapply(strsplit(germline_variants_final_df$Gene.refGene,";"), function(x){
        x[2]
    }))
    germline_variants_final_df$ensembl_gene_id_1 <- ensembl_gene_id_1
    germline_variants_final_df$ensembl_gene_id_2 <- ensembl_gene_id_2
    # Update ENSEMBLE GENE IDs (90%)
    missing_rows <- which(is.na(germline_variants_final_df$ensembl_gene_id_2))
    matching_genes <- ensembl_v88_gene_transcript_genesymbol[ensembl_v88_gene_transcript_genesymbol$gene_id %in% c(ensembl_gene_id_1,ensembl_gene_id_2),c("gene_id","gene_name")]
    matching_genes <- matching_genes[!duplicated(matching_genes),]
    germline_variants_final_df_1 <- germline_variants_final_df[missing_rows,]
    germline_variants_final_df_1 <- merge(germline_variants_final_df_1,matching_genes, by.x = "ensembl_gene_id_1", by.y = "gene_id", all.x = TRUE)
    germline_variants_final_df_2 <- germline_variants_final_df[-missing_rows,]
    germline_variants_final_df_2 <- merge(germline_variants_final_df_2,matching_genes, by.x = "ensembl_gene_id_2", by.y = "gene_id", all.x = TRUE)
    germline_variants_final_df_updated <- rbind(germline_variants_final_df_1,germline_variants_final_df_2)
    germline_variants_final_df_updated$gene_name.y <- NULL
    colnames(germline_variants_final_df_updated)[colnames(germline_variants_final_df_updated) == "gene_name.x"] <- "gene_name"
    germline_variants_final_df_updated <- germline_variants_final_df_updated %>% 
                                mutate(ensembl_gene_id = coalesce(ensembl_gene_id_1,ensembl_gene_id_2))
    return(germline_variants_final_df_updated)
}

# 1) Merge all variants files from all GTEx samples into one dataframe
GTEx_samples_variants_path <- paste0("/g/strcombio/fsupek_cancer1/gpalou/GTEx/WES_VCF/rare_germline_variants/samples/[X1]_rare_germline_variants.txt")
GTEx_samples <- read.table(file = "/g/strcombio/fsupek_cancer1/gpalou/GTEx/WES_VCF/list_samples.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)$V1
germline_variants_final_df <- c()
for (GTEx_sample in GTEx_samples) {
    print(GTEx_sample)
    GTEx_sample_variants_path <- gsub("\\[X1\\]",GTEx_sample,GTEx_samples_variants_path)
    # Check if variants file exists
    if (file.exists(GTEx_sample_variants_path)) {
        # Open sample CGC file
        GTEx_sample_variants <- read.csv(file = GTEx_sample_variants_path, header = TRUE, sep = "\t")
    } else {
        print("No germline variants file for the sample")
        next
    }
    # Add sample barcode
    GTEx_sample_variants$GTEx_sample <- GTEx_sample
    # Merge
    if ( length(germline_variants_final_df) == 0 ) {
        germline_variants_final_df <- GTEx_sample_variants
    } else {
        germline_variants_final_df <- rbind(germline_variants_final_df, GTEx_sample_variants)
    }
}
print(dim(germline_variants_final_df))

# 2) Add new columns to the final dataframe
conversor_tables_path <- "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/"
# Genotype
germline_variants_final_df$gt <- unlist(lapply(strsplit(germline_variants_final_df[,"Otherinfo13"],":"),function(x){x[[1]]}))
# Gene Symbols
ensembl_v88_gene_transcript_genesymbol <- read.table(file = paste0(conversor_tables_path,"ensembl_v88_gene_transcript_genesymbol.txt"), 
                                                    header = TRUE, sep = "\t", stringsAsFactors = FALSE)
germline_variants_final_df <- genes_update(germline_variants_final_df, ensembl_v88_gene_transcript_genesymbol)

# VARITY scores (better to use CADD scores because Mischan used them, replication...)

# varity_scores <- read.table(file = paste0(conversor_tables_path,"varity_predictions.txt"), 
#                             header = FALSE, stringsAsFactors = FALSE)
# colnames(varity_scores) <- c("chr","start_pos","varity_score","ref","alt")
# # to get rid of mapping to multiple VARITY score
# varity_scores_fixed <- varity_scores %>%
#   group_by(chr, start_pos, ref, alt) %>%
#   summarise(varity_score = sum(varity_score)/n())
# germline_variants_final_df <- merge(germline_variants_final_df, varity_scores, by.x = c("Chr","Start","Otherinfo7","Otherinfo8"), by.y = c("chr","start_pos","ref","alt"), all.x = TRUE)
# for (mutation_type in names(table(germline_variants_final_df$ExonicFunc.refGene))) {
#     print(mutation_type)
#     varity_scores_mutation <- germline_variants_final_df[germline_variants_final_df$ExonicFunc.refGene %in% mutation_type,"varity_score.y"]
#     print(length(varity_scores_mutation))
#     print(summary(varity_scores_mutation))
# }
# CADD scores
CADD_scores <- read.table(file = "/g/strcombio/fsupek_cancer1/gpalou/GTEx/WES_VCF/rare_germline_variants/CADD_scores/all_samples_liftOver_annotated_rare_variants_CADD.tsv.gz",
                            header = TRUE, sep = "\t", stringsAsFactors = FALSE, skip = 1, comment.char = "")
CADD_scores$X.Chrom <- paste0("chr",CADD_scores$X.Chrom)

germline_variants_final_df_final <- merge(germline_variants_final_df, CADD_scores, by.x = c("Chr","Start","Otherinfo7","Otherinfo8","ensembl_gene_id"), 
                                            by.y = c("X.Chrom","Pos","Ref","Alt","GeneID"), all.x = TRUE)
germline_variants_final_df_final$CADD_phred <- as.numeric(germline_variants_final_df_final$CADD_phred)
table(germline_variants_final_df_final$CADD_phred == germline_variants_final_df_final$PHRED)

# NMD rules (PTC evading vs triggering)
# GTEx_all_PTCs_NMDeff <- read.table(file = "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/PTCs_dataset/germline_PTCs_all_GTEx.txt", 
#                                     header = TRUE, sep = "\t")

for (mutation_type in names(table(germline_variants_final_df_final$ExonicFunc.refGene))) {
    print(mutation_type)
    CADD_scores_mutation <- germline_variants_final_df_final[germline_variants_final_df_final$ExonicFunc.refGene %in% mutation_type,"PHRED"]
    print(length(CADD_scores_mutation))
    print(summary(CADD_scores_mutation))
    CADD_scores_mutation <- as.numeric(germline_variants_final_df[germline_variants_final_df$ExonicFunc.refGene %in% mutation_type,"CADD_phred"])
    print(length(CADD_scores_mutation))
    print(summary(CADD_scores_mutation))
}

# Other of potential interest (Mischan columns)
#ref_depth,alt_depth,GQX,CADD,splicing_altering_effect,chr_lastExon,start_lastExon,end_lastExon,exonNo,
#gendeID,gene_length_transcribed,gene_symbol,variant_in_last_exon

# 3) Save raw dataset
output_path <- paste0("/g/strcombio/fsupek_cancer1/gpalou/GTEx/WES_VCF/rare_germline_variants/GTEx_all_rare_germline_variants.txt")
write.table(germline_variants_final_df_final, file = output_path, sep = "\t", quote = FALSE,
                      col.names = TRUE, row.names = FALSE)

# By Mischan paper --> 
# "PTVs comprised in this study frameshift deletions, frameshift insertions, stoploss variants, stopgain variants, 
# startloss variants and splicing variants. Splicing variants comprise the canonical splice variants annotated by 
# ANNOVAR131 (version 2019-10-24) and variants with a predicted donor loss or acceptor loss >0.8 by SpliceAI"

# 4) Create three datasets:
splicing_mut <- c("synonymous SNV","unknown",".")
PTV_mut <- c("frameshift deletion","frameshift insertion", "stopgain", "stoploss", "startloss")
# 4.1) #PTV_0.1perc
PTV_0.1perc <- germline_variants_final_df_final %>% 
                filter(( (ExonicFunc.refGene %in% splicing_mut) & (SpliceAI.acc.loss > 0.8 | SpliceAI.don.loss > 0.8) ) |
                        (ExonicFunc.refGene %in% PTV_mut))
output_path <- paste0("/g/strcombio/fsupek_cancer1/gpalou/GTEx/WES_VCF/rare_germline_variants/GTEx_germline_input_variants_PTV_0.1perc.txt")
write.table(PTV_0.1perc, file = output_path, sep = "\t", quote = FALSE,
                      col.names = TRUE, row.names = FALSE)
# 4.2) #PTV_Missense_CADD15_0.1perc
PTV_Missense_CADD15_0.1perc <- germline_variants_final_df_final %>% 
                filter(( (ExonicFunc.refGene %in% splicing_mut) & (SpliceAI.acc.loss > 0.8 | SpliceAI.don.loss > 0.8) ) |
                        (ExonicFunc.refGene %in% PTV_mut) |
                        ( (ExonicFunc.refGene %in% "nonsynonymous SNV") & (PHRED >= 15) )
                        )
output_path <- paste0("/g/strcombio/fsupek_cancer1/gpalou/GTEx/WES_VCF/rare_germline_variants/GTEx_germline_input_variants_PTV_Missense_CADD15_0.1perc.txt")
write.table(PTV_Missense_CADD15_0.1perc, file = output_path, sep = "\t", quote = FALSE,
                      col.names = TRUE, row.names = FALSE)
# 4.3) #PTV_Missense_CADD25_0.1perc
PTV_Missense_CADD25_0.1perc <- germline_variants_final_df_final %>% 
                filter(( (ExonicFunc.refGene %in% splicing_mut) & (SpliceAI.acc.loss > 0.8 | SpliceAI.don.loss > 0.8) ) |
                        (ExonicFunc.refGene %in% PTV_mut) |
                        ( (ExonicFunc.refGene %in% "nonsynonymous SNV") & (PHRED >= 25) )
                        )
output_path <- paste0("/g/strcombio/fsupek_cancer1/gpalou/GTEx/WES_VCF/rare_germline_variants/GTEx_germline_input_variants_PTV_Missense_CADD25_0.1perc.txt")
write.table(PTV_Missense_CADD25_0.1perc, file = output_path, sep = "\t", quote = FALSE,
                      col.names = TRUE, row.names = FALSE)




