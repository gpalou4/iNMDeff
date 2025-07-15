# 1) DATA
# 1.1) Open Cancer Gene Census from COSMIC database
CGC <- read.table("/g/strcombio/fsupek_cancer1/gpalou/COSMIC/cancer_gene_census.tsv", 
                    header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# 1.2) Open ENSEMBL v88
ensembl_gene_v88 <- read.table(file = "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/ensembl_v88_gene_transcript_genesymbol.txt", 
                                    header = TRUE, sep = "\t")
# 1.3) ENSEMBL transcripts IDs hg38 GTF
ensembl_v88_gtf <- rtracklayer::import("~/data/NMD_project/conversor_tables/ENSEMBL/gencode.v26.annotation.gtf")
ensembl_v88_gtf <- as.data.frame(ensembl_v88_gtf)
# Remove versions from IDs
ensembl_v88_gtf[,c("gene_id","transcript_id")] <- sapply(ensembl_v88_gtf[,c("gene_id","transcript_id")], function(col) { 
    gsub("(.*)\\..*","\\1", col)
})
ensembl_v88_gtf <- ensembl_v88_gtf[ensembl_v88_gtf$transcript_type %in% "protein_coding",]

# 2) Update ENSEMBL gene ID
# Obtain ENSEMBL gene ID
CGC$Synonyms <- gsub("\\w.*(ENSG.*[0-9]{11}).*","\\1",CGC$Synonyms)
# Check which genes do not have ENSEMBL gene ID and obtain its Gene Symbol
missing_ENSEMBL_index <- which(CGC$Synonyms == "")
missing_ENSEMBL <- CGC[missing_ENSEMBL_index,"Gene.Symbol"]
# There are 7 genes we do not have the ENSEMBL gene ID...
# Retrieve ENSEMBL gene ID from our ensembl_gene_v88 table and update the GCG
CGC_missing <- merge(CGC[missing_ENSEMBL_index,],ensembl_gene_v88[,c("gene_id","gene_name")], by.x = "Gene.Symbol", by.y="gene_name", all.x = TRUE)
CGC_missing <- CGC_missing[!duplicated(CGC_missing),]
CGC_missing$Synonyms <- CGC_missing$gene_id
CGC_missing$gene_id <- NULL
CGC_filt <- CGC[-missing_ENSEMBL_index,]
CGC_updated <- rbind(CGC_filt,CGC_missing)

# 3) Add NMD factors at CGC table
NMD_genes <- read.table(file = "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/NMD_genes.txt",
                header = TRUE, sep = "\t", stringsAsFactors = FALSE)
NMD_genes_filt <- NMD_genes[-which(NMD_genes$gene_symbol%in%CGC$Gene.Symbol),]
NMD_genes_df <- as.data.frame(matrix(nrow=nrow(NMD_genes_filt), ncol = length(colnames(CGC_updated))))
colnames(NMD_genes_df) <- colnames(CGC_updated)
NMD_genes_df$Gene.Symbol <- NMD_genes_filt$gene_symbol
NMD_genes_df$Synonyms <- NMD_genes_filt$ensembl_gene_id
# Add CGC
CGC_NMD <- rbind(CGC_updated,NMD_genes_df)
# Add NMD group
NMD_genes$gene_symbol_synonym <- NULL
CGC_NMD <- merge(CGC_NMD,NMD_genes, by.x = c("Gene.Symbol","Synonyms"), by.y = c("gene_symbol","ensembl_gene_id"), all.x = TRUE)

# 4) Add random set of genes as negative controls
# Coding
# 4.1) Random number from which we will select random genes
# n <- 5000
# random_genes <- sample(unique(ensembl_v88_gtf$gene_id),n)
# random_genes_filt <- random_genes[!random_genes%in%CGC_NMD$Synonyms]
# # 4.2) Remove those close to CGC genes
# ensembl_v88_gtf_filt <- ensembl_v88_gtf[ensembl_v88_gtf$gene_id %in% random_genes_filt,]
# ensembl_v88_gtf_CGC <- ensembl_v88_gtf[ensembl_v88_gtf$gene_id %in% CGC_NMD$Synonyms,]
# distance <- 50000
# overlaps <- c()
# overlapping_genes <- c()
# for (random_gene in random_genes_filt) {
#     print(random_gene)
#     overlaps <- c()
#     ensembl_v88_gtf_filt2 <- ensembl_v88_gtf_filt[ensembl_v88_gtf_filt$gene_id %in% random_gene,]
#     overlap_1 <- which(abs(ensembl_v88_gtf_CGC$start - min(unique(ensembl_v88_gtf_filt2$start))) < distance)
#     overlap_2 <- which(abs(ensembl_v88_gtf_CGC$start - max(unique(ensembl_v88_gtf_filt2$start))) < distance)
#     overlap_3 <- which(abs(ensembl_v88_gtf_CGC$start - min(unique(ensembl_v88_gtf_filt2$end))) < distance)
#     overlap_4 <- which(abs(ensembl_v88_gtf_CGC$start - max(unique(ensembl_v88_gtf_filt2$end))) < distance)
#     overlap_5 <- which(abs(ensembl_v88_gtf_CGC$end - min(unique(ensembl_v88_gtf_filt2$start))) < distance)
#     overlap_6 <- which(abs(ensembl_v88_gtf_CGC$end - max(unique(ensembl_v88_gtf_filt2$start))) < distance)
#     overlap_7 <- which(abs(ensembl_v88_gtf_CGC$end - min(unique(ensembl_v88_gtf_filt2$end))) < distance)
#     overlap_8 <- which(abs(ensembl_v88_gtf_CGC$end - max(unique(ensembl_v88_gtf_filt2$end))) < distance)
#     overlaps <- unique(c(overlaps,overlap_1,overlap_2,overlap_3,overlap_4,overlap_5,overlap_6,overlap_7,overlap_8))
#     if (length(overlaps >= 1)) {
#         overlapping_genes <- c(overlapping_genes,random_gene)
#     }
# }
# random_genes_final <- ensembl_v88_gtf_filt[!ensembl_v88_gtf_filt$gene_id %in% overlapping_genes,]
# random_genes_final <- random_genes_final[,c("gene_name","gene_id")]
# random_genes_final <- random_genes_final[!duplicated(random_genes_final),]
# number_random_genes <- nrow(CGC)
# n <- sample(1:nrow(random_genes_final),number_random_genes)
# random_genes_final <- random_genes_final[n,]

# Add all other genes into the CGC table
all_other_genes <- ensembl_v88_gtf[!ensembl_v88_gtf$gene_id %in% CGC_NMD$Synonyms, c("gene_name","gene_id")]
all_other_genes <- all_other_genes[!duplicated(all_other_genes),]

# Add to CGC
final_df <- matrix(nrow=nrow(all_other_genes), ncol=ncol(CGC_NMD))
colnames(final_df) <- colnames(CGC_NMD)
final_df <- as.data.frame(final_df)
final_df$Role.in.Cancer <- "non_CGC_NMD"
final_df[,c("Gene.Symbol","Synonyms")] <- all_other_genes
CGC_NMD_final <- rbind(CGC_NMD,final_df)

# 5) Add Gene Symbols for Immunoglobulin locus from CGC (FRAN SAYS NOT NECESSARY)
# Only those found in somatic CNV data
# IGH_genes <- unique(CNV_file[grep("^IGH.*",CNV_file$Gene.Symbol),"Gene.Symbol"])
# IGH_genes <- IGH_genes[-1]
# IGK_genes <- unique(CNV_file[grep("^IGK.*",CNV_file$Gene.Symbol),"Gene.Symbol"])
# IGL_genes <- unique(CNV_file[grep("^IGL.*",CNV_file$Gene.Symbol),"Gene.Symbol"])
# IGL_genes <- IGL_genes[-1]
# TRA_genes <- unique(CNV_file[grep("^TRA[CJV].*",CNV_file$Gene.Symbol),"Gene.Symbol"])
# TRB_genes <- unique(CNV_file[grep("^TRB.*",CNV_file$Gene.Symbol),"Gene.Symbol"])
# TRD_genes <- unique(CNV_file[grep("^TRD[DVJC].*",CNV_file$Gene.Symbol),"Gene.Symbol"])
# IG_TCR_genes <- c(IGH_genes,IGK_genes,IGL_genes,TRA_genes,TRB_genes,TRD_genes)

# 6) Change Gene Symbols to match somatic CNV data
# path <- "/g/strcombio/fsupek_cancer1/gpalou/TCGA_CNV_legacy/primary_tumor/gdac.broadinstitute.org_ESCA-TP.CopyNumber_Gistic2.Level_4.2016012800.0.0/all_thresholded.by_genes.txt"
# CNV_file <- read.table(file = path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# CNV_file[grep("CRLF",CNV_file$Gene.Symbol),1:3]
# CGC_NMD_final[!CGC_NMD_final$Gene.Symbol %in% CNV_file$Gene.Symbol,c("Gene.Symbol","Synonyms")]

CGC_NMD_final$Gene.Symbol <- as.character(CGC_NMD_final$Gene.Symbol)
CGC_NMD_final[CGC_NMD_final$Gene.Symbol %in% "AFDN","Gene.Symbol"] <- "MLLT4"
CGC_NMD_final[CGC_NMD_final$Gene.Symbol %in% "DUX4L1","Gene.Symbol"] <- "DUX4"
CGC_NMD_final[CGC_NMD_final$Gene.Symbol %in% "KNL1","Gene.Symbol"] <- "CASC5"
CGC_NMD_final[CGC_NMD_final$Gene.Symbol %in% "LHFPL6","Gene.Symbol"] <- "LHFP"
CGC_NMD_final[CGC_NMD_final$Gene.Symbol %in% "MRTFA","Gene.Symbol"] <- "MKL1"
CGC_NMD_final[CGC_NMD_final$Gene.Symbol %in% "NKX2-1","Gene.Symbol"] <- "TTF1"
CGC_NMD_final[CGC_NMD_final$Gene.Symbol %in% "NSD2","Gene.Symbol"] <- "WHSC1"
CGC_NMD_final[CGC_NMD_final$Gene.Symbol %in% "NSD3","Gene.Symbol"] <- "WHSC1L1"
CGC_NMD_final[CGC_NMD_final$Gene.Symbol %in% "SHTN1","Gene.Symbol"] <- "KIAA1598"
CGC_NMD_final[CGC_NMD_final$Gene.Symbol %in% "TENT5C","Gene.Symbol"] <- "FAM46C"
CGC_NMD_final[CGC_NMD_final$Gene.Symbol %in% "WDCP","Gene.Symbol"] <- "C2orf44"
CGC_NMD_final[CGC_NMD_final$Gene.Symbol %in% "ICE1","Gene.Symbol"] <- "KIAA0947"
CGC_NMD_final[CGC_NMD_final$Gene.Symbol %in% "PRELID3B","Gene.Symbol"] <- "SLMO2"
CGC_NMD_final[CGC_NMD_final$Gene.Symbol %in% "SF3B6","Gene.Symbol"] <- "SF3B14"
CGC_NMD_final[CGC_NMD_final$Gene.Symbol %in% "SNU13","Gene.Symbol"] <- "NHP2L1"
CGC_NMD_final[CGC_NMD_final$Gene.Symbol %in% "SEM1","Gene.Symbol"] <- "SHFM1"
CGC_NMD_final[CGC_NMD_final$Gene.Symbol %in% "RACK1","Gene.Symbol"] <- "GNB2L1"
CGC_NMD_final[CGC_NMD_final$Gene.Symbol %in% "SCX","Gene.Symbol"] <- "SCXA"

# Unkown genes symbols
# CGC_NMD_final[CGC_NMD_final$Gene.Symbol %in% "HLA-A","Gene.Symbol"]
# CGC_NMD_final[CGC_NMD_final$Gene.Symbol %in% "SNRPD2P1","Gene.Symbol"]
# CGC_NMD_final[CGC_NMD_final$Gene.Symbol %in% "SNRPEP2","Gene.Symbol"]
# CGC_NMD_final[CGC_NMD_final$Gene.Symbol %in% "RP11-402P6.15","Gene.Symbol"]
# CGC_NMD_final[CGC_NMD_final$Gene.Symbol %in% "AP003419.11","Gene.Symbol"]
# CGC_NMD_final[CGC_NMD_final$Gene.Symbol %in% "HNRNPCL4","Gene.Symbol"]

# Save
write.csv(CGC_NMD_final, file = "/g/strcombio/fsupek_cancer1/gpalou/COSMIC/cancer_gene_census_updated.tsv", quote = TRUE,
                      row.names = FALSE)













