# # facet_nested
library("ggh4x")
# library(cowplot)
# library(ggrepel)
library(scales)
library(biomaRt)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(plyr)

# 1) Data
# 1.1) GWAS somatic data
input_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/somatic_associations/final_associations/GWAS_somatic_mut_CNV.txt")
GWAS_final_df <- read.table(input_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# 1.2 ) Cytogenetic locations of genes
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
anno_gene <- biomaRt::getBM(attributes = c("ensembl_gene_id","external_gene_name","chromosome_name","start_position","end_position","band", "gene_biotype"),mart=ensembl, useCache = FALSE)

# 2) Check coefficients of genes of interest

# 2.1) CNA-PC3 and CNA-PC86 --> Chr 1q amplification in q21-23.1
GWAS_final_df_filt <- GWAS_final_df[GWAS_final_df$chr == 1,]
GWAS_final_df_filt <- merge(GWAS_final_df_filt, anno_gene[,c("ensembl_gene_id","band")], by.x = "ENSEMBL_gene", by.y = "ensembl_gene_id", all.x = TRUE)
#sig_band <- c("q21.1","q21.2","q21.3","q22","q23.1")
#sig_band <- c("q21.1","q22")
sig_band <- c("q21.1","q21.2","q21.3","q22","q23.1","q23.2","q23.3","q24.1","q24.2","q24.3","q25.1","q25.2",
                "q25.3","q31.1","q31.2","q31.3")
GWAS_final_df_filt2 <- GWAS_final_df_filt[GWAS_final_df_filt$band %in% sig_band,]
GWAS_final_df_filt2 <- GWAS_final_df_filt
selected_hits <- c("SMG5","SMG7","RBM8A","SF3B4","INTS3","SLC25A44","SEMA4A","SSR2","MEX3A","LAMTOR2","RAB25")
#selected_hits <- c("SMG5","SMG7","RBM8A","SF3B4","INTS3","SNORA26","SLC25A44","SEMA4A","SSR2","MEX3A","LAMTOR2","RAB25")
GWAS_selected_hits <- GWAS_final_df_filt2 %>%  
            #filter(genesets == "Selected") %>%
            #filter(Gene_symbol %in% selected_hits) %>%
            filter( dataset == "som_CNV_amp")

# ENSG00000212187 --> "SNORA26"
order_genome_location_by_chr <- function(genes_chr) {
  genes_chr$genome_location_chr <- as.numeric(gsub(".*:(.*)-.*","\\1",genes_chr$genome_location))
  genes_chr$chr <- as.numeric(genes_chr$chr)
#   df <- genes_chr %>%
#         group_by(NMD_method,cancer_type) %>%
#         arrange(chr,genome_location_chr)
  genes_chr <- genes_chr[order(genes_chr$chr, genes_chr$genome_location_chr, decreasing = FALSE),]
}

# remove extreme outliers

GWAS_selected_hits <- GWAS_selected_hits %>%  
            filter(NMD_method == "Endogenous")

GWAS_selected_hits <- GWAS_selected_hits[order(GWAS_selected_hits$cancer_type),]

GWAS_selected_hits_matrix <- data.frame(unstack(GWAS_selected_hits[,c("beta_coefficient","Gene_symbol")]))
GWAS_selected_hits_matrix <- GWAS_selected_hits_matrix[1:33,]
rownames(GWAS_selected_hits_matrix) <- unique(GWAS_selected_hits$cancer_type)
GWAS_selected_hits_matrix <- t(GWAS_selected_hits_matrix)
# Order by chr location
tmp <- unique(GWAS_selected_hits[!GWAS_selected_hits$Gene_symbol %in% c("RGS5","BTBD8"),c("Gene_symbol","genome_location","chr","band")])
GWAS_selected_hits_matrix <- merge(GWAS_selected_hits_matrix, tmp, 
                                    by.x = "row.names",by.y = "Gene_symbol", all.x = TRUE)
rownames(GWAS_selected_hits_matrix) <- GWAS_selected_hits_matrix$Row.names
GWAS_selected_hits_matrix_ordered <- order_genome_location_by_chr(GWAS_selected_hits_matrix)
GWAS_selected_hits_matrix_ordered$Row.names <- NULL
cols <- grep("TCGA-|pancancer",colnames(GWAS_selected_hits_matrix_ordered))
GWAS_selected_hits_matrix_ordered <- GWAS_selected_hits_matrix_ordered[!is.na(GWAS_selected_hits_matrix_ordered$genome_location),]

# Create a color for each band
GWAS_selected_hits_matrix_ordered$arm <- substr(GWAS_selected_hits_matrix_ordered$band,1,1)
n_colors <- length(unique(GWAS_selected_hits_matrix_ordered$band))
color_vector <- viridis(n_colors)
colors <- setNames(color_vector, unique(GWAS_selected_hits_matrix_ordered$band))
GWAS_selected_hits_matrix_ordered$color <- mapvalues(
  GWAS_selected_hits_matrix_ordered$band,
  from = names(colors),
  to = colors
)

GWAS_selected_hits_matrix_ordered$selected_genes <- NA
GWAS_selected_hits_matrix_ordered[rownames(GWAS_selected_hits_matrix_ordered) %in% selected_hits,"selected_genes"] <-  rownames(GWAS_selected_hits_matrix_ordered[rownames(GWAS_selected_hits_matrix_ordered) %in% selected_hits,])



# scale?
df <- GWAS_selected_hits_matrix_ordered
# Remove extreme outliers

df <- df %>%
  mutate_all(~ifelse(. > 10, NA, .)) %>%
  mutate_all(~ifelse(. < -10, NA, .))

df$selected_genes <- GWAS_selected_hits_matrix_ordered$selected_genes
# Remove Genes with many NAs
df <- df[-which(rowSums(is.na(df)) > 20),]
# remove cancer types
df <- df[,-which(colSums(is.na(df)) > 10)]
cols <- grep("TCGA-|pancancer",colnames(df))



# df[rownames(df) %in% "SMG5",]
# df2 <- GWAS_final_df[GWAS_final_df$Gene_symbol %in% "SMG5" & GWAS_final_df$dataset == "som_CNV_amp" & GWAS_final_df$NMD_method == "Endogenous",]
# df2[,c("beta_coefficient","cancer_type")]
# Plot
output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/CNV_associations/fine_mapping/selected_genes/chr1q_selected_genes_somatic_assoc_heatmap.png")
png(output_path, width = 5500, height = 3500, res = 300)
p <- heatmap(as.matrix(df[,cols]),  
            #Rowv = NA, 
            na.rm =TRUE,
            #Colv = NA,
            scale = c("col"),
            labRow = GWAS_selected_hits_matrix_ordered$selected_genes,
            RowSideColors = df$color,
            cexRow=1,
            cexCol=1.5)
print(p)
dev.off()



cols <- brewer.pal(n = 5, name = "RdBu")
p <- GWAS_selected_hits %>%
            ggplot(aes(x = Gene_symbol, y = cancer_type, fill = beta_coefficient)) +
            geom_tile(size = 2) +
            geom_text(aes(label = round(-beta_coefficient,2)), color = "black", size = 10) +
            theme_bw(base_size = 35) +
            facet_nested(. ~ NMD_method) +
            theme(axis.text.x = element_text(size = 43, angle = 45),
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
output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/CNV_associations/fine_mapping/selected_genes/chr1q_selected_genes_somatic_assoc.png")
png(output_path, width = 15000, height = 8000, res = 300)
print(p)
dev.off()
