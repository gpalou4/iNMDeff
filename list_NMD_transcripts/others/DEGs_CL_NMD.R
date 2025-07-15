rm(list=ls())

################################ FUNCTIONS ################################
################################ FUNCTIONS ################################
################################ FUNCTIONS ################################

################################ LIBRARIES ################################
################################ LIBRARIES ################################
################################ LIBRARIES ################################

#conda activate /home/gpalou/anaconda3_envs/general

library("edgeR")
library("xlsx")

################################ SCRIPT ################################
################################ SCRIPT ################################
################################ SCRIPT ################################

library("readxl")

# 1) Cell Line data

# 1.3.3) Cell Line UPF1 KD
CL_UPF1KD_path <- "/g/strcombio/fsupek_cancer1/gpalou/cell_lines_UPF1_KD/ENSEMBL/rsem.merged.transcript_counts.tsv"
CL_metadata_path <- "/g/strcombio/fsupek_cancer1/gpalou/cell_lines_UPF1_KD/CL_metadata_subset_and_controls.xlsx"
CL_UPF1KD <- read.table(file = CL_UPF1KD_path, header = TRUE, sep = "\t", row.names = 1)
CL_UPF1KD <- CL_UPF1KD[,-1]
CL_UPF1KD <- CL_UPF1KD[grep(".*PAR.*",rownames(CL_UPF1KD), invert = TRUE),]
rownames(CL_UPF1KD) <- sub("\\..*","", rownames(CL_UPF1KD))
CL_metadata <- data.frame(read_excel(CL_metadata_path))
# Same order as CL data
CL_metadata$ID <- factor(paste0(CL_metadata$GSE_project,"_R",CL_metadata$replicate))
CL_metadata$Type <- as.factor(CL_metadata$Type)
CL_metadata <- CL_metadata[order(CL_metadata$ID),]
CL_UPF1KD <- CL_UPF1KD[,order(colnames(CL_UPF1KD))]

# Change names
# CL 1: HeLa # DISCARD!! #
# cols <- grep("GSE152435_C",colnames(CL_UPF1KD))
# colnames(CL_UPF1KD)[cols] <- gsub("_C_","_WT1_",colnames(CL_UPF1KD)[cols])
# cols <- grep("GSE152435_R.*",colnames(CL_UPF1KD))
# colnames(CL_UPF1KD)[cols] <- gsub("_R","_KD1_R",colnames(CL_UPF1KD)[cols])
# # CL 2: HeLa (Colombo)
# cols <- grep("GSE86148_C.*",colnames(CL_UPF1KD))
# colnames(CL_UPF1KD)[cols] <- gsub("_C_","_WT2",colnames(CL_UPF1KD)[cols])
# cols <- grep("GSE86148_R.*",colnames(CL_UPF1KD))
# colnames(CL_UPF1KD)[cols] <- gsub("_R","_KD2_R",colnames(CL_UPF1KD)[cols])
# # CL 3: HepG2 (ENCODE)
# cols <- grep("GSE88148",colnames(CL_UPF1KD))
# colnames(CL_UPF1KD)[cols] <- gsub("_C_","_WT3",colnames(CL_UPF1KD)[cols])
# cols <- grep("GSE88466",colnames(CL_UPF1KD))
# colnames(CL_UPF1KD)[cols] <- gsub("_R","_KD3_R",colnames(CL_UPF1KD)[cols])
# # CL 4: K562 (ENCODE)
# cols <- grep("GSE88266",colnames(CL_UPF1KD))
# colnames(CL_UPF1KD)[cols] <- gsub("_C_","_WT4",colnames(CL_UPF1KD)[cols])
# cols <- grep("GSE88140",colnames(CL_UPF1KD))
# colnames(CL_UPF1KD)[cols] <- gsub("_R","_KD4_R",colnames(CL_UPF1KD)[cols])

# For each Cell Line separately

DGE_in_CL <- function(CL_UPF1KD, CL_metadata, cell_line) {

    # 1) Info for the cell line
    CL_metadata_filt <- CL_metadata[CL_metadata$CL %in% cell_line,]
    # Remove bad cell line
    if (cell_line == "HeLa") {
        CL_metadata_filt <- CL_metadata_filt[-grep("GSE152435",CL_metadata_filt$GSE_project),]
    }
    CL_UPF1KD_filt <- CL_UPF1KD[,colnames(CL_UPF1KD) %in% CL_metadata_filt$ID]
    # 2) Preprocessing
    # Create DGEList object
    CL_all_DGE <- DGEList(CL_UPF1KD_filt, sample = CL_metadata_filt, group = CL_metadata_filt$Type)
    # Calculate normalization factors
    CL_all_DGE <- calcNormFactors(CL_all_DGE)
    # Filter low-expressed genes
    cutoff <- 1
    drop <- which(apply(cpm(CL_all_DGE), 1, max) < cutoff)
    CL_all_DGE <- CL_all_DGE[-drop,] 
    dim(CL_all_DGE) # number of genes left
    # MDS plot
    MDS_output_path <- paste0("/g/strcombio/fsupek_cancer1/gpalou/cell_lines_UPF1_KD/DGE/MDS_CL_",cell_line,".png")
    png(MDS_output_path, width = 2500, height = 2000, res = 300)
    plotMDS(CL_all_DGE, col = as.numeric(CL_metadata_filt$Type))
    dev.off()
    # 3) Voom transformation and calculation of variance weights
    group <- CL_metadata_filt$Type
    mm <- model.matrix(~0 + group)
    voom_output_path <- paste0("/g/strcombio/fsupek_cancer1/gpalou/cell_lines_UPF1_KD/DGE/Voom_CL_",cell_line,".png")
    png(voom_output_path, width = 2500, height = 2000, res = 300)
    y <- voom(CL_all_DGE, mm, plot = T)
    print(y)
    dev.off()
    # 4) Fitting linear models in limma
    fit <- lmFit(y, mm)
    head(coef(fit))
    contr <- makeContrasts(groupcontrol - groupKD, levels = colnames(coef(fit)))
    contr
    tmp <- contrasts.fit(fit, contr)
    tmp <- eBayes(tmp)
    top_table <- topTable(tmp, sort.by = "P", n = Inf)
    head(top_table, 20)
    # Save and return
    output_path <- paste0("/g/strcombio/fsupek_cancer1/gpalou/cell_lines_UPF1_KD/DGE/top_hits_",cell_line,".txt")
    write.table(top_table, file = output_path, quote=FALSE, sep='\t',row.names = FALSE, col.names = TRUE)
    return(top_table)
}

CL_top_hits_HeLa <- DGE_in_CL(CL_UPF1KD, CL_metadata, cell_line = "HeLa")
CL_top_hits_HepG2 <- DGE_in_CL(CL_UPF1KD, CL_metadata, cell_line = "HepG2")
CL_top_hits_K562 <- DGE_in_CL(CL_UPF1KD, CL_metadata, cell_line = "K562")

# Overlap between Cell Lines at 10%
CL_top_hits_HeLa_transcripts <- rownames(CL_top_hits_HeLa[CL_top_hits_HeLa$adj.P.Val < 0.25,])
CL_top_hits_HepG2_transcripts <- rownames(CL_top_hits_HepG2[CL_top_hits_HepG2$adj.P.Val < 0.25,])
CL_top_hits_K562_transcripts <- rownames(CL_top_hits_K562[CL_top_hits_K562$adj.P.Val < 0.25,])

length(CL_top_hits_HeLa_transcripts)
length(CL_top_hits_HepG2_transcripts)
length(CL_top_hits_K562_transcripts)

final_transcripts <- intersect(CL_top_hits_HeLa_transcripts,CL_top_hits_HepG2_transcripts)
final_transcripts <- intersect(final_transcripts,CL_top_hits_K562_transcripts)

output_path <- paste0("/g/strcombio/fsupek_cancer1/gpalou/cell_lines_UPF1_KD/DGE/top_hits_overlap_cell_lines.txt")
write.table(final_transcripts, file = output_path, quote=FALSE, sep='\t',row.names = FALSE, col.names = TRUE)
