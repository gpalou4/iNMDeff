rm(list=ls())

################################################################################################

########################################## FUNCTIONS ###########################################

################################################################################################

transcripts_WT_expression <- function(PTCs.transcripts.df, RNAseq.TCGA, transcripts.germ.mut, TCGA.barcode, transcripts.to.remove, RNAseq_type) {
    PTC.transcripts <- PTCs.transcripts.df$transcript_id
    PTCs.transcripts.df[,paste0(RNAseq_type,"_WT_samples")] <- NA
    PTCs.transcripts.df[,paste0(RNAseq_type,"_transcript_expression_WT")] <- NA
    PTCs.transcripts.df[,paste0(RNAseq_type,"_median_sample_WT")] <- NA
    # 1) Obtain the median gene expression of the wild-type transcripts from the rest of the samples
    # Obtain germline mutations in those transcripts from the rest of the samples
    transcripts.germ.mut.filt <- transcripts.germ.mut[rownames(transcripts.germ.mut)%in%PTC.transcripts,!colnames(transcripts.germ.mut)%in%TCGA.barcode,drop = FALSE]
    # Obtain somatic CNV/mut in those transcripts from the rest of the samples
    transcripts.to.remove.filt <- transcripts.to.remove[rownames(transcripts.to.remove)%in%PTC.transcripts,!colnames(transcripts.to.remove)%in%TCGA.barcode,drop = FALSE]
    # For each PTC transcript obtain the sample names for those ones that contain wild-type transcripts
    # Obtain samples with same subtype
    sample.cluster <- cancer.subtypes[rownames(cancer.subtypes)%in%TCGA.barcode,"cluster"]
    samples.control.cluster <- cancer.subtypes[cancer.subtypes$cluster == sample.cluster,]
    sapply(time(rownames(transcripts.germ.mut.filt)), function(row.number) {
        PTC.transcript <- rownames(transcripts.germ.mut.filt)[row.number]
        PTC.row <- PTCs.transcripts.df$transcript_id %in% PTC.transcript
        # Remove samples with germline mut in the transcript
        mutations <- "stopgain|[^non]frameshift deletion|[^non]frameshift insertion|splicing|nonframeshift deletion|nonframeshift insertion|startloss|stoploss"
        samples.control.transcript <- colnames(transcripts.germ.mut.filt[rownames(transcripts.germ.mut.filt)==PTC.transcript,grep(mutations,transcripts.germ.mut.filt[rownames(transcripts.germ.mut.filt)==PTC.transcript,], invert = TRUE),])
        # Remove samples with overlapping somatic CNV/mut in the transcript        
        samples.control.transcript.somatic <- colnames(transcripts.to.remove.filt[,which(transcripts.to.remove.filt[rownames(transcripts.to.remove.filt)==PTC.transcript,] != 0)])
        samples.control.transcript <- samples.control.transcript[!samples.control.transcript %in% samples.control.transcript.somatic]
        # Keep only samples from the same subtype
        samples.control.transcript <- samples.control.transcript[samples.control.transcript%in%rownames(samples.control.cluster)]
        # Check that we have enough samples
        if (length(samples.control.transcript) == 0) {
            PTCs.transcripts.df[PTC.row,paste0(RNAseq_type,"_WT_samples")] <<- 0
        } else {
            # Obtain the median gene expression for the wild-type transcript
            PTC.transcript.wt.gene.exp <- RNAseq.TCGA[rownames(RNAseq.TCGA)==PTC.transcript,colnames(RNAseq.TCGA)%in%samples.control.transcript, drop = FALSE]
            # 1) If < 10 wild-type measurements (samples) per transcript per NMF cluster available we should not use the PTC
            num.WT.samples <- ncol(PTC.transcript.wt.gene.exp)
            PTCs.transcripts.df[PTC.row,paste0(RNAseq_type,"_WT_samples")] <<- num.WT.samples
            # 2) Median expression
            # Obtain TCGA barcode from the median sample
            if (num.WT.samples != 0 ){
                mid.value <- num.WT.samples/2
                if(mid.value == 0.5) {mid.value <- 1}
                # sort
                PTC.transcript.wt.gene.exp <- PTC.transcript.wt.gene.exp[,order(as.numeric(PTC.transcript.wt.gene.exp[1,])), drop = FALSE]
                tmp <- PTC.transcript.wt.gene.exp[,ceiling(mid.value), drop = FALSE]
                median.sample.name <- colnames(tmp)
                transcripts.wt.median <- as.numeric(tmp)
                PTCs.transcripts.df[PTC.row,paste0(RNAseq_type,"_transcript_expression_WT")] <<- transcripts.wt.median
                PTCs.transcripts.df[PTC.row,paste0(RNAseq_type,"_median_sample_WT")] <<- median.sample.name
            } else {
                PTCs.transcripts.df[PTC.row,paste0(RNAseq_type,"_transcript_expression_WT")] <<- NA
                PTCs.transcripts.df[PTC.row,paste0(RNAseq_type,"_median_sample_WT")] <<- NA
            }
        }
    })
    if (is.vector(PTCs.transcripts.df)) {
        return(NULL)
    } else {
        return(PTCs.transcripts.df)
    }
}

################################################################################################

########################################## LIBRARIES ###########################################

################################################################################################
.libPaths( rev( .libPaths() ) )
library("stringr")

##############################################################################################

########################################## SCRIPT ##############################################

################################################################################################

# 1) Load Data
# Arguments and paths

args <- commandArgs(trailingOnly=TRUE)
paths <- read.table(file = args[1], header = TRUE, sep = ",", stringsAsFactors = FALSE)
TCGA.cancer <- args[2]
TCGA.barcode <- args[3]
VCF.type <- args[4]
variant.type <- args[5]
TCGA.barcode <- gsub("-",".",TCGA.barcode)

#paths <- read.table(file = "/home/gpalou/projects/NMD/scripts/NMD_efficiency/PTC_NMD_rules/TCGA/NMD_rules_and_efficiency.txt", header = TRUE, sep = ",", stringsAsFactors = FALSE)

RNAseq.path <- paths[paths$folder_or_object=="RNAseq_path","path_or_filename"]
conversor.tables.path <- paths[paths$folder_or_object=="conversor_tables","path_or_filename"]
transcripts.to.remove.path <- paths[paths$folder_or_object=="transcripts_to_remove_path","path_or_filename"]
transcripts.germ.mut.path <- paths[paths$folder_or_object=="transcripts_germ_mut_path","path_or_filename"]
NMD.targets.path <- paths[paths$folder_or_object=="NMD_targets_path","path_or_filename"]
firehose.subtypes.path <- paths[paths$folder_or_object=="firehose_subtypes_path","path_or_filename"]
TCGA.RNAseq.PCA.path <- paths[paths$folder_or_object=="TCGA_RNAseq_PCA_path","path_or_filename"]
TCGA.metadata.path <- paths[paths$folder_or_object=="TCGA_metadata_path","path_or_filename"]
PTC.transcripts.output.path <- paths[paths$folder_or_object=="PTC_transcripts_output_path","path_or_filename"]

print(paste0("TCGA cancer --> ",TCGA.cancer))
print(paste0("TCGA sample  --> ",TCGA.barcode))

# 1.1) PTCs transcripts table
# Check if file exists (otherwise skip the script)
PTCs.df.path <- gsub("\\[X1\\]",TCGA.cancer, PTC.transcripts.output.path)
PTCs.df.path <- gsub("\\[X2\\]",gsub("\\.","-",TCGA.barcode), PTCs.df.path)
PTCs.df.path <- paste0(PTCs.df.path,paths[paths$folder_or_object==paste0(VCF.type,"_PTC_transcripts_output"),"path_or_filename"])
if (variant.type == "missense") {
  PTCs.df.path <- gsub(paste0(VCF.type,"_PTC_"),paste0(VCF.type,"_missense_syn_"),PTCs.df.path)
} 


if (!file.exists(PTCs.df.path)) {
    print("PTCs transcripts file does not exist")
} else {
        PTCs.df <- read.table(file = PTCs.df.path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
        print(PTCs.df[1:10,1:10])
        # 1.2.1) RNAseq TPM
        RNAseq.TCGA.TPM <- read.table(file = gsub("\\[X\\]",TCGA.cancer, paste0(RNAseq.path,paths[paths$folder_or_object=="RNAseq_TPM","path_or_filename"])), 
                                    header = TRUE, sep = "\t", row.names = 1)
        # 1.2.2) RNAseq raw
        RNAseq.TCGA.raw <- read.table(file = gsub("\\[X\\]",TCGA.cancer, paste0(RNAseq.path,paths[paths$folder_or_object=="RNAseq_raw","path_or_filename"])), 
                                    header = TRUE, sep = "\t", row.names = 1)
        if (TCGA.barcode %in% colnames(RNAseq.TCGA.TPM)) {
            # 1.3) Transcripts length
            ensembl.v88.transcripts.length <- read.table(file = paste0(conversor.tables.path,paths[paths$folder_or_object=="ensembl_v88_transcripts_length","path_or_filename"]), 
                                            header = TRUE, sep = "\t")
            # 1.4) Protein coding transcripts
            ensembl.v88.coding.transcripts <- read.table(file = paste0(conversor.tables.path,paths[paths$folder_or_object=="ensembl_v88_coding_transcripts","path_or_filename"]),
                                            header = FALSE, sep = "\t")
            ensembl.v88.coding.transcripts <- ensembl.v88.coding.transcripts[,1, drop = FALSE]
            ensembl.v88.coding.transcripts$protein_coding <- "yes"
            # 1.5) transcripts filtering matrix with somatic AND CNV
            transcripts.to.remove <- read.table(file = gsub("\\[X\\]",gsub("TCGA-","",TCGA.cancer), paste0(transcripts.to.remove.path,paths[paths$folder_or_object=="transcripts_to_remove","path_or_filename"])),
                                                header = TRUE, sep = "\t", row.names = 1)
            # 1.6) transcripts with germline mutations
            transcripts.germ.mut <- read.table(file = gsub("\\[X\\]",gsub("TCGA-","",TCGA.cancer), paste0(transcripts.germ.mut.path,paths[paths$folder_or_object=="transcripts_germ_mut","path_or_filename"])),
                                                header = TRUE, sep = "\t", row.names = 1)
            transcripts.germ.mut[] <- lapply(transcripts.germ.mut , as.character)
            # 1.7) ENSEMBL transcripts IDs hg38 GTF
            ensembl.v88.gtf <- rtracklayer::import(paste0(conversor.tables.path,paths[paths$folder_or_object=="ensembl_v88_gtf","path_or_filename"]))
            ensembl.v88.gtf <- as.data.frame(ensembl.v88.gtf)
            # Remove versions from IDs
            ensembl.v88.gtf[,c("gene_id","transcript_id")] <- sapply(ensembl.v88.gtf[,c("gene_id","transcript_id")], function(col) { 
                gsub("(.*)\\..*","\\1", col)
            })
            # 1.8) LOEUF score (negative selected genes)
            LOEUF.score.table <- read.table(file = paste0(conversor.tables.path,paths[paths$folder_or_object=="LOEUF_score","path_or_filename"]), header = TRUE, sep = "\t")
            LOEUF.score.ensembl <- LOEUF.score.table[,c("transcript","oe_lof_upper","oe_lof_upper_bin")]
            # 1.9) Positive selected genes (PSG)
            PSG.ensembl.gene <- read.table(file = paste0(conversor.tables.path,paths[paths$folder_or_object=="PSG_331","path_or_filename"]),
                                            header = TRUE, sep = "\t")
            # 1.10) NMD targets (for negative control)
            NMD.global.ensembl <- read.table(file = paste0(NMD.targets.path,paths[paths$folder_or_object=="NMD_global_ensembl","path_or_filename"]),
                                            header = TRUE, sep = "\t", stringsAsFactors = FALSE)
            NMD.global.ensembl.transcripts <- NMD.global.ensembl$ensembl_transcript_id
            # 1.11) Firehose subtypes (mRNA NMF) and CNV
            if (TCGA.cancer == "TCGA-SKCM") {
                cancer.subtypes.path <- gsub("-TP","-TM",paste0(firehose.subtypes.path,paths[paths$folder_or_object=="firehose_subtypes","path_or_filename"]))
            }  else if (TCGA.cancer == "TCGA-LAML")  {
                cancer.subtypes.path <- gsub("-TP","-TB",paste0(firehose.subtypes.path,paths[paths$folder_or_object=="firehose_subtypes","path_or_filename"]))
            }   else {
                cancer.subtypes.path <- paste0(firehose.subtypes.path,paths[paths$folder_or_object=="firehose_subtypes","path_or_filename"])
            }
            # Subtypes
            cancer.subtypes <- read.table(file = gsub("\\[X\\]",gsub("TCGA-","",TCGA.cancer), cancer.subtypes.path), 
                                                    header = TRUE, sep = "\t", skip = 1)
            cancer.subtypes$SampleName<- gsub("-",".", cancer.subtypes$SampleName, fixed = TRUE)
            rownames(cancer.subtypes) <- substr(cancer.subtypes$SampleName,1,12)
            # 1.12) Median transcript expression and coefficient of variation
            transcripts.median.exp <- read.table(file = gsub("\\[X\\]",TCGA.cancer, paste0(TCGA.metadata.path,paths[paths$folder_or_object=="transcripts_exp","path_or_filename"])),
                                                header = TRUE, sep = "\t", row.names = 1)
            
            # 2) Add metadata information for all PTC transcripts
            print("ADDING METADATA INFORMATION FOR ALL PTC TRANSCRIPTS")
            # 2.1) Transcript length (it includes UTRs)
            PTCs.df <- merge(PTCs.df,ensembl.v88.transcripts.length, by.x = "transcript_id", by.y = "transcript_id", all.x = TRUE)
            colnames(PTCs.df)[colnames(PTCs.df) %in% "length"] <- "transcript_length"
            # 2.2) Protein coding?
            # Genes with 1 exon are not in the list for some reason
            PTCs.df <- merge(PTCs.df,ensembl.v88.coding.transcripts, by.x = "transcript_id", by.y = "V1", all.x = TRUE)
            colnames(PTCs.df)[colnames(PTCs.df) %in% "length"] <- "transcript_length"
            # 2.3) Sample-specific somatic CNV/SNVs
            PTCs.df$somatic_CNV_SNV <- "no"
            if (TCGA.barcode%in%colnames(transcripts.to.remove)) {
                sample.transcripts.to.remove <- rownames(transcripts.to.remove[transcripts.to.remove[,TCGA.barcode]==1,])
                PTCs.df[PTCs.df$transcript_id %in% sample.transcripts.to.remove,"somatic_CNV_SNV"] <- "yes"
            }
            # 2.4) Sample-specific germline CNV/SNVs
            PTCs.df$germline_SNV <- "no"
            if (TCGA.barcode%in%colnames(transcripts.germ.mut)) {
                # 2.4) Sample-specific germline SNVs
                transcripts.germ.mut.sample <- transcripts.germ.mut[,TCGA.barcode, drop = FALSE]
                sample.transcripts.to.remove <- c()
                germ.muts <- c("stopgain","[^non]frameshift insertion","[^non]frameshift deletion","splicing")
                for (mutation in germ.muts) {
                    if (mutation == "splicing") {next}
                    # 1) Check which transcripts contains mutations
                    transcripts.mut.sample <- transcripts.germ.mut.sample[grep(mutation, transcripts.germ.mut.sample[,1]),, drop = FALSE]
                    # 2) Check if those transcripts contains additional - the same mutation (e.g. [exon9, exon8]_stopgain)
                    index <- grep(paste0("\\[exon[0-9]{1,4}\\]_",gsub("\\[\\^non\\]","",mutation),".*"),transcripts.mut.sample[,1], invert = TRUE)
                    additional.sameMut.transcripts <- rownames(transcripts.mut.sample[index,, drop = FALSE])
                    # 3) Check if those transcripts contains additional - two different mutations
                    germ.mut.tmp <- germ.muts[!germ.muts %in% mutation]
                    germ.mut.tmp <- paste0(germ.mut.tmp,collapse="|")
                    additional.diffMut.transcripts <- rownames(transcripts.mut.sample[grep(germ.mut.tmp, transcripts.mut.sample[,1]),, drop = FALSE])
                    # 4) Sum both sets
                    sample.transcripts.to.remove <- unique(c(sample.transcripts.to.remove,unique(c(additional.sameMut.transcripts,additional.diffMut.transcripts))))
                }
                PTCs.df[PTCs.df$transcript_id %in% sample.transcripts.to.remove,"germline_SNV"] <- "yes"
            }
            # 2.5) LOEUF score
            PTCs.df <- merge(PTCs.df,LOEUF.score.ensembl, by.x = "transcript_id", by.y = "transcript", all.x = TRUE)
            colnames(PTCs.df)[colnames(PTCs.df) %in% c("oe_lof_upper","oe_lof_upper_bin")] <- c("LOEUF_score","LOEUF_decile")
            # 2.6) PSG transcript
            PTCs.df$PSG <- "no"
            PTCs.df[PTCs.df$transcript_id %in% PSG.ensembl.gene$Ensembl.Transcript.ID,"PSG"] <- "yes"
            # 2.7) NMD target
            PTCs.df$NMD_target <- "no"
            PTCs.df[PTCs.df$transcript_id %in% NMD.global.ensembl.transcripts,"NMD_target"] <- "yes"
            # 2.8) Transcript expression
            # 2.8.1) TPM
            RNAseq.TCGA.TPM.sample <- RNAseq.TCGA.TPM[,TCGA.barcode, drop = FALSE]
            PTCs.df <- merge(PTCs.df,RNAseq.TCGA.TPM.sample, by.x = "transcript_id", by.y = "row.names", all.x = TRUE)
            colnames(PTCs.df)[colnames(PTCs.df) %in% TCGA.barcode] <- "TPM_transcript_expression"
            # 2.8.2) Raw
            RNAseq.TCGA.raw.sample <- RNAseq.TCGA.raw[,TCGA.barcode, drop = FALSE]
            PTCs.df <- merge(PTCs.df,RNAseq.TCGA.raw.sample, by.x = "transcript_id", by.y = "row.names", all.x = TRUE)
            colnames(PTCs.df)[colnames(PTCs.df) %in% TCGA.barcode] <- "Raw_transcript_expression"
            # 2.9) Median transcript expression and coefficient of variation (SD/mean)
            # Transcripts with large variation in expression among the tumor samples (coefficient of variation > 0.5)
            PTCs.df <- merge(PTCs.df, transcripts.median.exp, by.x = "transcript_id", by.y = "row.names", all.x = TRUE)
            # 2.10) Low-expressed transcript in the cancer?
            # At least 3 TPM or 20 raw counts
            RNAseq.TCGA.TPM.filt <- RNAseq.TCGA.TPM[rowSums(log2(RNAseq.TCGA.TPM) >= 1.58, na.rm = TRUE) >= round(length(colnames(RNAseq.TCGA.TPM)) * 0.50),]
            RNAseq.TCGA.raw.filt <- RNAseq.TCGA.raw[rowSums(log2(RNAseq.TCGA.raw) >= 4.32, na.rm = TRUE) >= round(length(colnames(RNAseq.TCGA.raw)) * 0.50),]
            PTCs.df$TPM_low_expressed_transcript <- "yes"
            PTCs.df[PTCs.df$transcript_id %in% rownames(RNAseq.TCGA.TPM.filt),"TPM_low_expressed_transcript"] <- "no"
            PTCs.df$Raw_low_expressed_transcript <- "yes"
            PTCs.df[PTCs.df$transcript_id %in% rownames(RNAseq.TCGA.raw.filt),"Raw_low_expressed_transcript"] <- "no"
            # 2.11) PCs for transcripts
            # 2.11.1) Tissue
            PCA.path <- paste0(TCGA.RNAseq.PCA.path,TCGA.cancer,"_tissue_PCA_var.txt")
            PCA.df <- read.table(file = PCA.path, header = TRUE, sep = "\t", row.names = 1)
            colnames(PCA.df) <- c("tissue_PC1","tissue_PC2","tissue_PC3","tissue_PC4")
            PTCs.df <- merge(PTCs.df, PCA.df, by.x = "transcript_id", by.y = "row.names", all.x = TRUE)
            # 2.11.2) Subtissue
            for (i in 1:max(cancer.subtypes$cluster)) {
                PCA.path <- paste0(TCGA.RNAseq.PCA.path,TCGA.cancer,"_subtissue_",i,"_PCA_var.txt")
                error <- FALSE
                tryCatch( { 
                    PCA.df <- read.table(file = PCA.path, header = TRUE, sep = "\t", row.names = 1)
                    colnames(PCA.df) <- paste("subtissue_",i,"_PC",1:4, sep = "")
                    PTCs.df <- merge(PTCs.df, PCA.df, by.x = "transcript_id", by.y = "row.names", all.x = TRUE)
                },error = function(e) {error <<- TRUE})
            }
            # 2.12) transcript expression of WT samples from same RNAseq NMF subtype
            # Also add library size from the WT median sample
            # 2.12.1) TPM
            PTCs.df$TPM_transcript_expression_WT <- NA
            df.tmp <- transcripts_WT_expression(PTCs.transcripts.df = PTCs.df[,c("transcript_id","start_pos","TPM_transcript_expression_WT")], TCGA.barcode = TCGA.barcode,
                                        RNAseq.TCGA = RNAseq.TCGA.TPM, transcripts.germ.mut = transcripts.germ.mut, transcripts.to.remove = transcripts.to.remove, RNAseq_type = "TPM")
            PTCs.df$TPM_transcript_expression_WT <- NULL
            PTCs.df <- merge(PTCs.df,df.tmp, by.x = c("start_pos","transcript_id"), by.y = c("start_pos","transcript_id"), all.x = TRUE)
            # 2.12.2) Raw
            PTCs.df$Raw_transcript_expression_WT <- NA
            df.tmp <- transcripts_WT_expression(PTCs.transcripts.df = PTCs.df[,c("transcript_id","start_pos","Raw_transcript_expression_WT")], TCGA.barcode = TCGA.barcode,
                                        RNAseq.TCGA = RNAseq.TCGA.raw, transcripts.germ.mut = transcripts.germ.mut, transcripts.to.remove = transcripts.to.remove, RNAseq_type = "Raw")
            PTCs.df$Raw_transcript_expression_WT <- NULL
            PTCs.df <- merge(PTCs.df,df.tmp, by.x = c("start_pos","transcript_id"), by.y = c("start_pos","transcript_id"), all.x = TRUE)
            # Remove duplicated
            PTCs.df <- PTCs.df[!duplicated(PTCs.df),]
            # 2.13) NMD efficiency
            PTCs.df$NMD_efficiency_TPM <- -log2(PTCs.df$TPM_transcript_expression/PTCs.df$TPM_transcript_expression_WT)
            PTCs.df$NMD_efficiency_Raw <- -log2(PTCs.df$Raw_transcript_expression/PTCs.df$Raw_transcript_expression_WT)
            # 3) Save results
            PTCs.df.metadata.path <- gsub("transcripts.txt","transcripts_metadata.txt",PTCs.df.path)
            write.table(PTCs.df, file = PTCs.df.metadata.path,
                        sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
    } else {
        print("Sample not found in RNAseq matrix, skipping...")
    }
}


