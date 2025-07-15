rm(list=ls())

################################################################################################

########################################## FUNCTIONS ###########################################

################################################################################################

transcripts_WT_expression <- function(PTCs.transcripts.df, RNAseq.GTEx, transcripts.germ.mut, GTEx.sample, RNAseq_type) {

        PTC.transcripts <- PTCs.transcripts.df$transcript_id
        PTCs.transcripts.df[,paste0(RNAseq_type,"_WT_samples")] <- NA
        PTCs.transcripts.df[,paste0(RNAseq_type,"_transcript_expression_WT")] <- NA
        PTCs.transcripts.df[,paste0(RNAseq_type,"_median_sample_WT")] <- NA
        # 1) Obtain the median gene expression of the wild-type transcripts from the rest of the samples
        # Obtain germline mutations in those transcripts from the rest of the samples
        transcripts.germ.mut.filt <- transcripts.germ.mut[rownames(transcripts.germ.mut)%in%PTC.transcripts,!colnames(transcripts.germ.mut)%in%GTEx.sample,drop = FALSE]
        colnames(transcripts.germ.mut.filt) <- gsub("(GTEX\\-\\w{4,5})\\-.*","\\1",colnames(transcripts.germ.mut.filt))        
        samples.control <- colnames(transcripts.germ.mut.filt)
        # For each PTC transcript obtain the sample names for those ones that contain wild-type transcripts
        sapply(time(rownames(transcripts.germ.mut.filt)), function(row.number) {
            PTC.transcript <- rownames(transcripts.germ.mut.filt)[row.number]
            PTC.row <- PTCs.transcripts.df$transcript_id %in% PTC.transcript
            # Remove samples with germline mut in the transcript
            mutations <- "stopgain|[^non]frameshift deletion|[^non]frameshift insertion|splicing|nonframeshift deletion|nonframeshift insertion|startloss|stoploss"
            samples.control.transcript <- colnames(transcripts.germ.mut.filt[rownames(transcripts.germ.mut.filt)%in%PTC.transcript,grep(mutations,transcripts.germ.mut.filt[rownames(transcripts.germ.mut.filt)==PTC.transcript,], invert = TRUE),])
            # Check that we have enough samples
            if (length(samples.control.transcript) == 0) {
                PTCs.transcripts.df[PTC.row,paste0(RNAseq_type,"_WT_samples")] <<- 0
            } else {
                # Obtain the median gene expression for the wild-type transcript
                PTC.transcript.wt.gene.exp <- RNAseq.GTEx[rownames(RNAseq.GTEx)==PTC.transcript,colnames(RNAseq.GTEx)%in%samples.control.transcript, drop = FALSE]
                # 2) If < 10 wild-type measurements (samples) per transcript per NMF cluster available we should not use the PTC
                num.WT.samples <- ncol(PTC.transcript.wt.gene.exp)
                PTCs.transcripts.df[PTC.row,paste0(RNAseq_type,"_WT_samples")] <<- num.WT.samples
                # 3) Median expression
                # Obtain GTEx barcode from the median sample
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

library("stringr")

##############################################################################################

########################################## SCRIPT ##############################################

################################################################################################

# 1) Load Data
# Arguments and paths

args <- commandArgs(trailingOnly=TRUE)
paths <- read.table(file = args[1], header = TRUE, sep = ",", stringsAsFactors = FALSE)
GTEx.tissue <- args[2]
GTEx.sample <- args[3]
VCF.type <- args[4]

#paths <- read.table(file = "/home/gpalou/projects/NMD/scripts/NMD_efficiency/PTC_NMD_rules/GTEx/NMD_rules_and_efficiency.txt", header = TRUE, sep = ",", stringsAsFactors = FALSE)

RNAseq.path <- paths[paths$folder_or_object=="RNAseq_path","path_or_filename"]
conversor.tables.path <- paths[paths$folder_or_object=="conversor_tables","path_or_filename"]
transcripts.germ.mut.tissue.path <- paths[paths$folder_or_object=="transcripts_germ_mut_tissue_path","path_or_filename"]
NMD.targets.path <- paths[paths$folder_or_object=="NMD_targets_path","path_or_filename"]
PTC.transcripts.output.path <- paths[paths$folder_or_object=="PTC_transcripts_output_path","path_or_filename"]
PTC.transcripts.metadata.output.path <- paths[paths$folder_or_object=="PTC_transcripts_metadata_output_path","path_or_filename"]
GTEx.metadata.path <- paths[paths$folder_or_object=="GTEx_metadata_path","path_or_filename"]
GTEx.RNAseq.PCA.path <- paths[paths$folder_or_object=="GTEx_RNAseq_PCA_path","path_or_filename"]

print(paste0("GTEx tissue --> ",GTEx.tissue))
print(paste0("GTEx sample  --> ",GTEx.sample))

# 1.1) PTCs transcripts table
# Check if file exists (otherwise skip the script)
GTEx.sample.filt <- gsub("(GTEX\\-\\w{4,5})\\-.*","\\1",GTEx.sample)
PTCs.df.path <- gsub("\\[X1\\]",GTEx.sample.filt, paste0(PTC.transcripts.output.path,paths[paths$folder_or_object=="germline_PTC_transcripts_output","path_or_filename"]))
PTCs.df.metadata.path <- gsub("\\[X1\\]",GTEx.tissue, paste0(PTC.transcripts.metadata.output.path,paths[paths$folder_or_object=="germline_PTC_transcripts_metadata_output","path_or_filename"]))
PTCs.df.metadata.path <- gsub("\\[X2\\]",GTEx.sample.filt, PTCs.df.metadata.path)

if (!file.exists(PTCs.df.path)) {
    print("PTCs transcripts file does not exist")
} else {
        PTCs.df <- read.table(file = PTCs.df.path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
        print(PTCs.df[1:10,1:10])
        # 1.2.1) RNAseq TPM
        RNAseq.GTEx.TPM <- read.table(file = gsub("\\[X\\]",GTEx.tissue, paste0(RNAseq.path,paths[paths$folder_or_object=="RNAseq_TPM","path_or_filename"])), 
                                    header = TRUE, row.names = 1)
        colnames(RNAseq.GTEx.TPM) <- gsub("\\.","-",gsub("(GTEX\\.\\w{4,5})\\..*","\\1",colnames(RNAseq.GTEx.TPM)))
        RNAseq.GTEx.TPM$gene_id <- NULL
        # Remove X-Y genes
        RNAseq.GTEx.TPM <- RNAseq.GTEx.TPM[-grep("PAR",rownames(RNAseq.GTEx.TPM)),]
        rownames(RNAseq.GTEx.TPM) <- gsub("(ENST.*)\\..*","\\1",rownames(RNAseq.GTEx.TPM))
        # Round values
        # 1.2.2) RNAseq raw
        RNAseq.GTEx.raw <- read.table(file = gsub("\\[X\\]",GTEx.tissue, paste0(RNAseq.path,paths[paths$folder_or_object=="RNAseq_raw","path_or_filename"])), 
                                    header = TRUE, row.names = 1)
        colnames(RNAseq.GTEx.raw) <- gsub("\\.","-",gsub("(GTEX\\.\\w{4,5})\\..*","\\1",colnames(RNAseq.GTEx.raw)))
        RNAseq.GTEx.raw$gene_id <- NULL
        # Remove X-Y genes
        RNAseq.GTEx.raw <- RNAseq.GTEx.raw[-grep("PAR",rownames(RNAseq.GTEx.raw)),]
        rownames(RNAseq.GTEx.raw) <- gsub("(ENST.*)\\..*","\\1",rownames(RNAseq.GTEx.raw))
        # Round values
        RNAseq.GTEx.raw <- round(RNAseq.GTEx.raw)
        # samples
        GTEx.samples <- colnames(RNAseq.GTEx.raw)
        if (GTEx.sample.filt %in% colnames(RNAseq.GTEx.TPM)) {
            # 1.3) Transcripts length
            ensembl.v88.transcripts.length <- read.table(file = paste0(conversor.tables.path,paths[paths$folder_or_object=="ensembl_v88_transcripts_length","path_or_filename"]), 
                                            header = TRUE, sep = "\t")
            # 1.4) Protein coding transcripts
            ensembl.v88.coding.transcripts <- read.table(file = paste0(conversor.tables.path,paths[paths$folder_or_object=="ensembl_v88_coding_transcripts","path_or_filename"]),
                                            header = FALSE, sep = "\t")
            ensembl.v88.coding.transcripts <- ensembl.v88.coding.transcripts[,1, drop = FALSE]
            ensembl.v88.coding.transcripts$protein_coding <- "yes"
            # 1.5.1) transcripts with germline mutations
            transcripts.germ.mut.tissue.path <-  gsub("\\[X\\]",GTEx.tissue,paste0(transcripts.germ.mut.tissue.path,paths[paths$folder_or_object=="transcripts_germ_mut_tissue","path_or_filename"]))
            if (!file.exists(transcripts.germ.mut.tissue.path)) {
                print("VCF germline file for tissue not found, skipping...")
            } else {
                transcripts.germ.mut.tissue <- read.table(file = transcripts.germ.mut.tissue.path, header = TRUE, sep = "\t", row.names = 1)
                transcripts.germ.mut.tissue[] <- lapply(transcripts.germ.mut.tissue , as.character)
                colnames(transcripts.germ.mut.tissue) <- gsub("\\.","-",colnames(transcripts.germ.mut.tissue))
                colnames(transcripts.germ.mut.tissue) <- gsub("(GTEX\\-\\w{4,5})\\-.*","\\1",colnames(transcripts.germ.mut.tissue))
            }
            # 1.6) ENSEMBL transcripts IDs hg38 GTF
            ensembl.v88.gtf <- rtracklayer::import(paste0(conversor.tables.path,paths[paths$folder_or_object=="ensembl_v88_gtf","path_or_filename"]))
            ensembl.v88.gtf <- as.data.frame(ensembl.v88.gtf)
            # Remove versions from IDs
            ensembl.v88.gtf[,c("gene_id","transcript_id")] <- sapply(ensembl.v88.gtf[,c("gene_id","transcript_id")], function(col) { 
                gsub("(.*)\\..*","\\1", col)
            })
            # 1.7) LOEUF score (negative selected genes)
            LOEUF.score.table <- read.table(file = paste0(conversor.tables.path,paths[paths$folder_or_object=="LOEUF_score","path_or_filename"]), header = TRUE, sep = "\t")
            LOEUF.score.ensembl <- LOEUF.score.table[,c("transcript","oe_lof_upper","oe_lof_upper_bin")]
            # 1.8) Positive selected genes (PSG)
            PSG.ensembl.gene <- read.table(file = paste0(conversor.tables.path,paths[paths$folder_or_object=="PSG_331","path_or_filename"]),
                                            header = TRUE, sep = "\t")
            # 1.9) NMD targets (for negative control)
            NMD.global.ensembl <- read.table(file = paste0(NMD.targets.path,paths[paths$folder_or_object=="NMD_global_ensembl","path_or_filename"]),
                                            header = TRUE, sep = "\t", stringsAsFactors = FALSE)
            NMD.global.ensembl.transcripts <- NMD.global.ensembl$ensembl_transcript_id
            # 1.10) Median transcript expression and coefficient of variation
            GTEx.metadata.path <- gsub("\\[X1\\]",GTEx.tissue, GTEx.metadata.path)
            transcripts.median.exp <- read.table(file = gsub("\\[X2\\]",GTEx.tissue, paste0(GTEx.metadata.path,paths[paths$folder_or_object=="GTEx_transcripts_exp","path_or_filename"])),
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
            # 2.3) Sample-specific germline SNVs
            PTCs.df$germline_SNV <- "no"
            if (GTEx.sample%in%colnames(transcripts.germ.mut.tissue)) {
                # 2.4) Sample-specific germline SNVs
                transcripts.germ.mut.sample <- transcripts.germ.mut.tissue[,GTEx.sample, drop = FALSE]
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
            # 2.4) LOEUF score
            PTCs.df <- merge(PTCs.df,LOEUF.score.ensembl, by.x = "transcript_id", by.y = "transcript", all.x = TRUE)
            colnames(PTCs.df)[colnames(PTCs.df) %in% c("oe_lof_upper","oe_lof_upper_bin")] <- c("LOEUF_score","LOEUF_decile")
            # 2.5) PSG transcript
            PTCs.df$PSG <- "no"
            PTCs.df[PTCs.df$transcript_id %in% PSG.ensembl.gene$Ensembl.Transcript.ID,"PSG"] <- "yes"
            # 2.6) NMD target
            PTCs.df$NMD_target <- "no"
            PTCs.df[PTCs.df$transcript_id %in% NMD.global.ensembl.transcripts,"NMD_target"] <- "yes"
            # 2.7) Transcript expression
            # 2.7.1) TPM
            RNAseq.GTEx.TPM.sample <- RNAseq.GTEx.TPM[,GTEx.sample.filt, drop = FALSE]
            PTCs.df <- merge(PTCs.df,RNAseq.GTEx.TPM.sample, by.x = "transcript_id", by.y = "row.names", all.x = TRUE)
            colnames(PTCs.df)[colnames(PTCs.df) %in% GTEx.sample.filt] <- "TPM_transcript_expression"
            # 2.7.2) Raw
            RNAseq.GTEx.raw.sample <- RNAseq.GTEx.raw[,GTEx.sample.filt, drop = FALSE]
            PTCs.df <- merge(PTCs.df,RNAseq.GTEx.raw.sample, by.x = "transcript_id", by.y = "row.names", all.x = TRUE)
            colnames(PTCs.df)[colnames(PTCs.df) %in% GTEx.sample.filt] <- "Raw_transcript_expression"
            # 2.8) Median transcript expression and coefficient of variation (SD/mean)
            # Transcripts with large variation in expression among the tumor samples (coefficient of variation > 0.5)
            PTCs.df <- merge(PTCs.df, transcripts.median.exp, by.x = "transcript_id", by.y = "row.names", all.x = TRUE)
            # 2.9) Low-expressed transcript in the tissue?
            # At least 3 TPM or 20 raw counts
            RNAseq.GTEx.TPM.filt <- RNAseq.GTEx.TPM[rowSums(log2(RNAseq.GTEx.TPM) >= 1.58, na.rm = TRUE) >= round(length(colnames(RNAseq.GTEx.TPM)) * 0.50),]
            RNAseq.GTEx.raw.filt <- RNAseq.GTEx.raw[rowSums(log2(RNAseq.GTEx.raw) >= 4.32, na.rm = TRUE) >= round(length(colnames(RNAseq.GTEx.raw)) * 0.50),]
            PTCs.df$TPM_low_expressed_transcript <- "yes"
            PTCs.df[PTCs.df$transcript_id %in% rownames(RNAseq.GTEx.TPM.filt),"TPM_low_expressed_transcript"] <- "no"
            PTCs.df$Raw_low_expressed_transcript <- "yes"
            PTCs.df[PTCs.df$transcript_id %in% rownames(RNAseq.GTEx.raw.filt),"Raw_low_expressed_transcript"] <- "no"
            # 2.10) PCs for transcripts, by tissue
            PCA.path <- paste0(GTEx.RNAseq.PCA.path,GTEx.tissue,"_tissue_PCA_var.txt")
            PCA.df <- read.table(file = PCA.path, header = TRUE, sep = "\t", row.names = 1)
            nPCs <- ncol(PCA.df)
            colnames(PCA.df) <- paste0("tissue_PC",1:nPCs)
            PTCs.df <- merge(PTCs.df, PCA.df, by.x = "transcript_id", by.y = "row.names", all.x = TRUE)
            # 2.11) transcript expression of WT samples from same RNAseq NMF subtype
            # Also add library size from the WT median sample
            # 2.11.1) TPM
            PTCs.df$TPM_transcript_expression_WT <- NA
            df.tmp <- transcripts_WT_expression(PTCs.transcripts.df = PTCs.df[,c("transcript_id","start_pos","TPM_transcript_expression_WT")], GTEx.sample = GTEx.sample,
                                        RNAseq.GTEx = RNAseq.GTEx.TPM, transcripts.germ.mut = transcripts.germ.mut.tissue, RNAseq_type = "TPM")
            PTCs.df$TPM_transcript_expression_WT <- NULL
            PTCs.df <- merge(PTCs.df,df.tmp, by.x = c("start_pos","transcript_id"), by.y = c("start_pos","transcript_id"), all.x = TRUE)
            # 2.11.2) Raw
            PTCs.df$Raw_transcript_expression_WT <- NA
            df.tmp <- transcripts_WT_expression(PTCs.transcripts.df = PTCs.df[,c("transcript_id","start_pos","Raw_transcript_expression_WT")], GTEx.sample = GTEx.sample,
                                        RNAseq.GTEx = RNAseq.GTEx.raw, transcripts.germ.mut = transcripts.germ.mut.tissue, RNAseq_type = "Raw")
            PTCs.df$Raw_transcript_expression_WT <- NULL
            PTCs.df <- merge(PTCs.df,df.tmp, by.x = c("start_pos","transcript_id"), by.y = c("start_pos","transcript_id"), all.x = TRUE)
            # Remove duplicated
            PTCs.df <- PTCs.df[!duplicated(PTCs.df),]
            # 2.12) NMD efficiency
            PTCs.df$NMD_efficiency_TPM <- -log2(PTCs.df$TPM_transcript_expression/PTCs.df$TPM_transcript_expression_WT)
            PTCs.df$NMD_efficiency_Raw <- -log2(PTCs.df$Raw_transcript_expression/PTCs.df$Raw_transcript_expression_WT)
            # 3) Save results
            print(paste0("Writting results to --> ",PTCs.df.metadata.path))
            write.table(PTCs.df, file = PTCs.df.metadata.path,
                        sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
    } else {
        print("Sample not found in RNAseq matrix, skipping...")
    }
}


