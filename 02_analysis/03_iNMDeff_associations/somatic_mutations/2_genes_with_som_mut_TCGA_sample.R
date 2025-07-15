rm(list=ls())

################################################################################################

########################################## FUNCTIONS ###########################################

################################################################################################

VCF_somatic_sample_file <- function(TCGA_cancer, TCGA_sample) {
    # SNVs
    tryCatch({
        sample_rand <- sample(1:2)[1]
        error_snvs <- FALSE
        vcf_somatic_snvs_DNA_file_path <- paste0("realpath /g/strcombio/fsupek_cancer1/gpalou/TCGA_somatic_variants/",gsub("TCGA-","",TCGA_cancer),"/",TCGA_sample,"/*/*/*/*_snvs_annotated.HG_multianno.txt.gz")
        vcf_somatic_snvs_DNA_file_path <- system(vcf_somatic_snvs_DNA_file_path, intern = TRUE)
        if (length(vcf_somatic_snvs_DNA_file_path) > 1) {
            vcf_somatic_snvs_DNA_file_path <- vcf_somatic_snvs_DNA_file_path[sample_rand]
        }
        vcf_somatic_snvs_DNA <- read.table(file = gzfile(vcf_somatic_snvs_DNA_file_path), header = TRUE, sep = "\t")
    },error = function(e) {
        print("sample with no original annotated somatic DNAseq VCF called...")
        error_snvs <- TRUE
    })
    # INDELs
    tryCatch({
        error_indels <- FALSE
        vcf_somatic_indels_DNA_file_path <- paste0("realpath /g/strcombio/fsupek_cancer1/gpalou/TCGA_somatic_variants/",gsub("TCGA-","",TCGA_cancer),"/",TCGA_sample,"/*/*/*/*_indels_annotated.HG_multianno.txt.gz")
        vcf_somatic_indels_DNA_file_path <- system(vcf_somatic_indels_DNA_file_path, intern = TRUE)
        if (length(vcf_somatic_indels_DNA_file_path) > 1) {
            vcf_somatic_indels_DNA_file_path <- vcf_somatic_indels_DNA_file_path[sample_rand]
        }
        vcf_somatic_indels_DNA <- read.table(file = gzfile(vcf_somatic_indels_DNA_file_path), header = TRUE, sep = "\t")
    },error = function(e) {
        print("sample with no original annotated somatic DNAseq VCF called...")
        error_indels <- TRUE
    })
    # # Merge INDELs and/or SNVs
    if (isTRUE(error_snvs)) {
        vcf_somatic_DNA <- vcf_somatic_indels_DNA
    } else if (isTRUE(error_indels)) {
        vcf_somatic_DNA <- vcf_somatic_snvs_DNA
    } else if ( isTRUE(error_snvs) & isTRUE(error_indels) ) {
        return(NA)
    } else {
        vcf_somatic_DNA <- rbind(vcf_somatic_snvs_DNA,vcf_somatic_indels_DNA)
    }
    # Fix VAF numeric
    VAFs <- vcf_somatic_DNA$AF
    VAFs <- ifelse(VAFs == "" | VAFs == ".",0,format(VAFs,scientific = TRUE))
    VAFs <- as.numeric(VAFs)
    vcf_somatic_DNA$AF <- VAFs
    # Fix ENSEMBL Genes IDs
    # Some transcript swill have different Genes IDs (for whatever reason), but bc it's only used on the regression as gene.id it's not problem
    vcf_somatic_DNA$Gene.refGene <- gsub("(ENSG\\d{11})\\.\\d*","\\1",vcf_somatic_DNA$Gene.refGene)
    return(vcf_somatic_DNA)
}

VCF_somatic_sample_filtering <- function(VCF_somatic_sample, mutations, VAF = NULL, filter = NULL) {
    splicing_filter <- which(VCF_somatic_sample$Func.refGene %in% c("splicing","exonic;splicing"))
    other_filter <- which(VCF_somatic_sample$ExonicFunc.refGene %in% mutations)
    mut_filter <- unique(c(splicing_filter,other_filter))
    VCF_somatic_sample_filt <- VCF_somatic_sample[mut_filter,]
    # MAF
    if (!is.null(VAF)) {
      VCF_somatic_sample_filt <- VCF_somatic_sample_filt[which(VCF_somatic_sample_filt$AF <= MAF),]
    }
    # PASS variants
    if (!is.null(filter)) {
      VCF_somatic_sample_filt <- VCF_somatic_sample_filt[VCF_somatic_sample_filt$Otherinfo10 == filter,]
    }
    return(VCF_somatic_sample_filt)
}

VCF_extract_info <- function(VCF) {
    variants.list <- list()
    for (i in 1:nrow(VCF)) {
      # Transcripts list from same variant
      transcripts.list <- strsplit(as.character(VCF[i,"AAChange.refGene"]), split=",")
      # Obtain the ID name and other info (exon location, PTC CDS location)
      transcripts.exons.info <- lapply(transcripts.list, function(transcript) {
        variant.CDS.exon.num <- gsub(".*ENST[0-9]{11}\\.[0-9]{1,3}:exon([0-9]{1,3}).*","\\1",transcript)
        transcript.id <- gsub(".*(ENST[0-9]{11})\\.[0-9]{1,3}:exon.*","\\1",transcript)
        variant.CDS.pos <- gsub(".*ENST[0-9]{11}\\.[0-9]{1,3}:exon[0-9]{1,3}:c\\.\\D*(\\d*)_?.*:p\\..*","\\1",transcript)
        data.frame(transcript_id=transcript.id,exon_num=variant.CDS.exon.num, variant_CDS_pos = variant.CDS.pos)
      })
      transcripts.exons.info <- transcripts.exons.info[[1]]
      transcripts.exons.info$start_pos <- VCF[i,"Start"]
      transcripts.exons.info$gene_id <- VCF[i,"Gene.refGene"]
      transcripts.exons.info$VAF <- VCF[i,"AF"]
      transcripts.exons.info$chr <- VCF[i,"Chr"]
      if ( VCF[i,"Func.refGene"] %in% c("exonic;splicing","splicing","ncRNA_exonic;splicing","ncRNA_splicing") ) {
        transcripts.exons.info$variant_classification <- "splicing"
      } else {
        transcripts.exons.info$variant_classification <- VCF[i,"ExonicFunc.refGene"]
      }
      transcripts.exons.info$filter <-VCF[i,"Otherinfo10"]
      transcripts.exons.info$EVS <- as.numeric(gsub(".*SomaticEVS=([0-9]*)","\\1",VCF[i,"Otherinfo11"]))
      variants.list[[i]] <- transcripts.exons.info
  }
  return(variants.list)
}

################################################################################################

########################################## LIBRARIES ###########################################

################################################################################################



################################################################################################

########################################## SCRIPT ##############################################

################################################################################################

# 1) Obtain data

# Arguments and paths

args <- commandArgs(trailingOnly=TRUE)
paths <- read.table(file = args[1], header = TRUE, sep = ",", stringsAsFactors = FALSE)
TCGA_cancer <- args[2]
TCGA_sample <- args[3]

print(paste0("TCGA cancer --> ",TCGA_cancer))

paths <- read.table(file = "/home/gpalou/projects/NMD/scripts/NMD_efficiency/analysis_results/somatic_associations/2_genes_with_som_mut_TCGA_sample.txt", header = TRUE, sep = ",", stringsAsFactors = FALSE)
conversor_tables_path <- paths[paths$folder_or_object=="conversor_tables","path_or_filename"]
# transcripts_to_remove_path <- paths[paths$folder_or_object=="transcripts_to_remove_path","path_or_filename"]
transcripts_germ_mut_path <- paths[paths$folder_or_object=="transcripts_germ_mut_path","path_or_filename"]
CNV_path <- paths[paths$folder_or_object=="CNV_file_path","path_or_filename"]
CGC_somatic_mut_CNV_path <- paths[paths$folder_or_object=="CGC_somatic_mut_CNV_path","path_or_filename"]
TMB_TCGA_path <- paths[paths$folder_or_object=="TMB_TCGA_path","path_or_filename"]

# 1.1) Obtain Cancer Gene Census table from COSMIC
CGC <- read.csv(file = "/g/strcombio/fsupek_cancer1/gpalou/COSMIC/cancer_gene_census_updated.tsv", 
                    header = TRUE, stringsAsFactors = FALSE)
# 1.2) Firehose subtypes (mRNA NMF) and CNV
if (TCGA_cancer == "TCGA-SKCM") {
  CNV_file_path <- gsub("-TP","-TM",paste0(CNV_path,paths[paths$folder_or_object=="CNV_file","path_or_filename"]))
} else {
  CNV_file_path <- paste0(CNV_path,paths[paths$folder_or_object=="CNV_file","path_or_filename"])
}
# CNV
CNV_file <- read.table(file = gsub("\\[X\\]",gsub("TCGA-","",TCGA_cancer), CNV_file_path ), 
                        header = TRUE, sep = "\t", stringsAsFactors = FALSE)
colnames(CNV_file) <- substr(colnames(CNV_file),1,12)
# Remove ENSG from Gene Symbols
CNV_file$Gene.Symbol <- gsub("(.*)\\|ENSG.*","\\1",CNV_file$Gene.Symbol)

# 2) Check somatic mutations in CGC
# 2.1) Obtain somatic VCF
print(paste0("Starting somatic mut search for TCGA sample --> ",TCGA_sample))
mutations <- c("synonymous SNV","stopgain","nonsynonymous SNV","startloss","stoploss","frameshift deletion",
                "frameshift insertion","nonframeshift deletion","nonframeshift insertion")
truncating_mut <- "stopgain|^frameshift deletion|^frameshift insertion|splicing"
miss_mut <- "nonsynonymous SNV|startloss|stoploss|nonframeshift deletion|nonframeshift insertion"
# Check if file exists
VCF_somatic_sample <- VCF_somatic_sample_file(TCGA_cancer = TCGA_cancer, TCGA_sample = TCGA_sample)
if(is.null(VCF_somatic_sample)){
    print("VCF somatic file not found, skipping...")
} else {
    # 2.2) Filter VCF
    VCF_somatic_sample_filt <- VCF_somatic_sample_filtering(VCF_somatic_sample = VCF_somatic_sample, 
                                mutations = mutations, VAF = 1, filter = "PASS")
    if (nrow(VCF_somatic_sample_filt) == 0) {
        print("No variants in this sample")
    } else {
        # 2.3) Obtain variant-per-transcript information 
        variants_transcript <- VCF_extract_info(VCF = VCF_somatic_sample_filt)
        variants_transcript_df <- do.call(rbind,variants_transcript)
        #genes_mut <- unique(ensembl_gene_v88[ensembl_gene_v88$transcript_id %in% variants_transcript_df$transcript_id,"gene_id"])
        # 2.4) Obtain mutated ENSEMBL genes from the VCF
        list_genes_mut <- strsplit(variants_transcript_df$gene_id,";")
        genes_mut <- unlist(list_genes_mut)
        # 2.5) Check which CGC genes overlap with somatic mutations
        filter <- CGC$Synonyms %in% genes_mut
        if (sum(filter) != 0) {
            CGC_mut_sample <- CGC[filter,]
            # Add transcripts information for the mutated genes
            CGC_mut_sample[,c("somatic_EVS","somatic_filter","somatic_VAF","somatic_variant_classification","somatic_mut_transcript")] <- NA
        } else {
            CGC_mut_sample <- data.frame()
        }
        if (nrow(CGC_mut_sample) == 0) {
            print("sample with no somatic mutations in CGC genes")
        } else {
            for (i in 1:nrow(CGC_mut_sample)) {
                # Obtain gene ID
                gene_mutated <- CGC_mut_sample[i,"Synonyms"]
                # Check which transcripts have somatic mutations
                genes_index <- grep(gene_mutated,list_genes_mut)
                variants_transcript_df_filt <- variants_transcript_df[genes_index,]
                # Check type of mutations
                truncating_mut_index <- grep(truncating_mut,variants_transcript_df_filt$variant_classification)
                miss_mut_index <- grep(miss_mut,variants_transcript_df_filt$variant_classification)
                # Save transcript information in CGC_mut_sample table
                CGC_mut_sample[i,"somatic_VAF"] <- paste0(variants_transcript_df_filt$VAF,collapse=";")
                CGC_mut_sample[i,"somatic_EVS"] <- paste0(variants_transcript_df_filt$EVS,collapse=";")
                CGC_mut_sample[i,"somatic_filter"] <- paste0(variants_transcript_df_filt$filter,collapse=";")
                CGC_mut_sample[i,"somatic_variant_classification"] <- paste0(variants_transcript_df_filt$variant_classification,collapse=";")
                CGC_mut_sample[i,"somatic_mut_transcript"] <- paste0(variants_transcript_df_filt$transcript_id,collapse=";")
                if (length(truncating_mut_index) != 0) {
                    CGC_mut_sample[i,"somatic_mut_truncating"] <- "yes"
                }
                if (length(miss_mut_index) != 0) {
                    CGC_mut_sample[i,"somatic_mut_missense"] <- "yes"
                }  
            }
        }
    }
}

# 3) Check somatic CNV in CGC genes
# 3.1) Obtain Genes Symbols with somatic CNV
print(paste0("Starting somatic CNV search for TCGA sample --> ",TCGA_sample))
TCGA_sample_updated <- gsub("-",".",TCGA_sample)
col_index <- which(colnames(CNV_file) %in% TCGA_sample_updated)
if (length(col_index) == 0 ) {
    print("Sample with no CNV file")
    CGC_CNV_sample <- data.frame()
} else {
    CNV_file_sample <- CNV_file[,c(1:3,col_index), drop = FALSE]
    genes_mut_index <- which(CNV_file_sample[,4] != 0)
    genes_mut_amp_index <- which(CNV_file_sample[,4] > 0)
    genes_mut_del_index <- which(CNV_file_sample[,4] < 0)
    genes_mut <- CNV_file_sample[genes_mut_index, "Gene.Symbol"]
    genes_mut_amp <- CNV_file_sample[genes_mut_amp_index, "Gene.Symbol"]
    genes_mut_del <- CNV_file_sample[genes_mut_del_index, "Gene.Symbol"]
    # 3.2) Check which CGC Genes Symbols overlap with somatic CNV
    filter <- CGC$Gene.Symbol %in% genes_mut
    CGC_CNV_sample <- CGC[filter,,drop = FALSE]
    if (nrow(CGC_CNV_sample) == 0) {
        print("sample with no somatic CNV in CGC genes")
    } else {
        # Amplifications
        CGC_CNV_sample[which(CGC_CNV_sample$Gene.Symbol%in%genes_mut_amp),"somatic_CNV"] <- "amp"
        # Deletions
        CGC_CNV_sample[which(CGC_CNV_sample$Gene.Symbol%in%genes_mut_del),"somatic_CNV"] <- "del"
    }
}

# 4) Merge CGC somatic CNV and mut
if (nrow(CGC_CNV_sample) == 0 & nrow(CGC_mut_sample) != 0) {
    CGC_somatic_mut_CNV_sample <- CGC_mut_sample
} else if (nrow(CGC_CNV_sample) != 0 & nrow(CGC_mut_sample) != 0) {
    common_col <- colnames(CGC_mut_sample)[which(colnames(CGC_mut_sample)%in%colnames(CGC_CNV_sample))]
    CGC_somatic_mut_CNV_sample <- merge(CGC_mut_sample, CGC_CNV_sample, by.x = common_col, by.y = common_col, all = TRUE)
} else if (nrow(CGC_CNV_sample) != 0 & nrow(CGC_mut_sample) == 0) {
    CGC_somatic_mut_CNV_sample <- CGC_CNV_sample
} else if (nrow(CGC_CNV_sample) == 0 & nrow(CGC_mut_sample) == 0) {
    CGC_somatic_mut_CNV_sample <- NULL
}

# Save output table
CGC_somatic_mut_CNV_path <- gsub("\\[X1\\]",gsub("TCGA-","",TCGA_cancer),CGC_somatic_mut_CNV_path)
CGC_somatic_mut_CNV_path <- paste0(CGC_somatic_mut_CNV_path,paths[paths$folder_or_object=="CGC_somatic_mut_CNV","path_or_filename"])
CGC_somatic_mut_CNV_path <- gsub("\\[X2\\]",TCGA_sample,CGC_somatic_mut_CNV_path)
if (!is.null(CGC_somatic_mut_CNV_sample)) {
    print(paste0("Writting output at --> ",CGC_somatic_mut_CNV_path))
    write.table(CGC_somatic_mut_CNV_sample, file = CGC_somatic_mut_CNV_path,
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
}

