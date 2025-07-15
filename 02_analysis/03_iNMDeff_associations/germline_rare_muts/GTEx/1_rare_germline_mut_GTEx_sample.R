rm(list=ls())

################################################################################################

########################################## FUNCTIONS ###########################################

################################################################################################

VCF_germline_sample_file <- function(GTEx_sample) {
    # SNVs
    tryCatch({
        error <- FALSE
        vcf_germline_DNA_file_path <- paste0("realpath /g/strcombio/fsupek_cancer1/gpalou/GTEx/WES_VCF/VCFs/",GTEx_sample,"_liftOver_annotated.HG_multianno.txt.gz")
        vcf_germline_DNA_file_path <- system(vcf_germline_DNA_file_path, intern = TRUE)
        vcf_germline_DNA <- read.csv(file = gzfile(vcf_germline_DNA_file_path), header = TRUE, sep = "\t")
    },error = function(e) {
        print("sample with no original annotated germline DNAseq VCF called...")
        error <- TRUE
    })
    # Fix MAF numeric
    MAFs <- vcf_germline_DNA$AF
    MAFs <- ifelse(MAFs == "" | MAFs == ".",NA,format(MAFs,scientific = TRUE))
    MAFs <- as.numeric(MAFs)
    vcf_germline_DNA$AF <- MAFs
    # Fix ENSEMBL Genes IDs
    # Some transcript swill have different Genes IDs (for whatever reason), but bc it's only used on the regression as gene.id it's not problem
    vcf_germline_DNA$Gene.refGene <- gsub("(ENSG\\d{11})\\.\\d*","\\1",vcf_germline_DNA$Gene.refGene)
    return(vcf_germline_DNA)
}

VCF_germline_sample_filtering <- function(VCF_germline_sample, mutations, MAF = NULL, filter = NULL) {
    splicing_filter <- which(VCF_germline_sample$Func.refGene %in% c("splicing","exonic;splicing"))
    other_filter <- which(VCF_germline_sample$ExonicFunc.refGene %in% mutations)
    mut_filter <- unique(c(splicing_filter,other_filter))
    VCF_germline_sample_filt <- VCF_germline_sample[mut_filter,]
    # MAF
    if (!is.null(MAF)) {
      VCF_germline_sample_filt <- VCF_germline_sample_filt[which(VCF_germline_sample_filt$AF <= MAF),]
    }
    # PASS variants
    if (!is.null(filter)) {
      VCF_germline_sample_filt <- VCF_germline_sample_filt[VCF_germline_sample_filt$Otherinfo10 == filter,]
    }
    return(VCF_germline_sample_filt)
}

# VCF_extract_info <- function(VCF) {
#     variants.list <- list()
#     for (i in 1:nrow(VCF)) {
#       # Transcripts list from same variant
#       transcripts.list <- strsplit(as.character(VCF[i,"AAChange.refGene"]), split=",")
#       # Obtain the ID name and other info (exon location, PTC CDS location)
#       transcripts.exons.info <- lapply(transcripts.list, function(transcript) {
#         variant.CDS.exon.num <- gsub(".*ENST[0-9]{11}\\.[0-9]{1,3}:exon([0-9]{1,3}).*","\\1",transcript)
#         transcript.id <- gsub(".*(ENST[0-9]{11})\\.[0-9]{1,3}:exon.*","\\1",transcript)
#         variant.CDS.pos <- gsub(".*ENST[0-9]{11}\\.[0-9]{1,3}:exon[0-9]{1,3}:c\\.\\D*(\\d*)_?.*:p\\..*","\\1",transcript)
#         data.frame(transcript_id=transcript.id,exon_num=variant.CDS.exon.num, variant_CDS_pos = variant.CDS.pos)
#       })
#       transcripts.exons.info <- transcripts.exons.info[[1]]
#       transcripts.exons.info$start_pos <- VCF[i,"Start"]
#       transcripts.exons.info$gene_id <- VCF[i,"Gene.refGene"]
#       transcripts.exons.info$MAF <- VCF[i,"AF"]
#       transcripts.exons.info$chr <- VCF[i,"Chr"]
#       if ( VCF[i,"Func.refGene"] %in% c("exonic;splicing","splicing","ncRNA_exonic;splicing","ncRNA_splicing") ) {
#         transcripts.exons.info$variant_classification <- "splicing"
#       } else {
#         transcripts.exons.info$variant_classification <- VCF[i,"ExonicFunc.refGene"]
#       }
#       transcripts.exons.info$filter <- VCF[i,"Otherinfo10"]
#       variants.list[[i]] <- transcripts.exons.info
#   }
#   return(variants.list)
# }

################################################################################################

########################################## LIBRARIES ###########################################

################################################################################################

################################################################################################

########################################## SCRIPT ##############################################

################################################################################################

# 1) Obtain data

# Arguments and paths

args <- commandArgs(trailingOnly=TRUE)
#paths <- read.table(file = args[1], header = TRUE, sep = ",", stringsAsFactors = FALSE)
GTEx_sample <- args[1]

# 1) Check rare germline mutations
# 1.1) Obtain germline VCF
print(paste0("Starting germline mut search for TCGA sample --> ",GTEx_sample))
mutations <- c("stopgain","nonsynonymous SNV","startloss","stoploss","frameshift deletion",
                "frameshift insertion","nonframeshift deletion","nonframeshift insertion")
truncating_mut <- "stopgain|^frameshift deletion|^frameshift insertion|splicing"
miss_mut <- "nonsynonymous SNV|startloss|stoploss|nonframeshift deletion|nonframeshift insertion"
# Check if file exists
VCF_germline_sample <- VCF_germline_sample_file(GTEx_sample = GTEx_sample)
if(is.null(VCF_germline_sample)){
    print("VCF germline file not found, skipping...")
} else {
    # 2.2) Filter VCF
    VCF_germline_sample_filt <- VCF_germline_sample_filtering(VCF_germline_sample = VCF_germline_sample, 
                                                                mutations = mutations, MAF = 0.001, filter = "PASS")
    if (nrow(VCF_germline_sample_filt) == 0) {
        print("No variants in this sample")
    } else {
        # 2.3) Obtain variant-per-transcript information 
        # variants_transcript <- VCF_extract_info(VCF = VCF_germline_sample_filt)
        # variants_transcript_df <- do.call(rbind,variants_transcript)
    }
}

# Save output table
sample_germline_mut_path <- paste0("/g/strcombio/fsupek_cancer1/gpalou/GTEx/WES_VCF/rare_germline_variants/samples/",GTEx_sample,"_rare_germline_variants.txt")
if (!is.null(sample_germline_mut_path)) {
    print(paste0("Writting output at --> ",sample_germline_mut_path))
    write.table(VCF_germline_sample_filt, file = sample_germline_mut_path,
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
}

# Check if in VCF we have the same rare variants as in our TXT filtered manually file

VCF <- read.table(file = paste0("/g/strcombio/fsupek_cancer1/gpalou/GTEx/WES_VCF/VCFs/",GTEx_sample,"_liftOver_annotated.HG_multianno_rare_variants.vcf.gz"),
            comment.char = "", header = TRUE, skip = 825, sep = "\t")
dim(VCF)
dim(VCF_germline_sample_filt)
VCF_germline_sample_filt[!VCF_germline_sample_filt$Start %in% VCF$POS,"ExonicFunc.refGene"]

if (nrow(VCF) != nrow(VCF_germline_sample_filt)) {
    df <- data.frame(VCF = nrow(VCF), TXT = nrow(VCF_germline_sample_filt))
    write.table(paste0("/g/strcombio/fsupek_cancer1/gpalou/GTEx/WES_VCF/rare_germline_variants/samples/",GTEx_sample,"_error.log"), file = df,
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
}

