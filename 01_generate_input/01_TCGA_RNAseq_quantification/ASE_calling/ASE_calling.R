library("vcfR")
library("stringr")

ASE_calling <- function(TCGA_cancer,TCGA_sample,type) {

    if (type == "germline") {
        # Check if ASE_calling already done
        file_path <- paste0("/g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/primary_tumor/",TCGA_cancer,"/",TCGA_sample,"/ASE_germline/ASE_germline_strelka_calls.txt")
        if (file.exists(file_path)) {
            return(NA)
        }
        # Open VCF germline RNASeq-strelka file
        vcf_RNA_file <- paste0("/g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/primary_tumor/",TCGA_cancer,"/",TCGA_sample,"/ASE_germline/VCF_germline_strelka_recall.vcf.gz")
        # Open original VCF DNAseq-Strelka germline variants file that includes the Annovar annotations
        vcf_germline_DNA_file_path <- paste0("realpath /g/strcombio/fsupek_cancer1/gpalou/TCGA_germline_variants/",gsub("TCGA-","",TCGA_cancer),"/",TCGA_sample,"/*/*/*/*/*_annotated.HG_multianno.txt.gz")
        tryCatch({
            error <- FALSE
            vcf_germline_DNA_file_path <- system(vcf_germline_DNA_file_path, intern = TRUE)
            vcf_germline_DNA <- read.table(file = gzfile(vcf_germline_DNA_file_path), header = TRUE, sep = "\t")
        },error = function(e) {
            error <<- TRUE
            print("sample with no original annotated germline DNAseq VCF called...")
        })
        if (isTRUE(error)) {
            return(NA)
        }
        # if (file.exists(vcf_germline_DNA_file_path)) {
        #     vcf_germline_DNA <- read.table(file = gzfile(vcf_germline_DNA_file_path), header = TRUE, sep = "\t")
        # } else {
        #     print("sample with no original annotated germline DNAseq VCF called...")
        #     return(NA)
        # }
    } else if (type == "somatic") {
        # Check if ASE_calling already done
        file_path <- paste0("/g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/primary_tumor/",TCGA_cancer,"/",TCGA_sample,"/ASE_somatic/ASE_somatic_strelka_calls.txt")
        # if (file.exists(file_path)) {
        #     return(NA)
        # }
        # Open VCF somatic RNASeq-strelka file
        vcf_RNA_file <- paste0("/g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/primary_tumor/",TCGA_cancer,"/",TCGA_sample,"/ASE_somatic/VCF_somatic_indels_strelka_recall.vcf.gz")
        # Open original somatic VCF DNAseq-Strelka variants SNVs and INDELs files

        tryCatch({
            error <- FALSE
            vcf_somatic_snvs_DNA_file_path <- paste0("realpath /g/strcombio/fsupek_cancer1/gpalou/TCGA_somatic_variants/",gsub("TCGA-","",TCGA_cancer),"/",TCGA_sample,"/*/*/*/*_snvs_annotated.HG_multianno.txt.gz")
            vcf_somatic_snvs_DNA_file_path <- system(vcf_somatic_snvs_DNA_file_path, intern = TRUE)
            vcf_somatic_snvs_DNA <- read.table(file = gzfile(vcf_somatic_snvs_DNA_file_path), header = TRUE, sep = "\t")
            vcf_somatic_indels_DNA_file_path <- paste0("realpath /g/strcombio/fsupek_cancer1/gpalou/TCGA_somatic_variants/",gsub("TCGA-","",TCGA_cancer),"/",TCGA_sample,"/*/*/*/*_indels_annotated.HG_multianno.txt.gz")
            vcf_somatic_indels_DNA_file_path <- system(vcf_somatic_indels_DNA_file_path, intern = TRUE)
            vcf_somatic_indels_DNA <- read.table(file = gzfile(vcf_somatic_indels_DNA_file_path), header = TRUE, sep = "\t")
        },error = function(e) {
            print("sample with no original annotated somatic DNAseq VCF called...")
            error <- TRUE
        })
        if (isTRUE(error)) {
            return(NA)
        }
        # if (file.exists(vcf_somatic_snvs_DNA_file_path)) {
        #     vcf_somatic_snvs_DNA <- read.table(file = gzfile(vcf_somatic_snvs_DNA_file_path), header = TRUE, sep = "\t")
        # } else {
        #     print("sample with no original annotated somatic SNVs DNAseq VCF called...")
        #     return(NA)
        # }
        # if (file.exists(vcf_somatic_indels_DNA_file_path)) {
        #     vcf_somatic_indels_DNA <- read.table(file = gzfile(vcf_somatic_indels_DNA_file_path), header = TRUE, sep = "\t")
        # } else {
        #     print("sample with no original annotated somatic INDELs DNAseq VCF called...")
        #     return(NA)
        # }
        # Merge
        vcf_somatic_DNA <- rbind(vcf_somatic_snvs_DNA,vcf_somatic_indels_DNA)
    }

    # 1) Calculate ASE calls for RNA-seq
    if (file.exists(vcf_RNA_file)) {
        tryCatch({
        error <- FALSE
        vcf_RNA <- read.vcfR(vcf_RNA_file, verbose = FALSE ) 
        },error = function(e) {
            print(e)
            error <<- TRUE
        })
        if (isTRUE(error)) {
            return(NA)
        }
    } else {
        print("sample with no germline RNAseq VCF called...")
        return(NA)
    }
    # Convert to dataframes
    vcf_RNA_fix <- data.frame(vcf_RNA@fix)
    vcf_RNA_gt <- data.frame(vcf_RNA@gt)

    # Obtain Allelic Depth and Genome Quality for each variant
    # Calculate ASE using AD
    count <- 0
    error <- c()

    df_ASE_RNAseq <- apply(vcf_RNA_gt,1,function(variant) {
        count <<- count + 1
        error <- variant
        vcf_format_variant <- unlist(strsplit(as.character(variant[["FORMAT"]]),":"))
        vcf_gt_variant <- unlist(strsplit(as.character(variant[["sm1"]]),":"))
        # Check position of AD, GT, GQ column (is different for SNVs than indels)
        vcf_AD_pos <- which(vcf_format_variant=="AD")
        vcf_GQ_pos <- which(vcf_format_variant=="GQ")
        vcf_GT_pos <- which(vcf_format_variant=="GT")
        # Obtain AD, GQ, GT
        vcf_AD_reads <- vcf_gt_variant[vcf_AD_pos]
        vcf_GQ <- as.numeric(vcf_gt_variant[vcf_GQ_pos])
        vcf_GT <- vcf_gt_variant[vcf_GT_pos]
        reads <- unlist(strsplit(vcf_AD_reads,","))
        if (is.null(reads)) {
            variant[["refCount"]] <- NA
            variant[["altCount"]] <- NA
            variant[["totalCount"]] <- NA
        } else {
            variant[["refCount"]] <- as.numeric(reads[1])
            variant[["altCount"]] <- as.numeric(reads[2])
            variant[["totalCount"]] <- as.numeric(reads[1]) + as.numeric(reads[2])
        }
        if (is.na(vcf_GQ) || length(vcf_GQ) == 0) {
            variant[["GQ"]] <- NA
        } else {
            variant[["GQ"]] <- vcf_GQ
        }
        if (is.na(vcf_GT) || length(vcf_GT) == 0) {
            variant[["GT"]] <- NA
        } else {
            variant[["GT"]] <- vcf_GT
        }
        variant
    })

    df_ASE_RNAseq_final <- data.frame(t(df_ASE_RNAseq))
    #df_ASE_final <- data.frame(do.call(rbind,df_ASE))
    # Add rest of VCF information
    VCF_ASE_RNAseq <- cbind(vcf_RNA_fix,df_ASE_RNAseq_final)
    # Keep only heterozygous variants for ASE
    VCF_ASE_RNAseq <- VCF_ASE_RNAseq[VCF_ASE_RNAseq$GT == "0/1",]

    # 2) Calculate ASE calls for DNA-seq
    if (type == "germline") {
        # Obtain Allelic Depth and Genome Quality for each variant
        # Calculate ASE using AD
        df_ASE_DNAseq <- apply(vcf_germline_DNA,1,function(variant) {
            vcf_format_variant <- unlist(strsplit(as.character(variant[["Otherinfo12"]]),":"))
            vcf_gt_variant <- unlist(strsplit(as.character(variant[["Otherinfo13"]]),":"))
            # Check position of AD, GT, GQ column (is different for SNVs than indels)
            vcf_AD_pos <- which(vcf_format_variant=="AD")
            vcf_GQ_pos <- which(vcf_format_variant=="GQ")
            vcf_GT_pos <- which(vcf_format_variant=="GT")
            # Obtain AD, GQ and GT
            vcf_AD_reads <- vcf_gt_variant[vcf_AD_pos]
            vcf_GQ <- as.numeric(vcf_gt_variant[vcf_GQ_pos])
            vcf_GT <- vcf_gt_variant[vcf_GT_pos]
            reads <- unlist(strsplit(vcf_AD_reads,","))
            if (is.null(reads)) {
                variant[["refCount"]] <- NA
                variant[["altCount"]] <- NA
                variant[["totalCount"]] <- NA
            } else {
                variant[["refCount"]] <- as.numeric(reads[1])
                variant[["altCount"]] <- as.numeric(reads[2])
                variant[["totalCount"]] <- as.numeric(reads[1]) + as.numeric(reads[2])
            }
            if (is.na(vcf_GQ)) {
                variant[["GQ"]] <- NA
            } else {
                variant[["GQ"]] <- vcf_GQ
            }
            if (is.na(vcf_GT)) {
                variant[["GT"]] <- NA
            } else {
                variant[["GT"]] <- vcf_GT
            }
            variant
        })
        VCF_ASE_DNAseq <- data.frame(t(df_ASE_DNAseq), stringsAsFactors = FALSE)
        #df_ASE_final <- data.frame(do.call(rbind,df_ASE))
        colnames(VCF_ASE_DNAseq) <- c("CHROM","START","END","REF_DNA","ALT_DNA","Func.refGene","Gene.refGene","GeneDetail.refGene","ExonicFunc.refGene","AAChange.refGene",
            "AF","AF_popmax","AF_male","AF_female","AF_raw","AF_afr","AF_sas","AF_amr","AF_eas","AF_nfe","AF_fin","AF_asj","AF_oth","non_topmed_AF_popmax",
            "non_neuro_AF_popmax","non_cancer_AF_popmax","controls_AF_popmax","ExAC_ALL","ExAC_AFR","ExAC_AMR","ExAC_EAS","ExAC_FIN","ExAC_NFE","ExAC_OTH","ExAC_SAS",
            "Otherinfo1_DNA","Otherinfo2_DNA","Otherinfo3_DNA","Otherinfo4_DNA","Otherinfo5_DNA","Otherinfo6_DNA","Otherinfo7_DNA","Otherinfo8_DNA","Otherinfo9_DNA",
            "Otherinfo10_DNA","Otherinfo11_DNA","Otherinfo12_DNA","Otherinfo13_DNA","refCount_DNA","altCount_DNA","totalCount_DNA","GQ_DNA","GT_DNA")
        VCF_ASE_DNAseq[,c("START","END","Otherinfo5_DNA")] <- lapply(VCF_ASE_DNAseq[,c("START","END","Otherinfo5_DNA")] , as.character)
        VCF_ASE_DNAseq[,c("START","END","Otherinfo5_DNA")] <- lapply(VCF_ASE_DNAseq[,c("START","END","Otherinfo5_DNA")] , as.numeric)
        # Keep only heterozygous variants for ASE
        VCF_ASE_DNAseq <- VCF_ASE_DNAseq[VCF_ASE_DNAseq$GT == "0/1",]
    } else if (type == "somatic") {
        nt_index <- function(nt){
            if (nt == "A"){
                return(1)
            } else if (nt == "C") {
                return(2)
            } else if (nt == "G") {
                return(3)
            } else {
                return(4)
            }
        }
        # Retrieve EVS score for each variant
        df_ASE_DNAseq <- apply(vcf_somatic_DNA,1,function(variant) {
            vcf_info_variant <- unlist(strsplit(as.character(variant[["Otherinfo11"]]),";"))
            vcf_format_variant <- unlist(strsplit(as.character(variant[["Otherinfo12"]]),":"))
            # Otherinfo13 is NORMAL, Otherinfo14 is TUMOR
            vcf_gt_variant <- unlist(strsplit(as.character(variant[["Otherinfo14"]]),":"))
            # Check position of GT and EVS score
            EVS_score_pos <- grep("SomaticEVS.*",vcf_info_variant)
            vcf_GT_pos <- which(vcf_format_variant=="GT")
            # Obtain GT and EVS score
            EVS_score <- as.numeric(gsub("SomaticEVS=","",vcf_info_variant[EVS_score_pos]))
            vcf_GT <- vcf_gt_variant[vcf_GT_pos]
            if (is.null(EVS_score)) {
                variant[["EVS_score"]] <- NA
            } else {
                variant[["EVS_score"]] <- EVS_score
            }
            if (is.na(vcf_GT)) {
                variant[["GT"]] <- NA
            } else {
                variant[["GT"]] <- vcf_GT
            }
            # Obtain DNA ASE calls
            refCount_DNA <- as.numeric(str_split(str_split(variant[["Otherinfo14"]], ":")[[1]][5+nt_index(as.character(variant[["Ref"]]))], ",")[[1]][1])
            altCount_DNA <- as.numeric(str_split(str_split(variant[["Otherinfo14"]], ":")[[1]][5+nt_index(as.character(variant[["Alt"]]))], ",")[[1]][1])
            variant[["refCount_DNA"]] <- as.numeric(refCount_DNA)
            variant[["altCount_DNA"]] <- as.numeric(altCount_DNA)
            variant[["totalCount_DNA"]] <- as.numeric(refCount_DNA) + as.numeric(altCount_DNA)
            # Output
            variant
        })
        VCF_ASE_DNAseq <- data.frame(t(df_ASE_DNAseq), stringsAsFactors = FALSE)
        colnames(VCF_ASE_DNAseq) <- c("CHROM","START","END","REF_DNA","ALT_DNA","Func.refGene","Gene.refGene","GeneDetail.refGene","ExonicFunc.refGene","AAChange.refGene",
            "AF","AF_popmax","AF_male","AF_female","AF_raw","AF_afr","AF_sas","AF_amr","AF_eas","AF_nfe","AF_fin","AF_asj","AF_oth","non_topmed_AF_popmax",
            "non_neuro_AF_popmax","non_cancer_AF_popmax","controls_AF_popmax","ExAC_ALL","ExAC_AFR","ExAC_AMR","ExAC_EAS","ExAC_FIN","ExAC_NFE","ExAC_OTH","ExAC_SAS",
            "Otherinfo1_DNA","Otherinfo2_DNA","Otherinfo3_DNA","Otherinfo4_DNA","Otherinfo5_DNA","Otherinfo6_DNA","Otherinfo7_DNA","Otherinfo8_DNA","Otherinfo9_DNA",
            "Otherinfo10_DNA","Otherinfo11_DNA","Otherinfo12_DNA","Otherinfo13_DNA","Otherinfo14_DNA","EVS_score","GT_DNA","refCount_DNA","altCount_DNA","totalCount_DNA")
        VCF_ASE_DNAseq[,c("START","END","Otherinfo5_DNA")] <- lapply(VCF_ASE_DNAseq[,c("START","END","Otherinfo5_DNA")] , as.character)
        VCF_ASE_DNAseq[,c("START","END","Otherinfo5_DNA")] <- lapply(VCF_ASE_DNAseq[,c("START","END","Otherinfo5_DNA")] , as.numeric)
        # Keep only heterozygous variants for ASE
        VCF_ASE_DNAseq <- VCF_ASE_DNAseq[VCF_ASE_DNAseq$GT_DNA == "0/1",]
    }
    # 3) Filter
    colnames(VCF_ASE_RNAseq) <- c("CHROM","POS","ID_RNA","REF","ALT","QUAL_RNA","FILTER_RNA","INFO_RNA","FORMAT_RNA","sample_GT_RNA","refCount_RNA","altCount_RNA","totalCount_RNA","GQ_RNA","GT_RNA")
    VCF_ASE_RNAseq$POS <- as.numeric(as.character(VCF_ASE_RNAseq$POS))
    # Filter by positions from the original VCF strelka file
    VCFs_merged <- merge(VCF_ASE_RNAseq, VCF_ASE_DNAseq, by.x = c("CHROM","POS","REF","ALT"), by.y = c("CHROM","Otherinfo5_DNA","Otherinfo7_DNA","Otherinfo8_DNA"))

    return(VCFs_merged)
    
}

args <- commandArgs(trailingOnly=TRUE)
TCGA_cancer <- args[1]
TCGA_sample <- args[2]
type <- args[3]
print(paste0("TCGA CANCER --> ",TCGA_cancer))
print(paste0("TCGA sample --> ",TCGA_sample))
print(paste0("VCF TYPE --> ",type))

# 1) GERMLINE
if (type == "germline") {
    VCF_germline_ASE <- ASE_calling(TCGA_cancer=TCGA_cancer, TCGA_sample=TCGA_sample, type=type)
    if (is.na(VCF_germline_ASE)) {
        print("Sample germline ASE calling already done or error...")
    } else {
        # Save output
        output_file <- paste0("/g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/primary_tumor/",TCGA_cancer,"/",TCGA_sample,"/ASE_germline/ASE_germline_strelka_calls.txt")
        print(paste0("Output file --> ",output_file))
        write.table(VCF_germline_ASE, file = output_file, 
                            quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
    } 
} else if (type == "somatic") { # 2) SOMATIC
    VCF_somatic_ASE <- ASE_calling(TCGA_cancer = TCGA_cancer, TCGA_sample = TCGA_sample, type=type)
    if (is.na(VCF_somatic_ASE)) {
        print("Sample somatic ASE calling already done or error...")
    } else {
        # Save output
        output_file <- paste0("/g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/primary_tumor/",TCGA_cancer,"/",TCGA_sample,"/ASE_somatic/ASE_somatic_strelka_calls.txt")
        print(paste0("Output file --> ",output_file))
        write.table(VCF_somatic_ASE, file = output_file,
                            quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
    }
}




