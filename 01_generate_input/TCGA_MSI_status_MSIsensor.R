# TCGA cancer tissues
TCGA_cancer_tissues <- read.table(file = "/g/strcombio/fsupek_home/gpalou/data/TCGA_RNAseq_quantification/TCGA_projects_names.txt", header = FALSE, stringsAsFactors = FALSE)$V1
TCGA_cancer_tissues <- gsub("TCGA-","",TCGA_cancer_tissues)


# Create output Data frame

df_res <- data.frame(TCGA_sample = NA, TCGA_cancer = NA, MSI_score = NA)

# For each cancer
for (cancer in TCGA_cancer_tissues) {

    # 1) Find where the cancer is located on the cluster
    # and obtain the sample names
    fsupek <- 0
    for (i in 1:3) {
        error <- FALSE
        tryCatch( {
            system(paste0("/g/strcombio/fsupek_cancer",i,"/TCGA_bam/",cancer), intern = TRUE)
        },error = function(e) {
            error <<- TRUE
            }
        )
        if (!isTRUE(error)) {
            fsupek <- i
        }
    }
    cancer_path <- paste0("/g/strcombio/fsupek_cancer",fsupek,"/TCGA_bam/",cancer)
    samples <- system(paste0("ls ",cancer_path), intern = TRUE)
    samples <- samples[grep("^TCGA-",samples)]

    # 2) For each sample obtain its MSIsensor score
    # First check if there are more than one MSIsensor folders
    # Means same sample with >1 vial or replicate
    # Then take the one that corresponds with the VCF somatic
    # Otherwise take directly the unique MSIsensor folder that exists
    count <- 0
    for (sample in samples) {
        folders <- system(paste0("ls ",cancer_path,"/",sample), intern = TRUE)
        MSIsensor_folders <- folders[grep("MSIsensor", folders)]
        if (length(MSIsensor_folders) > 1) {
            VCF_somatic_path <- system(paste0("readlink /g/strcombio/fsupek_cancer1/TCGA_bam/strelka_EVS/tissue_folders/",cancer,"/",sample,".vcf.gz"), intern = TRUE)
            MSIsensor_folder <- gsub("(.*Strelka2).*","\\1",VCF_somatic_path)
            MSIsensor_folder <- gsub("Strelka2","MSIsensor",MSIsensor_folder)
        } else if (length(MSIsensor_folders) == 1) {
            MSIsensor_folder <- paste0(cancer_path,"/",sample,"/",MSIsensor_folders)
        } else {
            print("Sample with no MSIsensor folder... skipping!")
            next
        }
        MSIsensor_path <- paste0(MSIsensor_folder,"/MSIsensor")
        if (file.exists(MSIsensor_path)) {
            error <- FALSE
            tryCatch({
                MSIsensor_output <- read.table(paste0(MSIsensor_folder,"/MSIsensor"), header = TRUE)
            },error = function(e) {
            error <<- TRUE
            })
        } else {
            print("Sample with no MSIsensor folder... skipping!")
            next   
        }
        if (!isTRUE(error)) {
            count <- count +1 
            df_res <- rbind(df_res,c(sample,cancer,as.numeric(MSIsensor_output[,3])))
        }
    }
    print(paste0("Cancer ",cancer," has ",count," samples with MSIsensor output"))
}

df_res <- na.omit(df_res)

# MSI status
df_res$MSI_score <- as.numeric(df_res$MSI_score)
df_res$MSI_status <- "MSS"
df_res[df_res$MSI_score >= 20,"MSI_status"] <- "MSI-H"

table(df_res$MSI_status,df_res$TCGA_cancer)

write.table(df_res, "/g/strcombio/fsupek_home/gpalou/data/TCGA_clinical/MSI_status_MSIsensor.txt", sep = "\t", quote = FALSE,
            col.names = TRUE, row.names = FALSE)



# For each sample check MSI sensor folder

# If more than one

#  readlink /g/strcombio/fsupek_cancer1/TCGA_bam/strelka_EVS/tissue_folders/COAD/TCGA-AA-3833.vcf.gz

#  Check which vial is used 
#  /g/strcombio/fsupek_cancer3/TCGA_bam/COAD/TCGA-AA-3833/d6e178d7_VS_a30ca210_Strelka2/results/variants/somatic.snvs.vcf.gz

# gsub("d6e178d7_VS_a30ca210_Strelka2","d6e178d7_VS_a30ca210_MSIsensor")
