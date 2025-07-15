library("vcfR")
library("ggplot2")
library("reshape2")

# args <- commandArgs(trailingOnly=TRUE)
# TCGA_cancer <- args[1]
# TCGA_sample <- args[2]

# Input data
sample_batches_ASE <- read.csv("/g/strcombio/fsupek_home/gpalou/data/TCGA_RNAseq_quantification/sample_batches_ASE/prueba.csv")

df_res <- data.frame(samples = sample_batches_ASE$submitter_id)

for (i in 1:nrow(sample_batches_ASE)) {

    TCGA_cancer <- gsub("TCGA-","",as.character(sample_batches_ASE[i,"project_id"]))
    TCGA_sample <- as.character(sample_batches_ASE[i,"submitter_id"])

    # 1) VCFs Germline
    input.path <- system(paste0("grep ", TCGA_sample," /g/strcombio/fsupek_home/gpalou/data/TCGA_germline_variants/",TCGA_cancer,"_VCFs_samples_list.txt"), intern = TRUE)
    input <- read.vcfR(input.path, verbose = FALSE)
    # Strelka output
    output_strelka <- read.vcfR(paste0("~/data/TCGA_RNAseq_quantification/primary_tumor/TCGA-",TCGA_cancer,"/",TCGA_sample,"/ASE_germline/VCF_germline_strelka_recall.vcf.gz"), verbose = FALSE )

    input_df <- as.data.frame(input@fix)
    output_strelka_df <- as.data.frame(output_strelka@fix)
    a <- merge(input_df,output_strelka_df,by.x=c("POS","REF","ALT"), by.y=c("POS","REF","ALT"))
    # SNVs
    a.snvs <- a[which(nchar(as.character(a$REF)) == 1 & nchar(as.character(a$ALT)) == 1),]
    input_df_snvs <- input_df[which(nchar(as.character(input_df$REF)) == 1 & nchar(as.character(input_df$ALT)) == 1),]
    output_strelka_df_snvs <- output_strelka_df[which(nchar(as.character(output_strelka_df$REF)) == 1 & nchar(as.character(output_strelka_df$ALT)) == 1),]
    # Indels
    a.indels <- a[which(nchar(as.character(a$REF)) >1 | nchar(as.character(a$ALT)) >1),]
    input_df_indels <- input_df[which(nchar(as.character(input_df$REF)) >1 | nchar(as.character(input_df$ALT)) >1),]
    output_strelka_df_indels <- output_strelka_df[which(nchar(as.character(output_strelka_df$REF)) > 1 | nchar(as.character(output_strelka_df$ALT)) > 1),]

    df_res[which(df_res$samples%in%TCGA_sample),"VCF_germline_input"] <- nrow(input_df)
    df_res[which(df_res$samples%in%TCGA_sample),"strelka_germline_output"] <- nrow(output_strelka_df)
    df_res[which(df_res$samples%in%TCGA_sample),"VCF_germline_snvs"] <- nrow(input_df_snvs)
    df_res[which(df_res$samples%in%TCGA_sample),"strelka_germline_snvs"] <- nrow(output_strelka_df_snvs)
    df_res[which(df_res$samples%in%TCGA_sample),"VCF_strelka_germline_snvs_shared"] <- nrow(a.snvs)
    df_res[which(df_res$samples%in%TCGA_sample),"VCF_germline_indels"] <- nrow(input_df_indels)
    df_res[which(df_res$samples%in%TCGA_sample),"strelka_germline_indels"] <- nrow(output_strelka_df_indels)
    df_res[which(df_res$samples%in%TCGA_sample),"VCF_strelka_germline_indels_shared"] <- nrow(a.indels)

    # Allelecounter (Samtools mpileup) output
    output_samtools <- read.table(paste0("/g/strcombio/fsupek_home/gpalou/data/TCGA_RNAseq_quantification/primary_tumor/TCGA-",TCGA_cancer,"/",TCGA_sample,"/ASE_germline/ase_germline_counts.txt"), 
                        header = TRUE, row.names = NULL)

    b <- merge(input_df,output_samtools,by.x=c("POS","REF","ALT"), by.y=c("position","refAllele","altAllele"))

    # SNVs
    b.snvs <- b[which(nchar(as.character(b$REF)) == 1 & nchar(as.character(b$ALT)) == 1),]
    output_samtools_df_snvs <- output_samtools[which(nchar(as.character(output_samtools$refAllele)) == 1 & nchar(as.character(output_samtools$altAllele)) == 1),]

    # Indels
    b.indels <- b[which(nchar(as.character(b$REF)) >1 | nchar(as.character(b$ALT)) >1),]
    output_samtools_df_indels <- output_samtools[which(nchar(as.character(output_samtools$refAllele)) > 1 | nchar(as.character(output_samtools$altAllele)) > 1),]

    df_res[which(df_res$samples%in%TCGA_sample),"samtools_germline_output"] <- nrow(output_samtools)
    df_res[which(df_res$samples%in%TCGA_sample),"samtools_germline_snvs"] <- nrow(output_samtools_df_snvs)
    df_res[which(df_res$samples%in%TCGA_sample),"VCF_samtools_germline_snvs_shared"] <- nrow(b.snvs)
    df_res[which(df_res$samples%in%TCGA_sample),"samtools_germline_indels"] <- nrow(output_samtools_df_indels)
    df_res[which(df_res$samples%in%TCGA_sample),"VCF_samtools_germline_indels_shared"] <- nrow(b.indels)

    # Shared Samtools and strelka
    strelka_samtools_snvs <- merge(output_strelka_df_snvs,output_samtools_df_snvs,by.x=c("POS","REF","ALT"), by.y=c("position","refAllele","altAllele"))
    strelka_samtools_indels <- merge(output_strelka_df_indels,output_samtools_df_indels,by.x=c("POS","REF","ALT"), by.y=c("position","refAllele","altAllele"))

    df_res[which(df_res$samples%in%TCGA_sample),"samtools_strelka_germline_snvs_shared"] <- nrow(strelka_samtools_snvs)
    df_res[which(df_res$samples%in%TCGA_sample),"samtools_strelka_germline_indels_shared"] <- nrow(strelka_samtools_indels)

    # 2) VCFs somatic SNVs Strelka output
    input.path <- system(paste0("readlink -f /g/strcombio/fsupek_cancer1/TCGA_bam/strelka_EVS/tissue_folders/",TCGA_cancer,"/",TCGA_sample,".vcf.gz"), intern = TRUE)
    input <- read.vcfR(input.path, verbose = FALSE)

    output_strelka <- read.vcfR(paste0("~/data/TCGA_RNAseq_quantification/primary_tumor/TCGA-",TCGA_cancer,"/",TCGA_sample,"/ASE_somatic/VCF_somatic_snvs_strelka_recall.vcf.gz"), verbose = FALSE )

    input_df <- as.data.frame(input@fix)
    output_strelka_df <- as.data.frame(output_strelka@fix)
    a <- merge(input_df,output_strelka_df,by.x=c("POS","REF","ALT"), by.y=c("POS","REF","ALT"))

    # SNVs
    a.snvs <- a[which(nchar(as.character(a$REF)) == 1 & nchar(as.character(a$ALT)) == 1),]
    input_df_snvs <- input_df[which(nchar(as.character(input_df$REF)) == 1 & nchar(as.character(input_df$ALT)) == 1),]
    output_strelka_df_snvs <- output_strelka_df[which(nchar(as.character(output_strelka_df$REF)) == 1 & nchar(as.character(output_strelka_df$ALT)) == 1),]

    df_res[which(df_res$samples%in%TCGA_sample),"VCF_somatic_snvs"] <- nrow(input_df_snvs)
    df_res[which(df_res$samples%in%TCGA_sample),"strelka_somatic_snvs"] <- nrow(output_strelka_df_snvs)
    df_res[which(df_res$samples%in%TCGA_sample),"VCF_strelka_somatic_snvs_shared"] <- nrow(a.snvs)

    # 3) VCFs somatic indels Strelka output

    input.path <- gsub("somatic.snvs.vcf.gz","somatic.indels.vcf.gz",input.path)
    input <- read.vcfR(input.path, verbose = FALSE)

    # Strelka output
    output_strelka <- read.vcfR(paste0("~/data/TCGA_RNAseq_quantification/primary_tumor/TCGA-",TCGA_cancer,"/",TCGA_sample,"/ASE_somatic/VCF_somatic_indels_strelka_recall.vcf.gz"), verbose = FALSE )

    input_df <- as.data.frame(input@fix)
    output_strelka_df <- as.data.frame(output_strelka@fix)
    b <- merge(input_df,output_strelka_df,by.x=c("POS","REF","ALT"), by.y=c("POS","REF","ALT"))
    
    # Indels
    b.indels <- b[which(nchar(as.character(b$REF)) >1 | nchar(as.character(b$ALT)) >1),]
    input_df_indels <- input_df[which(nchar(as.character(input_df$REF)) >1 | nchar(as.character(input_df$ALT)) >1),]
    output_strelka_df_indels <- output_strelka_df[which(nchar(as.character(output_strelka_df$REF)) > 1 | nchar(as.character(output_strelka_df$ALT)) > 1),]

    df_res[which(df_res$samples%in%TCGA_sample),"VCF_somatic_indels"] <- nrow(input_df_indels)
    df_res[which(df_res$samples%in%TCGA_sample),"strelka_somatic_indels"] <- nrow(output_strelka_df_indels)
    df_res[which(df_res$samples%in%TCGA_sample),"VCF_strelka_somatic_indels_shared"] <- nrow(b.indels)

    # 4) VCFs somatic SNVs Allelecounter (Samtools mpileup) output

    output_samtools <- read.table(paste0("/g/strcombio/fsupek_home/gpalou/data/TCGA_RNAseq_quantification/primary_tumor/TCGA-",TCGA_cancer,"/",TCGA_sample,"/ASE_somatic/ase_somatic_snvs_counts.txt"), 
                        header = TRUE, row.names = NULL)

    a <- merge(input_df_snvs,output_samtools,by.x=c("POS","REF","ALT"), by.y=c("position","refAllele","altAllele"))

    # SNVs
    a.snvs <- a[which(nchar(as.character(a$REF)) == 1 & nchar(as.character(a$ALT)) == 1),]
    output_samtools_df_snvs <- output_samtools[which(nchar(as.character(output_samtools$refAllele)) == 1 & nchar(as.character(output_samtools$altAllele)) == 1),]

    df_res[which(df_res$samples%in%TCGA_sample),"samtools_somatic_snvs"] <- nrow(output_samtools_df_snvs)
    df_res[which(df_res$samples%in%TCGA_sample),"VCF_samtools_somatic_snvs_shared"] <- nrow(a.snvs)

    # 5) VCFs somatic indels Allelecounter (Samtools mpileup) output

    output_samtools <- read.table(paste0("/g/strcombio/fsupek_home/gpalou/data/TCGA_RNAseq_quantification/primary_tumor/TCGA-",TCGA_cancer,"/",TCGA_sample,"/ASE_somatic/ase_somatic_indels_counts.txt"), 
                        header = TRUE, row.names = NULL)

    b <- merge(input_df_indels,output_samtools,by.x=c("POS","REF","ALT"), by.y=c("position","refAllele","altAllele"))

    # Indels
    b.indels <- b[which(nchar(as.character(b$REF)) > 1 | nchar(as.character(b$ALT)) > 1),]
    output_samtools_df_indels <- output_samtools[which(nchar(as.character(output_samtools$refAllele)) > 1 | nchar(as.character(output_samtools$altAllele)) > 1),]

    df_res[which(df_res$samples%in%TCGA_sample),"samtools_somatic_indels"] <- nrow(output_samtools_df_indels)
    df_res[which(df_res$samples%in%TCGA_sample),"VCF_samtools_somatic_indels_shared"] <- nrow(b.indels)

    # Shared Samtools and strelka
    strelka_samtools_snvs <- merge(output_strelka_df_snvs,output_samtools_df_snvs,by.x=c("POS","REF","ALT"), by.y=c("position","refAllele","altAllele"))
    strelka_samtools_indels <- merge(output_strelka_df_indels,output_samtools_df_indels,by.x=c("POS","REF","ALT"), by.y=c("position","refAllele","altAllele"))

    df_res[which(df_res$samples%in%TCGA_sample),"samtools_strelka_somatic_snvs_shared"] <- nrow(strelka_samtools_snvs)
    df_res[which(df_res$samples%in%TCGA_sample),"samtools_strelka_somatic_indels_shared"] <- nrow(strelka_samtools_indels)
}

write.table(df_res, file = "/home/gpalou/projects/NMD/scripts/TCGA_RNAseq_quantification/ASE_calling/check_input_output_df.txt", quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)

# Barplot
# df_res_rescale <- apply(df_res[,-grep("samples|shared",colnames(df_res))],2,function(col) {
#                             rescale(as.numeric(col))
#                     })
# df_res_rescale <- cbind(df_res_rescale,df_res[,grep("samples|shared",colnames(df_res))])


for (type in c("germline","somatic")) {

    for (mut in c("snvs","indels")) {

        #df_res_filt <- df_res[,c("samples",paste0("VCF_",type,"_",mut),paste0(tool,"_",type,"_",mut),paste0("VCF_",tool,"_",type,"_",mut,"_shared"))]

        df_res_filt <- df_res[,c("samples",paste0("VCF_",type,"_",mut),paste0("strelka_",type,"_",mut),paste0("VCF_strelka_",type,"_",mut,"_shared"),
                            paste0("samtools_",type,"_",mut),paste0("VCF_samtools_",type,"_",mut,"_shared"), paste0("samtools_strelka_",type,"_",mut,"_shared"))]
        df_res_melt <- melt(df_res_filt)

        png(file = paste0("/home/gpalou/projects/NMD/scripts/TCGA_RNAseq_quantification/ASE_calling/checks/",type,"_",mut,".png"), width = 4000, height = 3000, res = 300)
        p <- ggplot(data=df_res_melt, aes(x=value, y=samples, fill = variable)) +
            geom_bar(stat="identity", position=position_dodge()) + ggtitle(paste0(type,"_",mut))
        print(p)
        dev.off()

    }

}

df_res_filt <- df_res[,c("samples","VCF_germline_snvs","strelka_germline_snvs","VCF_strelka_germline_snvs_shared")]
df_res_melt <- melt(df_res_filt)




