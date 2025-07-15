
# 1) Arguments and data
CGC <- read.csv(file = "/g/strcombio/fsupek_cancer1/gpalou/COSMIC/cancer_gene_census_updated.tsv", 
                    header = TRUE, stringsAsFactors = FALSE)
CGC <- CGC[,c("Gene.Symbol","Role.in.Cancer")]

CGC_genes <- CGC[which(CGC$Role.in.Cancer != "random_control_gene"),"Gene.Symbol"]
random_genes <- CGC[which(CGC$Role.in.Cancer == "random_control_gene"),"Gene.Symbol"]

alpha <- "1e-04"
robust_SPCA <- "no"
CNV_PCs <- 1000

for (CNV_PCs in c(300,350,400,450,500)) {

    print(paste0("PC ---> ",CNV_PCs))
    # 1) Open association results
    output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/somatic_associations/optimizing_PCs/associations/pancancer_endogenous_CGC_somatic_mut_CNV_PCs_",CNV_PCs,"_test_1_",alpha,"_robust_",robust_SPCA,".txt") 
    df_ass <- read.table(file = output_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    
    # 2) Calculate Lambda for CGC genes
    df_ass_filt <- df_ass[df_ass$Gene_symbol %in% CGC_genes,]
    df_lambda_CGC_genes <- data.frame(CNV_amp = NA,CNV_del = NA, mut_missense = NA, mut_truncating = NA , mut_synonymous = NA)
    for (mutation in c("CNV_amp","CNV_del","mut_missense","mut_truncating","mut_synonymous")) {
        p_values <- df_ass_filt[,paste0("som_",mutation,"_pval")]
        # Lambda calculation
        chisq <- qchisq(1-p_values,1)
        lambda = round(median(chisq,na.rm=TRUE)/qchisq(0.5,1),2)
        df_lambda_CGC_genes[,mutation] <- lambda
    }
    print(df_lambda_CGC_genes)
    # 3) Calculate Lambda for Random genes
    df_ass_filt <- df_ass[df_ass$Gene_symbol %in% random_genes,]
    df_lambda_random_genes <- data.frame(CNV_amp = NA,CNV_del = NA, mut_missense = NA, mut_truncating = NA , mut_synonymous = NA)
    for (mutation in c("CNV_amp","CNV_del","mut_missense","mut_truncating","mut_synonymous")) {
        p_values <- df_ass_filt[,paste0("som_",mutation,"_pval")]
        # Lambda calculation
        chisq <- qchisq(1-p_values,1)
        lambda = round(median(chisq,na.rm=TRUE)/qchisq(0.5,1),2)
        df_lambda_random_genes[,mutation] <- lambda
    }
    print(df_lambda_random_genes)


}


