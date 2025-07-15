glm_nb_regression <- function(NMD.final.df, genes.size = 50, offset1 = NULL, offset2 = NULL, offset3 = NULL, PCs.names = NULL) {

  # Subsample if dataframe has more than 100 rows
  if (nrow(NMD.final.df) > (genes.size*2) ) {
    random.genes <- sample(unique(NMD.final.df$ensembl_gene_id),genes.size)
    NMD.final.df.filt <- NMD.final.df[NMD.final.df$ensembl_gene_id %in% random.genes,]
    NMD.final.df.filt <- NMD.final.df.filt[order(NMD.final.df.filt$ensembl_gene_id),]
  } else {
    NMD.final.df.filt <- NMD.final.df
  }
  # NMD target variable converted to numeric   
  NMD.final.df.filt$NMD_target <- ifelse(NMD.final.df.filt$final_consensus == "NMD_target", 1, 0 )
  df <- NMD.final.df.filt
  print(dim( NMD.final.df.filt))
  # Bayesian NB regression model
  glm.nb.model <- "stan_glm(raw_gene_exp ~ as.factor(NMD_target) + as.factor(ensembl_gene_id) + length"
  # Raw model
  if (is.null(PCs.names) & is.null(offset1) & is.null(offset2)) {glm.nb.model = paste(glm.nb.model, ", data = df, family = neg_binomial_2, verbose = FALSE)") }
  # Add Principal Components
  if (!is.null(PCs.names)) {for (PC in PCs.names) glm.nb.model = paste(glm.nb.model, "+", PC) }
  # Add offset 1 (transcript length)
  if (!is.null(offset1) & is.null(offset2)) {glm.nb.model = paste(glm.nb.model, ",", "offset = log(df$offset1), data = df, family = neg_binomial_2, verbose = FALSE)") }
  # Add offset 2 (sample library size)
  if (!is.null(offset1) & !is.null(offset2) & is.null(offset3)) {glm.nb.model = paste(glm.nb.model, ",", "offset = log(df$offset1*df$offset2), data = df, family = neg_binomial_2, verbose = FALSE)") }
  # Add offset 3 (sample tumor purity)
  if (!is.null(offset1) & !is.null(offset2) & !is.null(offset3)) {glm.nb.model = paste(glm.nb.model, ",", "offset = log(df$offset1*df$offset2*df$offset3), data = df, family = neg_binomial_2, verbose = FALSE)") }

  # Perform the regression
  tryCatch( { 
    nb.res <- NA
    nb.res <- eval(parse(text = glm.nb.model))
  }, error = function(e) {
    print("################################################################################################")
    print("ERROR in NB regression. Model:")
    print(e)
    nb.error.3 <<- TRUE
    print(glm.nb.model)
    print(head(df))
    print(paste0("purity: ",offset3))
    print("################################################################################################")
  })
  return(nb.res)
}
.libPaths( rev( .libPaths() ) )
# Bayesian for NB regression
library("V8")
# library("shinystan")
library("rstan")
library("rstanarm")

args <- commandArgs(trailingOnly=TRUE)
paths <- read.table(file = args[1], header = TRUE, sep = ",", stringsAsFactors = FALSE)
GTEx.tissue <- args[2]
GTEx.sample <- args[3]
NMD.gene.set.name <- args[4]

print(GTEx.tissue)
print(GTEx.sample)
print(NMD.gene.set.name)

NMD.endogenous.transcripts.path <- paths[paths$folder_or_object=="NMD_endogenous_transcripts_path","path_or_filename"]
res.path <- gsub("\\[X1\\]",GTEx.tissue, paste0(NMD.endogenous.transcripts.path))
res.path <- gsub("\\[X2\\]",gsub("\\.","-",GTEx.sample), res.path)
# 1) Open NMD final table
NMD.transcripts.res.path <- paste0(res.path,paths[paths$folder_or_object=="NMD_endogenous_transcripts","path_or_filename"])
NMD.transcripts.res.path <- gsub("\\[X\\]",gsub("\\.","_",NMD.gene.set.name), NMD.transcripts.res.path)
NMD.final.df <- read.table(file = NMD.transcripts.res.path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# 2) Open NB reg results table
NB.reg.res.path <- paste0(res.path,paths[paths$folder_or_object=="NB_regression_results","path_or_filename"])
NB.reg.res.path <- gsub("\\[X\\]",gsub("\\.","_",NMD.gene.set.name), NB.reg.res.path)
nb.coeff.res <- read.table(file = NB.reg.res.path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
genes.size <- nb.coeff.res$num_NMD_targets
# 3) NB regression
nb.res <- glm_nb_regression(NMD.final.df = NMD.final.df, genes.size = genes.size,
                            offset1 = NULL, offset2 = NULL, offset3 = NULL, PCs.names = NULL)
# 4) Save output
if (is.na(nb.res)) {
    nb.res.pval <- NA
    nb.res.coeff <- NA
    #if (NMD.gene.set.name == "NMD.global.2.shared") {
    nb.coeff.res[1,"error"] <- "NB regression error" 
    #}
    nb.coeff.res[1,NMD.gene.set.name] <- nb.res.coeff
    nb.coeff.res[1,paste0(NMD.gene.set.name,".sd")] <- NA
    nb.coeff.res[1,paste0(NMD.gene.set.name,".CI_2.5")] <- NA
    nb.coeff.res[1,paste0(NMD.gene.set.name,".CI_97.5")] <- NA
} else if (!is.na(nb.res)) {
    nb.res.coeff <- as.numeric(coef(nb.res)[2])
    # Store NB coefficient
    if (NMD.gene.set.name == "NMD.global.2.shared") {
        print(paste0("NMD efficiency --> ", nb.res.coeff))
    }
    nb.coeff.res[1,NMD.gene.set.name] <- nb.res.coeff
    nb.coeff.res[1,paste0(NMD.gene.set.name,".sd")] <- as.numeric(nb.res$stan_summary[2,"sd"])
    #nb.coeff.res[1,paste0(NMD.gene.set.name,".pval")] <- nb.res.coeff
    nb.coeff.res[1,paste0(NMD.gene.set.name,".CI_2.5")] <- as.numeric(nb.res$stan_summary[2,"2.5%"])
    nb.coeff.res[1,paste0(NMD.gene.set.name,".CI_97.5")] <- as.numeric(nb.res$stan_summary[2,"97.5%"])
}
write.table(nb.coeff.res, file = NB.reg.res.path, sep = "\t", quote = FALSE,
        col.names = TRUE, row.names = FALSE)












