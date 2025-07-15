rm(list=ls())

# 1) Load Data

# Arguments and paths

args <- commandArgs(trailingOnly=TRUE)
paths <- read.table(file = args[1], header = TRUE, sep = ",", stringsAsFactors = FALSE)
GTEx.names.path <- paths[paths$folder_or_object=="GTEx_names_path","path_or_filename"]
GTEx.NB.res.tissue.path <- paths[paths$folder_or_object=="GTEx_NB_res_tissue_path","path_or_filename"]
GTEx.NB.res.list.path <- paths[paths$folder_or_object=="GTEx_NB_res_list_path","path_or_filename"]

# 1.2) GTEx NB results

GTEx.tissues <- read.table(file = paste0(GTEx.names.path,paths[paths$folder_or_object=="GTEx_names","path_or_filename"]))$V1
# Acronyms of tissues
acronyms <- c("ADPSBQ","ADPVSC","ADRNLG","ARTAORT","ARTCRN","ARTTBL","BLADDER","BRNAMY","BRNACC","BRNCDT","BRNCHB","BRNCHA","BRNCTXA","BRNCTXB",
              "BRNHPP","BRNHPT","BRNNCC","BRNPTM","BRNSCP","BRNSNG","BREAST","FIBRBLS","LCL","CML","CVXECT","CVXEND","CLNSGM","CLNTRN","ESPGEJ","ESPMCS","ESPMSL",
              "FLLPNT","HRTAA","HRTLV","KDNCTX","KDNMDL","LIVER","LUNG","SLVRYG","MSCLSK","NERVET","OVARY","PNCREAS","PTTARY","PRSTTE","SKINNS","SKINS",
              "SNTTRM","SPLEEN","STMACH","TESTIS","THYROID","UTERUS","VAGINA","WHLBLD")
tissues.acronyms <- data.frame(tissues = GTEx.tissues, acronyms = acronyms )

GTEx.NB.res.list <- list()
for (i in seq(1:length(GTEx.tissues))) {
  GTEx.tissue <- as.character(GTEx.tissues[i])
  acronym <- as.character(tissues.acronyms[tissues.acronyms$tissues==GTEx.tissue,"acronyms"])
  print(GTEx.tissue)
  print(acronym)
  nb.error <- FALSE
  tryCatch( {   NB.res <- read.table(file = gsub("\\[X\\]", GTEx.tissue, paste0(GTEx.NB.res.tissue.path,paths[paths$folder_or_object=="GTEx_NB_res_tissue","path_or_filename"])),
                                     header = TRUE, sep = "\t", row.names = 1)
  GTEx.NB.res.list[[acronym]] <- NB.res }
  ,error = function(e) {nb.error <<- TRUE})
  if (isTRUE(nb.error)) {next}
}

save(GTEx.NB.res.list, file=paste0(GTEx.NB.res.list.path,paths[paths$folder_or_object=="GTEx_NB_res_list","path_or_filename"]))

# Save for Nextflow
save(GTEx.NB.res.list, file=paths[paths$folder_or_object=="GTEx_NB_res_list","path_or_filename"])
