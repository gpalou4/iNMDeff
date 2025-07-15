# conda activate /home/gpalou/.conda/envs/biomaRt
library("biomaRt")

listEnsemblArchives()
ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = "apr2018.archive.ensembl.org", dataset = "hsapiens_gene_ensembl")
#ensembl = useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", version = 92)
#listEnsembl(version = 92)

attributes = listAttributes(ensembl, page = "structure")
attributes[grep("transcript", attributes$description, ignore.case = TRUE), ]

transcript_start_end_coords <- getBM(attributes = c("transcription_start_site","chromosome_name",
                                                    "transcript_start", "transcript_end",
                                                    "strand",  "ensembl_gene_id",
                                                    "ensembl_transcript_id", "external_gene_name"),
                                    mart = ensembl)
head(transcript_start_end_coords)

# Create Transcription End Site (TES)
transcript_start_end_coords$transcription_end_site <- NA
transcript_start_end_coords[which(transcript_start_end_coords$strand %in% 1),"transcription_end_site"] <- transcript_start_end_coords[which(transcript_start_end_coords$strand %in% 1),"transcript_end"]
transcript_start_end_coords[which(transcript_start_end_coords$strand %in% -1),"transcription_end_site"] <- transcript_start_end_coords[which(transcript_start_end_coords$strand %in% -1),"transcript_start"]

output_path <- "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/ensembl_v92_transcript_start_end_coords.txt"
write.table(transcript_start_end_coords, file = output_path, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
