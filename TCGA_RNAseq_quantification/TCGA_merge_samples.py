#!/usr/bin/python3

import os
import sys
import pandas

def merge_samples_cancer(TCGA_samples_folders, type, level):
	first = True
	sample_barcodes = {}
	# Raw values (expected_counts)
	if type == "raw":
		column_number = 4
	# TPM values
	elif type == "TPM":
		column_number = 5
	if level == "transcript":
		level_char = "isoforms"
		id_col_char = "transcript_id"
	elif level =="gene":
		level_char = "genes"
		id_col_char = "gene_id"
	for TCGA_barcode in TCGA_samples_folders:
		# Obtain the path for the RSEM quantification file
		file_name_path = path + "/" + TCGA_project + "/" + TCGA_barcode + "/RSEM_quantification/" + TCGA_barcode +".rsem."+level_char+".results"
		if os.path.isfile(file_name_path):
			# Open the file
			file = pandas.read_csv(filepath_or_buffer = file_name_path, sep = "\t")
			# Dictionary with TCGA sample barcode
			if TCGA_barcode in sample_barcodes:
				continue
			else:
				# Keep Transcripts and Raw/TMP columns
				file = file.iloc[:,[0,column_number]]
				# Change the column name to TCGA barcode
				if type == "TPM":
					file.columns.values[1] = TCGA_barcode
				# Merge the file to the merged matrix
				# Check if there is a matrix or we have to create one (first iteration)
				if first:
					matrix = file
					first = False
				else:
					matrix = pandas.merge(matrix, file, on = id_col_char, how = "left")
				sample_barcodes[TCGA_barcode] = 1
	return matrix

path = "/g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/primary_tumor"
TCGA_project = sys.argv[1]
level = sys.argv[2]
if level == "transcript":
	id_col_char = "transcript_id"
elif level =="gene":
	id_col_char = "gene_id"
TCGA_samples_folders = os.listdir(path + "/"+ TCGA_project)

# 1) RNA-seq matrix raw values

RNAseq_matrix = merge_samples_cancer(TCGA_samples_folders, type = "raw", level = level)
# Remove genes with PAR_Y
RNAseq_matrix_filt = RNAseq_matrix[~RNAseq_matrix[id_col_char].str.contains(r'_PAR_Y')]
# Remove ID version
RNAseq_matrix_filt[id_col_char] = RNAseq_matrix_filt[id_col_char].str.replace(r'\..*','')
RNAseq_matrix_filt.to_csv(path + "/" + TCGA_project + "/" + TCGA_project + "_RNAseq_matrix_raw_"+level+".txt", index=False, sep='\t')

# 2) RNA-seq matrix TPM values

RNAseq_matrix = merge_samples_cancer(TCGA_samples_folders, type = "TPM", level = level)
# Remove genes with PAR_Y
RNAseq_matrix_filt = RNAseq_matrix[~RNAseq_matrix[id_col_char].str.contains(r'_PAR_Y')]
# Remove ID version
RNAseq_matrix_filt[id_col_char] = RNAseq_matrix_filt[id_col_char].str.replace(r'\..*','')

RNAseq_matrix_filt.to_csv(path + "/" + TCGA_project + "/" + TCGA_project + "_RNAseq_matrix_TPM_"+level+".txt", index=False, sep='\t')
