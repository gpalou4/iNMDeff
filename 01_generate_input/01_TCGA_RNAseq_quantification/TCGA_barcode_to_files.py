#!/usr/bin/python

import os
import sys

def add_barcode_to_header(file_name, header):
	""" Modify the header of a file """
	# name of temporary file
	tmp_file = file_name + '.tmp'
	# open original file in read mode and temp file in write mode
	with open(file_name, 'r') as read_obj, open(tmp_file, 'w') as write_obj:
		# Write given header to the temp file
		write_obj.write(header + '\n')
		# Read lines from original file one by one and append them to the temp file
		#next(read_obj)
		first_line = read_obj.readline()
		for line in read_obj:
			write_obj.write(line)
	# remove original file
	os.remove(file_name)
	# Rename file as the original file
	os.rename(tmp_file, file_name)

path = "/g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/primary_tumor"
TCGA_project = sys.argv[1]
level = sys.argv[2]
TCGA_samples_folders = os.listdir(path + "/"+ TCGA_project)

if level == "transcript":
	level_char = "isoforms"
elif level =="gene":
	level_char = "genes"

for TCGA_barcode in TCGA_samples_folders:
	# Obtain the path for the RSEM quantification file
	file_name_path = path + "/" + TCGA_project + "/" + TCGA_barcode + "/RSEM_quantification/" + TCGA_barcode +".rsem."+level_char+".results"
	if os.path.isfile(file_name_path):
		# New header name
		if level == "transcript":
			header = "transcript_id\tgene_id\tlength\teffective_length\t"+TCGA_barcode+"\tTPM\tFPKM\tIsoPct"
		elif level =="gene":
			header = "gene_id\ttranscript_id(s)\tlength\teffective_length\t"+TCGA_barcode+"\tTPM\tFPKM"
		# Use the file path and the barcode to add the barcode as a header in the file
		add_barcode_to_header(file_name = file_name_path, header = header)
	else:
		print(TCGA_barcode+" sample without RSEM file")
		continue

