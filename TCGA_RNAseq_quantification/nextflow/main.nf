// Parameters
parameters = ""
params.eachWithIndex { entry, i ->
	parameters = "$entry.key=$entry.value,"+"$parameters"
}
print params
dir_scripts=params.dir_scripts

// Open input file

Channel
	.fromPath(params.input_csv)
	.splitCsv(header:true)
	.map{ row -> tuple(row.file_id, row.submitter_id, row.project_id) }
	.set { samples_ch }

//samples_ch.into{samples_ch_1; samples_ch_2; samples_ch_3}
//samples_ch.subscribe { println it }
print samples_ch

process GDC_download {

	input:
	tuple val(sample_UUID), val(tcga_barcode), val(tcga_project) from samples_ch

	output:
	tuple val(sample_UUID), val(tcga_barcode), val(tcga_project), file("${sample_UUID}/*tar.gz") into fastQ_tar_gz

	"""
	#!/usr/bin/env bash

	bash ${dir_scripts}/gdc_download.sh ${task.cpus*3} ${dir_scripts}/$params.gdc_token $params.gdc_download_tries $sample_UUID

	"""
}

process fastQ_decompress {

	publishDir "/g/strcombio/fsupek_decider/gpalou/tmp/${tcga_project}/${tcga_barcode}/fastq", mode: 'copy'

	input:
	tuple val(sample_UUID), val(tcga_barcode), val(tcga_project), file(fastq_file) from fastQ_tar_gz

	output:
	tuple val(sample_UUID), val(tcga_barcode), val(tcga_project), file("*_1.fastq"), file("*_2.fastq") into fastQ_files

	"""
	#!/usr/bin/env bash

	bash ${dir_scripts}/fastq_decompress.sh ${fastq_file}

	"""
}

process STAR_alignment {

	input:
	tuple val(sample_UUID), val(tcga_barcode), val(tcga_project), file(fastq_1), file(fastq_2) from fastQ_files

	output:
	tuple val(sample_UUID), val(tcga_barcode), val(tcga_project), file("${tcga_barcode}.Aligned.toTranscriptome.out.bam") into STAR_alignment_bam
	tuple val(sample_UUID), val(tcga_barcode), val(tcga_project), file("${tcga_barcode}.Aligned.sortedByCoord.out.bam") into STAR_alignment_bam_sorted

	"""
	#!/usr/bin/env bash

	bash ${dir_scripts}/STAR_alignment.sh ${fastq_1} ${fastq_2} ${task.cpus} $tcga_barcode || exit 1

	"""

}

STAR_alignment_bam_sorted.into{STAR_alignment_bam_sorted_1; STAR_alignment_bam_sorted_2}

process RSEM_quantification {

	publishDir "${params.out_dir}/${tcga_project}/${tcga_barcode}/RSEM_quantification", mode: 'move'

	input:
	tuple val(sample_UUID), val(tcga_barcode), val(tcga_project), file(STAR_bam_file) from STAR_alignment_bam

	output:
	tuple val(sample_UUID), val(tcga_barcode), val(tcga_project), file("${tcga_barcode}.rsem.isoforms.results") into RSEM_isoforms_quantification_file
	tuple val(sample_UUID), val(tcga_barcode), val(tcga_project), file("${tcga_barcode}.rsem.genes.results") into RSEM_genes_quantification_file
	"""
	
	#!/usr/bin/env bash

	bash ${dir_scripts}/RSEM_quantification.sh ${task.cpus*3} $tcga_barcode || exit 1
	
	"""

}

