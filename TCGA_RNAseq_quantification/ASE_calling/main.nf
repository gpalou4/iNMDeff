// Parameters
parameters = ""
params.eachWithIndex { entry, i ->
	parameters = "$entry.key=$entry.value,"+"$parameters"
}
print params
dir_scripts=params.dir_scripts

Channel
	.fromPath(params.input_csv)
	.splitCsv(header:true)
	.map{ row -> tuple(row.file_id, row.submitter_id, row.project_id, file("/g/strcombio/fsupek_decider/gpalou/tmp/${row.project_id}/${row.submitter_id}/fastq/fastq.tar.gz")) }
	.set { samples_ch }

print samples_ch

process fastQ_decompress {

	input:
	tuple val(sample_UUID), val(tcga_barcode), val(tcga_project), file(fastq) from samples_ch

	output:
	tuple val(sample_UUID), val(tcga_barcode), val(tcga_project), file("*fastq*") optional true into fastQ_files

	"""
	#!/usr/bin/env bash

	if [[ ! -f /g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/primary_tumor/${tcga_project}/${tcga_barcode}/ASE_germline/VCF_germline_strelka_recall.vcf.gz ]]
	then
		bash ${dir_scripts}/fastq_decompress.sh ${fastq}
	fi
	"""
}

process STAR_alignment {

	input:
	tuple val(sample_UUID), val(tcga_barcode), val(tcga_project), file(fastq) from fastQ_files

	output:
	tuple val(sample_UUID), val(tcga_barcode), val(tcga_project), file("${tcga_barcode}.Aligned.toTranscriptome.out.bam"), val(fastq_var) into STAR_alignment_bam
	tuple val(sample_UUID), val(tcga_barcode), val(tcga_project), file("${tcga_barcode}.Aligned.sortedByCoord.out.bam"), val(fastq_var) into STAR_alignment_bam_sorted 

	script:

	if (fastq[1] == null) {
		fastq1 = fastq.grep(~/.*.fastq/)[0]
		fastq2 = "single_end"
		fastq_var = "single_end"
	} else {
		fastq1 = fastq.grep(~/.*_1.fastq/)[0]
		fastq2 = fastq.grep(~/.*_2.fastq/)[0]
		fastq_var = "paired_end"
	}
	"""
	#!/usr/bin/env bash

	if [[ ! -f /g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/primary_tumor/${tcga_project}/${tcga_barcode}/ASE_germline/VCF_germline_strelka_recall.vcf.gz ]]
	then
		bash ${dir_scripts}/STAR_alignment.sh ${fastq1} ${fastq2} ${task.cpus} $tcga_barcode "./" || exit 1			
	fi

	"""
}

process bam_QC_control_1 {

	input:
	tuple val(sample_UUID), val(tcga_barcode), val(tcga_project), file(STAR_bam_sorted_file), val(fastq_var) from STAR_alignment_bam_sorted
 
	output:
 	tuple val(sample_UUID), val(tcga_barcode), val(tcga_project), file("find_intersecting_snps/${tcga_barcode}_q1.remap.fq*.gz"), file("find_intersecting_snps/${tcga_barcode}_q1.to.remap.bam"), file("find_intersecting_snps/${tcga_barcode}_q1.keep.bam") into bam_QC_control_1_res
	
	"""
	#!/usr/bin/env bash	
	
	if [[ ! -f /g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/primary_tumor/${tcga_project}/${tcga_barcode}/ASE_germline/VCF_germline_strelka_recall.vcf.gz ]]
	then	
		bash ${dir_scripts}/bam_QC_control_1.sh ${STAR_bam_sorted_file} $tcga_barcode $tcga_project ${task.cpus} $fastq_var
	fi
	"""

}

process STAR_alignment_remap {

	input:
	tuple val(sample_UUID), val(tcga_barcode), val(tcga_project), file(fastq), file(bam_to_remap), file(bam_to_keep) from bam_QC_control_1_res

	output:
	tuple val(sample_UUID), val(tcga_barcode), val(tcga_project), file("STAR_remap/${tcga_barcode}.Aligned.sortedByCoord.out.bam"), file(bam_to_remap), file(bam_to_keep), val(fastq_var) into STAR_alignment_bam_sorted_remap 

	script:

	if (fastq[1] == null) {
		fastq1 = fastq.grep(~/.*fq.*/)[0]
		fastq2 = "single_end"
		fastq_var = "single_end"
	} else {
		fastq1 = fastq.grep(~/.*fq1.*/)[0]
		fastq2 = fastq.grep(~/.*fq2.*/)[0]
		fastq_var = "paired_end"
	}
	"""
	#!/usr/bin/env bash

	if [[ ! -f /g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/primary_tumor/${tcga_project}/${tcga_barcode}/ASE_germline/VCF_germline_strelka_recall.vcf.gz ]]
	then
		bash ${dir_scripts}/STAR_alignment.sh ${fastq1} ${fastq2} ${task.cpus} $tcga_barcode "STAR_remap" || exit 1			
	fi

	"""
}

process bam_QC_control_2 {

	input:
	tuple val(sample_UUID), val(tcga_barcode), val(tcga_project), file(STAR_bam_sorted_file_remap), file(bam_to_remap), file(bam_to_keep), val(fastq_var) from STAR_alignment_bam_sorted_remap
 
	output:
 	tuple val(sample_UUID), val(tcga_barcode), val(tcga_project), file("${tcga_barcode}.keep.rmdup.merged.sorted.final.bam"), val(fastq_var) into bam_QC_control_2_res
	
	"""
	#!/usr/bin/env bash	
	
	if [[ ! -f /g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/primary_tumor/${tcga_project}/${tcga_barcode}/ASE_germline/VCF_germline_strelka_recall.vcf.gz ]]
	then	
		bash ${dir_scripts}/bam_QC_control_2.sh ${STAR_bam_sorted_file_remap} ${bam_to_remap} ${bam_to_keep} $tcga_barcode $tcga_project ${task.cpus} $fastq_var
	fi
	"""

}

bam_QC_control_2_res.into{bam_QC_control_2_res_1; bam_QC_control_2_res_2}

process strelka_germline {

	publishDir "${params.out_dir}/${tcga_project}/${tcga_barcode}/ASE_germline", mode: 'copy'

	input:
	tuple val(sample_UUID), val(tcga_barcode), val(tcga_project), file(STAR_bam_sorted_file_fixed), val(fastq_var) from bam_QC_control_2_res_1
 
	output:
 	tuple val(sample_UUID), val(tcga_barcode), val(tcga_project), file("VCF_germline_strelka_recall.vcf.gz") optional true into strelka_germline
	
	"""
	#!/usr/bin/env bash	
	
	if [[ ! -f /g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/primary_tumor/${tcga_project}/${tcga_barcode}/ASE_germline/VCF_germline_strelka_recall.vcf.gz ]]
	then	
		bash ${dir_scripts}/strelka2_germline.sh ${STAR_bam_sorted_file_fixed} $tcga_barcode $tcga_project ${task.cpus}
	fi
	"""

}

process strelka_somatic {

	publishDir "${params.out_dir}/${tcga_project}/${tcga_barcode}/ASE_somatic", mode: 'copy'

	input:
	tuple val(sample_UUID), val(tcga_barcode), val(tcga_project), file(STAR_bam_sorted_file_fixed), val(fastq_var) from bam_QC_control_2_res_2

	output:
	tuple val(sample_UUID), val(tcga_barcode), val(tcga_project), file("VCF_somatic_snvs_strelka_recall.vcf.gz"), file("VCF_somatic_indels_strelka_recall.vcf.gz") optional true into strelka_somatic

	"""
	#!/usr/bin/env bash
	
	if [[ ! -f /g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/primary_tumor/${tcga_project}/${tcga_barcode}/ASE_somatic/VCF_somatic_indels_strelka_recall.vcf.gz ]]
	then
		bash ${dir_scripts}/strelka2_somatic.sh ${STAR_bam_sorted_file_fixed} $tcga_barcode $tcga_project ${task.cpus}
	fi
	"""

}

