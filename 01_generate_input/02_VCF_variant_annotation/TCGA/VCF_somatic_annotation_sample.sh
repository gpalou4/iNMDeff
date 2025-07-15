#!/usr/bin/bash

#SBATCH --partition=normal_prio
##SBATCH --account=high_prio
#SBATCH --job-name=annovar_sample
#SBATCH --output=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.out
#SBATCH --error=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.err
#SBATCH --ntasks=1
##SBATCH --cpus-per-task=X
##SBATCH --nodes=X
##SBATCH --mem-per-cpu=X
#SBATCH --mem=16GB
#SBATCH --time=8:00:00
#SBATCH --mail-type=ALL
##SBATCH --mail-user=guillermo.palou@irbbarcelona.org

# Activate base environment
conda activate base

TCGA_project=$1
VCF_somatic_snvs_sample=$2
VCF_somatic_indels_sample=$3
VCF_somatic_annotation_path=$4
TCGA_barcode=$5

DIR_SCRIPTS="/g/strcombio/fsupek_home/gpalou/projects/NMD/scripts/TCGA_RNAseq_quantification"
tmp_dir="/scratch/2/gpalou/trash"

# Check if sample already exists

# SNVs
if [ ! -f "${VCF_somatic_annotation_path}/${TCGA_barcode}_snvs_annotated.HG_multianno.vcf.gz" ]; then
	
	# Fix Genotype, GT, tag (it doesn't exist in somatic Strelka2 VCF output)
	# We will create the tag with a modified script found in bcbio-nextgen
	gunzip -c $VCF_somatic_snvs_sample > ${tmp_dir}/${TCGA_barcode}_snvs_somatic.vcf
	python3 ${DIR_SCRIPTS}/strelka2_somatic_tumor_normal_genotypes.py ${tmp_dir}/${TCGA_barcode}_snvs_somatic.vcf
	# Annotate sample VCFs with annovar (hg38 + gnomAD and exac Allele Frequencies)
	# Beware the refGene is actually GENCODE/ENSEMBL, I built the custom gene annotation as stated here --> /g/strcombio/fsupek_cancer1/gpalou/human_genome/humandb_annovar_hg38/commands.txt
	~/software/annovar/table_annovar.pl -buildver HG -out ${VCF_somatic_annotation_path}/${TCGA_barcode}_snvs_annotated -remove -protocol refGene,gnomad211_exome,exac03 -operation g,f,f -nastring . -polish -vcfinput ${tmp_dir}/${TCGA_barcode}_snvs_somatic_GT_fixed.vcf /g/strcombio/fsupek_cancer1/gpalou/human_genome/humandb_annovar_hg38
	# Zip and Fix the output files
	sed 's/\\x3b/-/g'  ${VCF_somatic_annotation_path}/${TCGA_barcode}_snvs_annotated.HG_multianno.vcf | sed 's/\\x3d/:/g' | gzip -c > ${VCF_somatic_annotation_path}/${TCGA_barcode}_snvs_annotated.HG_multianno.vcf.gz
	sed 's/\\x3b/-/g'  ${VCF_somatic_annotation_path}/${TCGA_barcode}_snvs_annotated.HG_multianno.txt | sed 's/\\x3d/:/g' | gzip -c > ${VCF_somatic_annotation_path}/${TCGA_barcode}_snvs_annotated.HG_multianno.txt.gz
	# Remove temp files
	rm ${VCF_somatic_annotation_path}/${TCGA_barcode}_snvs_annotated.HG_multianno.vcf ${VCF_somatic_annotation_path}/${TCGA_barcode}_snvs_annotated.HG_multianno.txt ${tmp_dir}/${TCGA_barcode}_snvs_somatic.vcf ${tmp_dir}/${TCGA_barcode}_snvs_somatic_GT_fixed.vcf 
	gzip ${VCF_somatic_annotation_path}/${TCGA_barcode}_snvs_annotated.avinput
fi

# Indels
if [ ! -f "${VCF_somatic_annotation_path}/${TCGA_barcode}_indels_annotated.HG_multianno.vcf.gz" ]; then
	# Fix Genotype, GT, tag (it doesn't exist in somatic Strelka2 VCF output)
	# We will create the tag with a modified script found in bcbio-nextgen
	gunzip -c $VCF_somatic_indels_sample > ${tmp_dir}/${TCGA_barcode}_indels_somatic.vcf
	python3 ${DIR_SCRIPTS}/strelka2_somatic_tumor_normal_genotypes.py ${tmp_dir}/${TCGA_barcode}_indels_somatic.vcf
	# Annotate sample VCFs with annovar (hg38 + gnomAD and exac Allele Frequencies)
	~/software/annovar/table_annovar.pl -buildver HG -out ${VCF_somatic_annotation_path}/${TCGA_barcode}_indels_annotated -remove -protocol refGene,gnomad211_exome,exac03 -operation g,f,f -nastring . -polish -vcfinput ${tmp_dir}/${TCGA_barcode}_indels_somatic_GT_fixed.vcf /g/strcombio/fsupek_cancer1/gpalou/human_genome/humandb_annovar_hg38
	# Zip and Fix the output files
	sed 's/\\x3b/-/g'  ${VCF_somatic_annotation_path}/${TCGA_barcode}_indels_annotated.HG_multianno.vcf | sed 's/\\x3d/:/g' | gzip -c > ${VCF_somatic_annotation_path}/${TCGA_barcode}_indels_annotated.HG_multianno.vcf.gz
	sed 's/\\x3b/-/g'  ${VCF_somatic_annotation_path}/${TCGA_barcode}_indels_annotated.HG_multianno.txt | sed 's/\\x3d/:/g' | gzip -c > ${VCF_somatic_annotation_path}/${TCGA_barcode}_indels_annotated.HG_multianno.txt.gz
	# Remove temp files
	rm ${VCF_somatic_annotation_path}/${TCGA_barcode}_indels_annotated.HG_multianno.vcf ${VCF_somatic_annotation_path}/${TCGA_barcode}_indels_annotated.HG_multianno.txt ${tmp_dir}/${TCGA_barcode}_indels_somatic.vcf ${tmp_dir}/${TCGA_barcode}_indels_somatic_GT_fixed.vcf 
	gzip ${VCF_somatic_annotation_path}/${TCGA_barcode}_indels_annotated.avinput
fi
