#!/usr/bin/bash

#SBATCH --partition=normal_prio
#SBATCH --job-name=annovar_sample
#SBATCH --output=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.out
#SBATCH --error=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.err
#SBATCH --ntasks=1
##SBATCH --cpus-per-task=X
##SBATCH --nodes=X
##SBATCH --mem-per-cpu=X
#SBATCH --mem=30GB
#SBATCH --time=5:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=guillermo.palou@irbbarcelona.org

# Activate base environment
conda activate base

TCGA_project=$1
VCF_germline_sample=$2
VCF_germline_annotation_path=$3

# Check if sample already exists

if [ ! -f "${VCF_germline_annotation_path}/${TCGA_project}_annotated.HG_multianno.vcf.gz" ]; then
	# Annotate sample VCF with annovar (hg38 + gnomAD and exac Allele Frequencies)
	# Beware the refGene is actually GENCODE/ENSEMBL, I built the custom gene annotation as stated here --> /g/strcombio/fsupek_cancer1/gpalou/human_genome/humandb_annovar_hg38/commands.txt
	~/software/annovar/table_annovar.pl -buildver HG -out ${VCF_germline_annotation_path}/${TCGA_project}_annotated -remove -protocol refGene,gnomad211_exome,exac03 -operation g,f,f -nastring . -polish -vcfinput ${VCF_germline_sample} /g/strcombio/fsupek_cancer1/gpalou/human_genome/humandb_annovar_hg38
	# Zip and Fix the output files
	sed 's/\\x3b/-/g'  ${VCF_germline_annotation_path}/${TCGA_project}_annotated.HG_multianno.vcf | sed 's/\\x3d/:/g' | gzip -c > ${VCF_germline_annotation_path}/${TCGA_project}_annotated.HG_multianno.vcf.gz
	sed 's/\\x3b/-/g'  ${VCF_germline_annotation_path}/${TCGA_project}_annotated.HG_multianno.txt | sed 's/\\x3d/:/g' | gzip -c > ${VCF_germline_annotation_path}/${TCGA_project}_annotated.HG_multianno.txt.gz
	rm ${VCF_germline_annotation_path}/${TCGA_project}_annotated.HG_multianno.vcf ${VCF_germline_annotation_path}/${TCGA_project}_annotated.HG_multianno.txt
	gzip ${VCF_germline_annotation_path}/${TCGA_project}_annotated.avinput
	# Obtain information from INFO field and create a new tab file (not needed because it's already done in the .txt file)
	#vcftools --vcf ${VCF_germline_annotation_path}/${TCGA_project}_annotated.HG_multianno.vcf.gz --get-INFO Func.refGene --get-INFO Gene.refGene --get-INFO GeneDetail.refGene --get-INFO ExonicFunc.refGene --get-INFO AAChange.refGene --get-INFO AF --get-INFO ExAC_ALL --out ${VCF_germline_annotation_path}/${TCGA_project}_annotated.HG_multianno_filt.txt

fi
