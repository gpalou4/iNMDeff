#!/usr/bin/bash

#SBATCH --partition=normal_prio
#SBATCH --job-name=GTEx_filter_sample
#SBATCH --output=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.out
#SBATCH --error=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.err
#SBATCH --ntasks=1
##SBATCH --cpus-per-task=X
##SBATCH --nodes=X
##SBATCH --mem-per-cpu=X
#SBATCH --mem=4GB
#SBATCH --time=1:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=guillermo.palou@irbbarcelona.org

conda activate /home/gpalou/.conda/envs/pandas_bcftools

GTEx_sample=$1
MAF=$2
GTEx_samples_path="/g/strcombio/fsupek_cancer1/gpalou/GTEx/WES_VCF/VCFs"

# Check if sample already exists
bcftools filter -e 'INFO/AF="."' ${GTEx_samples_path}/${GTEx_sample}_liftOver_annotated.HG_multianno.vcf.gz | \
bcftools filter -i 'FILTER="PASS"' | \
bcftools filter -i '
	INFO/ExonicFunc.refGene="stopgain" | INFO/ExonicFunc.refGene="nonsynonymous_SNV" 
	| INFO/ExonicFunc.refGene="startloss" | INFO/ExonicFunc.refGene="stoploss"
	| INFO/ExonicFunc.refGene="frameshift_deletion" | INFO/ExonicFunc.refGene="frameshift_insertion"
	| INFO/ExonicFunc.refGene="nonframeshift_deletion" | INFO/ExonicFunc.refGene="nonframeshift_insertion"
	| INFO/Func.refGene="splicing" | INFO/Func.refGene="exonic-splicing" ' |  \
	sed -r 's/##INFO=<ID=(AF.*),Number=.,Type=String,/##INFO=<ID=\1,Number=.,Type=Float,/' > \
	${GTEx_samples_path}/${GTEx_sample}_liftOver_annotated.HG_multianno_tmp.vcf
gzip ${GTEx_samples_path}/${GTEx_sample}_liftOver_annotated.HG_multianno_tmp.vcf
bcftools filter -i "INFO/AF<=${MAF}" ${GTEx_samples_path}/${GTEx_sample}_liftOver_annotated.HG_multianno_tmp.vcf.gz \
	-O z -o ${GTEx_samples_path}/${GTEx_sample}_liftOver_annotated.HG_multianno_rare_variants.vcf.gz
rm ${GTEx_samples_path}/${GTEx_sample}_liftOver_annotated.HG_multianno_tmp.vcf.gz



