#!/usr/bin/bash

#SBATCH --partition=normal_prio
#SBATCH --job-name=VCF_merge_GTEx
#SBATCH --output=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.out
#SBATCH --error=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
##SBATCH --nodes=X
##SBATCH --nodelist=fsupeksvr
##SBATCH --mem-per-cpu=X
#SBATCH --mem=50GB
#SBATCH --time=10:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=guillermo.palou@irbbarcelona.org

output_path="/g/strcombio/fsupek_cancer1/gpalou/GTEx/WES_VCF/rare_germline_variants/CADD_scores"
input_path="/g/strcombio/fsupek_cancer1/gpalou/GTEx/WES_VCF/rare_germline_variants"
/g/strcombio/fsupek_cancer1/gpalou/CADD_scores/CADD.sh -a -g GRCh38 -o ${output_path}/all_samples_liftOver_annotated_rare_variants_CADD.tsv.gz ${input_path}/all_samples_liftOver_annotated_rare_variants.vcf.gz
