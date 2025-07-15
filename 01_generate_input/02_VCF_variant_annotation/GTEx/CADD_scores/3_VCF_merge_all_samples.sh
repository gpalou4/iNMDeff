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

conda activate /home/gpalou/.conda/envs/pandas_bcftools

# Generate list of samples
# > samples_list.txt
# while read GTEx_sample
# do    
#     echo "/g/strcombio/fsupek_cancer1/gpalou/GTEx/WES_VCF/VCFs/${GTEx_sample}_liftOver_annotated.HG_multianno_rare_variants.vcf.gz" >> samples_list.txt 
# done < /g/strcombio/fsupek_cancer1/gpalou/GTEx/WES_VCF/list_samples.txt

# Fixing the bgzip files -->

# while read GTEx_sample
# do    
#     echo $GTEx_sample
#     bgzip -d "/g/strcombio/fsupek_cancer1/gpalou/GTEx/WES_VCF/VCFs/${GTEx_sample}_liftOver_annotated.HG_multianno_rare_variants.vcf.gz"
#     bgzip /g/strcombio/fsupek_cancer1/gpalou/GTEx/WES_VCF/VCFs/${GTEx_sample}_liftOver_annotated.HG_multianno_rare_variants.vcf
#     tabix /g/strcombio/fsupek_cancer1/gpalou/GTEx/WES_VCF/VCFs/${GTEx_sample}_liftOver_annotated.HG_multianno_rare_variants.vcf.gz
# done < /g/strcombio/fsupek_cancer1/gpalou/GTEx/WES_VCF/list_samples.txt

output_path="/g/strcombio/fsupek_cancer1/gpalou/GTEx/WES_VCF/rare_germline_variants"

bcftools merge -l samples_list.txt -Oz -o ${output_path}/all_samples_liftOver_annotated_rare_variants.vcf.tmp.gz

# Change chromosome notation
bcftools annotate --rename-chrs chr_notation.txt ${output_path}/all_samples_liftOver_annotated_rare_variants.vcf.tmp.gz -O z -o ${output_path}/all_samples_liftOver_annotated_rare_variants.vcf.gz
rm ${output_path}/all_samples_liftOver_annotated_rare_variants.vcf.tmp.gz



