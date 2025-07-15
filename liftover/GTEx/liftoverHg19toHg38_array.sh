#!/usr/bin/bash

#SBATCH --partition=normal_prio
#SBATCH --job-name=GTEx_liftover_array
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
#SBATCH --array=1-979%10

#GTEx_samples_list="/g/strcombio/fsupek_cancer1/gpalou/GTEx/WES_VCF/list_samples.txt"
GTEx_samples_list=$1

GTEx_sample=$(head -n${SLURM_ARRAY_TASK_ID} ${GTEx_samples_list} | tail -n1)

sbatch /home/gpalou/projects/NMD/scripts/liftover/GTEx/liftoverHg19toHg38_sample_vcf.sh ${GTEx_sample}