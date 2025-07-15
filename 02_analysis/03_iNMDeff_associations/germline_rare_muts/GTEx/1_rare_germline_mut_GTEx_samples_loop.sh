#!/usr/bin/bash

#SBATCH --partition=normal_prio
##SBATCH --account=high_prio
#SBATCH --job-name=rare_germline_GTEx
#SBATCH --output=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.out
#SBATCH --error=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.err
#SBATCH --ntasks=1
##SBATCH --cpus-per-task=X
##SBATCH --nodes=X
##SBATCH --mem-per-cpu=X
#SBATCH --mem=4GB
#SBATCH --time=2:00:00
#SBATCH --mail-type=ALL
##SBATCH --mail-user=guillermo.palou@irbbarcelona.org

# Activate base environment

while read GTEx_sample;
do
    echo $GTEx_sample;
    # Check rare germline variants
    sbatch 1_rare_germline_mut_GTEx_sample.sh ${GTEx_sample}
done < /g/strcombio/fsupek_cancer1/gpalou/GTEx/WES_VCF/list_samples.txt


