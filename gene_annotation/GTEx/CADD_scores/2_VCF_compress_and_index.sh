#!/usr/bin/bash

#SBATCH --partition=normal_prio
#SBATCH --job-name=GTEx_annovar
#SBATCH --output=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.out
#SBATCH --error=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
##SBATCH --nodes=X
##SBATCH --nodelist=fsupeksvr
##SBATCH --mem-per-cpu=X
#SBATCH --mem=2GB
#SBATCH --time=10:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=guillermo.palou@irbbarcelona.org

while read VCF_germline_sample
do
        sbatch VCF_compress_and_index_sample.sh $VCF_germline_sample
done < /g/strcombio/fsupek_cancer1/gpalou/GTEx/WES_VCF/list_samples.txt

