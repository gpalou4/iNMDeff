#!/usr/bin/bash

#SBATCH --partition=downloads
#SBATCH --job-name=gdc_download_TCGA
#SBATCH --output=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.out
#SBATCH --error=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nice=100000
##SBATCH --nodes=X
##SBATCH --nodelist=fsupeksvr
##SBATCH --mem-per-cpu=X
#SBATCH --mem=12GB
#SBATCH --time=5:00:00
#SBATCH --mail-type=ALL
##SBATCH --mail-user=guillermo.palou@irbbarcelona.org

## It says how much time it will take the job to be scheduled
##SBATCH --test-only

# Activate base environment
conda activate /opt/anaconda3/envs/gdc-client

CPUS=$1
GDC_TOKEN=$2
GDC_DOWNLOAD_TRIES=$3
UUID=$4

gdc-client download -n $CPUS -t $GDC_TOKEN --retry-amount $GDC_DOWNLOAD_TRIES $UUID || exit 1
#gdc-client download -n !{task.cpus*3} -t !{gdctok} -d ./step_output_dir --retry-amount !{gdc_download_tries} !{gdcid_normalbam} || exit 1
