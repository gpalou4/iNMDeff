#!/usr/bin/bash

#SBATCH --partition=normal_prio
#SBATCH --job-name=fastq_decompress
##SBATCH --output=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.out
##SBATCH --error=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
##SBATCH --nodes=X
##SBATCH --nodelist=fsupeksvr
##SBATCH --mem-per-cpu=X
#SBATCH --mem=3GB
#SBATCH --time=1:00:00
#SBATCH --mail-type=ALL
##SBATCH --mail-user=guillermo.palou@irbbarcelona.org

## It says how much time it will take the job to be scheduled
##SBATCH --test-only

FASTQ_FILE=$1

if [[ $FASTQ_FILE == *.gz ]]; then
    tar -zxvf $FASTQ_FILE
else
    tar -xvf $FASTQ_FILE; gunzip -d *fastq.gz
fi

