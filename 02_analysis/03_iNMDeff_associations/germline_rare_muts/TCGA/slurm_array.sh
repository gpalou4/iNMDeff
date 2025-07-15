#!/usr/bin/bash

#SBATCH --partition=normal_prio
#SBATCH --job-name=RGWAS_TCGA_array
#SBATCH --output=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.out
#SBATCH --error=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.err
#SBATCH --ntasks=1
##SBATCH --cpus-per-task=X
##SBATCH --nodes=X
##SBATCH --mem-per-cpu=X
##SBATCH --exclude=fsupeksvr5
#SBATCH --mem=20GB
#SBATCH --time=5:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=guillermo.palou@irbbarcelona.org

commands_file=$1

l=`head -n${SLURM_ARRAY_TASK_ID} ${commands_file} | tail -n1`

eval $l

