#!/usr/bin/bash

#SBATCH --partition=normal_prio
#SBATCH --job-name=send_jobs
#SBATCH --output=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.out
#SBATCH --error=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
##SBATCH --nodes=X
##SBATCH --mem-per-cpu=X
#SBATCH --nodelist=fsupeksvr2
##SBATCH --excludenode=fsupeksvr5
#SBATCH --mem=2GB
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=guillermo.palou@irbbarcelona.org

path="/g/strcombio/fsupek_cancer1/gpalou/NMD_project/NMD_targets/ENSEMBL_transcripts/NMD_global_2_shared_ensembl_final_old.txt"
awk 'NR > 1 {print $2}' $path | while read NMD_gene
do
    echo $NMD_gene
    sbatch /home/gpalou/projects/NMD/scripts/NMD_efficiency/analysis_results/Estimate_cNMDeff/Estimate_cNMDeff.sh "${NMD_gene}"
done