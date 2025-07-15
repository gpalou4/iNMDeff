#!/usr/bin/bash

#SBATCH --partition=normal_prio
#SBATCH --job-name=RSEM_reference_TCGA
#SBATCH --output=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.out
#SBATCH --error=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
##SBATCH --nodes=X
##SBATCH --nodelist=fsupeksvr
##SBATCH --mem-per-cpu=X
#SBATCH --mem=12GB
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=guillermo.palou@irbbarcelona.org

## It says how much time it will take the job to be scheduled
##SBATCH --test-only

# build the RSEM reference:
singularity run /g/strcombio/fsupek_home/gpalou/singularity/gtex_rnaseq_V8.sif \
    /bin/bash -c "rsem-prepare-reference \
        /home/gpalou/data/human_genome/TCGA/GRCh38.d1.vd1.fa \
        /home/gpalou/data/TCGA_RNAseq_quantification/rsem_reference/rsem_reference \
        --gtf /home/gpalou/data/NMD_project/conversor_tables/ENSEMBL/gencode.v26.annotation.gtf \
        --num-threads 4"
