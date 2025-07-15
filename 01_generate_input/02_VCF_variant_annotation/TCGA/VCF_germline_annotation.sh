#!/usr/bin/bash

#SBATCH --partition=normal_prio
#SBATCH --job-name=annovar_tissue
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

TCGA_project=$1

while read VCF_germline_sample
do
        VCF_germline_annotation_path=$(echo $VCF_germline_sample | sed "s/.*\/${TCGA_project}\(.*\)\/variants.vcf.gz/\1/g")
        VCF_germline_annotation_path="/g/strcombio/fsupek_cancer1/gpalou/TCGA_germline_variants/"${TCGA_project}$VCF_germline_annotation_path
        mkdir -p $VCF_germline_annotation_path
	if [ -f "${VCF_germline_annotation_path}/${TCGA_project}_annotated.HG_multianno.vcf.gz" ]; then
                echo "Already done, skipping -->" ${VCF_germline_annotation_path}
        else
		sbatch VCF_germline_annotation_sample.sh $TCGA_project $VCF_germline_sample $VCF_germline_annotation_path
                sleep 2
        fi
done < ~/data/TCGA_germline_variants/${TCGA_project}_VCFs_samples_list.txt

