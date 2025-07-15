#!/usr/bin/bash

WD="/home/gpalou/projects/NMD/scripts/NMD_efficiency/analysis_results/rare_germline_variants_associations/TCGA"

> ${WD}/sbatch_commands.txt

# With gene list
# while read cancer
# do
#     echo $cancer
#     for randomization in {no,yes}
#     do
#         for randomGenes in {no,yes}
#         do
#             for NMD_method in {ASE,endogenous}
#             do
#                 for dataset in {PTV_Missense_CADD15_0.1perc,PTV_Missense_CADD25_0.1perc,PTV_0.1perc,MRT_25thCentile_0.1perc,CCR_90orhigher_0.1perc}
#                 do
#                     echo "${WD}/TCGA_rare_germline_variants_dataset.sh ${cancer:5} $dataset $NMD_method 0.2 yes 2 $randomGenes $randomization" >> ${WD}/sbatch_commands.txt
#                 done
#             done
#         done
#     done
# done < /g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/TCGA_projects_names.txt

# Without gene list
while read cancer
do
    echo $cancer
    for randomization in {no,yes}
    do
        for NMD_method in {ASE,endogenous}
        do
            for dataset in {PTV_Missense_CADD15_0.1perc,PTV_Missense_CADD25_0.1perc,PTV_0.1perc,MRT_25thCentile_0.1perc,CCR_90orhigher_0.1perc}
            do
                echo "${WD}/TCGA_rare_germline_variants_dataset.sh ${cancer:5} $dataset $NMD_method 0.2 no 2 no $randomization" >> ${WD}/sbatch_commands.txt
            done
        done
    done
done < /g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/TCGA_projects_names.txt