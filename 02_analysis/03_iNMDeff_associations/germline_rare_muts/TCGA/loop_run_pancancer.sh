#!/usr/bin/bash

# With gene list
# for randomization in {no,yes}
# do
#     for randomGenes in {no,yes}
#     do
#         for NMD_method in {ASE,endogenous}
#         do
#             for dataset in {PTV_Missense_CADD15_0.1perc,PTV_Missense_CADD25_0.1perc,PTV_0.1perc,CCR_90orhigher_0.1perc,MRT_25thCentile_0.1perc}
#             do
#                 sbatch TCGA_rare_germline_variants_dataset.sh pancancer $dataset $NMD_method 0.2 yes 2 $randomGenes $randomization
#             done
#         done
#     done
# done

# Without gene list
for randomization in {no,yes}
do
    for NMD_method in {ASE,endogenous}
    do
        for dataset in {PTV_Missense_CADD15_0.1perc,PTV_Missense_CADD25_0.1perc,PTV_0.1perc,CCR_90orhigher_0.1perc,MRT_25thCentile_0.1perc}
        do
            sbatch TCGA_rare_germline_variants_dataset.sh pancancer $dataset $NMD_method 0.2 no 2 no $randomization
        done
    done
done