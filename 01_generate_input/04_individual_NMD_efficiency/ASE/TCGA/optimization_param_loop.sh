
# To exclude on variable just put xxx where the variable should be

cancer=$1
somatic_transcripts_filt=$2

methods=( strelka samtools )

for method in "${methods[@]}"
do
	while read VAF
	do
		sbatch ASE_NMD_efficiency.sh ${cancer} ${somatic_transcripts_filt} $VAF $method
	done < VAFs.txt
done
