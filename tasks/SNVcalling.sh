#!/bin/sh

#########################
### TASK: SNV calling ###
#########################

# Haplotipe Caller
# Genotype Caller


###############
## ARGUMENTS ##
###############

local=$1
run=$2
MDAP=$3
sample=$4
REF=$5
bamfile=$6
intervals=$7
panel=$8
cvcf=$9
removebam=$10





#####################
## MODULES AND DBs ##
#####################


if [ "$local" != "True" ]; then

	module load samtools/1.9
	module load picard/2.18.9
	module load gatk/4.1.2.0
	module load bedtools/2.27.0
	module load R
	alias picard='java -jar /usr/local/bioinfo/picard-tools/2.18.9/picard.jar'
	alias gatk='java -jar /usr/local/bioinfo/gatk/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar'	


	softwareFile="${MDAP}/software_${run}.txt"
	title="SNV CALLING"
	if [ ! -f $softwareFile ] || [ `grep -q $title $softwareFile` ] ; then 

		printf "SNV CALLING:\n" >> ${softwareFile}
		module list 2>> ${softwareFile}
	
	fi


	MAF="/home/proyectos/bioinfo/fjd/MAF_FJD"


else


	export SFT=/mnt/genetica3/marius/pipeline_practicas_marius/software
	alias samtools='$SFT/samtools/samtools'
	alias picard='java -Xmx10g -jar $SFT/picard/build/libs/picard.jar'
	alias gatk='java  -Xmx10g -jar $SFT/gatk/build/libs/gatk-package-4.0.6.0-22-g9d9484f-SNAPSHOT-local.jar'

	softwareFile="${MDAP}/software_${run}.txt"
	title="SNV CALLING"
	
	if [ ! -f $softwareFile ] || [ `grep -q $title $softwareFile` ] ; then 

		printf "SNV CALLING:\n" >> ${softwareFile}
		
		printf "\nSamtools VERSION\n" >> ${softwareFile}
		samtools --version 2>&1 | head -n2 >> ${softwareFile}
		
		printf "\nPicard VERSION\n" >> ${softwareFile}
		picard MarkDuplicates --version  2>&1 | head -n2 >> ${softwareFile}
		
		printf "\nGATK VERSION\n" >> ${softwareFile}
		gatk ApplyBQSR 2>&1 | head -n4 | tail -n1 >> ${softwareFile}


	fi




fi




###############
## VARIABLES ##
###############



# Ref Fasta
fasta=$REF


# genotyping folder
GEN=$MDAP/genotyping 
mkdir $GEN


# Temporal Folder
TMP=$MDAP/${sample}_tmp
mkdir $TMP





##############
## COMMANDS ##
##############




printf "\n\n\n- HAPLOTYPECALLER (GATK)"
printf "\n--------------------------\n"





printf '\nGATK HaplotypeCallerGVCF for '${sample}' STARTS'
start=`date +%s`

if [ "$cvcf" != "True" ]; then 

	if [ "$intervals" != "True" ]; then
		gatk HaplotypeCaller --tmp-dir=$TMP \
		-R $fasta \
		-I $bamfile \
		-bamout $GEN/${sample}_bamout.bam \
		-O $GEN/${sample}.vcf.gz \
		-G AS_StandardAnnotation \
		-A FisherStrand -A StrandOddsRatio -A RMSMappingQuality -A MappingQualityRankSumTest -A ReadPosRankSumTest -A DepthPerSampleHC -A BaseQualityRankSumTest -A ExcessHet  \
		--annotate-with-num-discovered-alleles=true 
	else
		gatk HaplotypeCaller --tmp-dir=$TMP \
		-R $fasta \
		-I $bamfile \
		-bamout $GEN/${sample}_bamout.bam \
		-O $GEN/${sample}.vcf.gz \
		-G AS_StandardAnnotation \
		-A FisherStrand -A StrandOddsRatio -A RMSMappingQuality -A MappingQualityRankSumTest -A ReadPosRankSumTest -A DepthPerSampleHC -A BaseQualityRankSumTest -A ExcessHet \
		--annotate-with-num-discovered-alleles=true \
		-L $panel -ip 1000
	fi


else

	if [ "$intervals" != "True" ]; then
		gatk HaplotypeCaller --tmp-dir=$TMP \
		-R $fasta \
		-I $bamfile \
		-ERC GVCF \
		-bamout $GEN/${sample}_bamout.bam \
		-O $GEN/${sample}.g.vcf \
		-G AS_StandardAnnotation \
		-A FisherStrand -A StrandOddsRatio -A RMSMappingQuality -A MappingQualityRankSumTest -A ReadPosRankSumTest -A DepthPerSampleHC -A BaseQualityRankSumTest -A ExcessHet  \
		--annotate-with-num-discovered-alleles=true 
	else
		gatk HaplotypeCaller --tmp-dir=$TMP \
		-R $fasta \
		-I $bamfile \
		-ERC GVCF \
		-bamout $GEN/${sample}_bamout.bam \
		-O $GEN/${sample}.g.vcf \
		-G AS_StandardAnnotation \
		-A FisherStrand -A StrandOddsRatio -A RMSMappingQuality -A MappingQualityRankSumTest -A ReadPosRankSumTest -A DepthPerSampleHC -A BaseQualityRankSumTest -A ExcessHet \
		--annotate-with-num-discovered-alleles=true \
		-L $panel -ip 1000
	fi

	printf '\n'$GEN/${sample}.g.vcf >> $GEN/my_list_of_gvcfs_files_to_combine_$run.list

fi



if [ "$?" = "0" ]; then
	printf '\nEXIT STATUS: 0'
	printf '\nGATK HaplotypeCaller for '${sample}' DONE\n' 

	if [ "$removebam" = "single" ]; then
	 	basename=${bamfile%.*}
	 	rm ${basename}.bam
	 	rm ${basename}.bai
	fi

else
	printf "\nERROR: PROBLEMS WITH HAPLOTYPECALLER"
	exit 1
fi


end=`date +%s`
runtime=$((end-start))
printf '\nExecuting time: '$runtime 












printf "\n\n\n- MOSDEPTH FOR BAMOUT"
printf "\n-----------------------\n"


module load miniconda/3.6
module load bedtools


mosdepth --quantize 10: -n -x  $MAF/coverage/${sample}_cov ${GEN}/${sample}_bamout.bam

if [ "$?" = "0" ]; then
	printf '\nEXIT STATUS: 0'
	printf '\nMOSDEPTH for '${sample}' DONE\n' 
	rm ${MAF}/coverage/${sample}_cov.mosdepth.global.dist.txt
else
	printf "\nERROR: PROBLEMS WITH MOSDEPTH"
	exit 1
fi


# make padding file and intersect depth if intervals are given


if [ "$intervals" != "False" ]; then

	basenamePanel="$(basename -- $panel)"
	awk '{print $1"\t"$2-1000"\t"$3+1000}' $panel > ${MDAP}/${basenamePanel}_padding1000_tmp.bed
	bedtools merge -i ${MDAP}/${basenamePanel}_padding1000_tmp.bed > ${MDAP}/${basenamePanel}_padding1000_tmp_merged.bed
	bedtools intersect -a ${MDAP}/${basenamePanel}_padding1000_tmp_merged.bed -b $MAF/coverage/${sample}_cov.quantized.bed.gz -wb | awk '{print $1"\t"$2"\t"$3"\t"$7}' >  $MAF/coverage/${sample}_cov.quantized.bed_Genotyped.bed


	if [ "$?" = "0" ]; then
		printf '\nEXIT STATUS: 0'
		printf '\nBEDTOOLS INTERSECT for '${sample}' DONE\n' 
		rm ${panel}_padding1000_tmp.bed
		rm ${panel}_padding1000_tmp_merged.bed
	else
		printf "\nERROR: PROBLEMS WITH BEDTOOLS INTERSECT"
		exit 1
	fi


module unload miniconda/3.6
module unload bedtools







# printf "\n\n\n- GENOTYPECALLER (GATK) "
# printf "\n--------------------------\n"


# #GenotypeGVCFs into final VCF

# printf  '\nUsing GATK GenotypeGVCFs for final VCF'
# start=`date +%s`

# gatk GenotypeGVCFs --tmp-dir=$TMP \
# 	-R $fasta \
# 	-V $GEN/${sample}.g.vcf \
# 	-G StandardAnnotation \
# 	-O $GEN/${sample}.vcf 

# if [ "$?" = "0" ]; then
# 	printf '\nEXIT STATUS: 0'
# 	printf  '\nGATK GenotypeGVCFs for '${sample}' DONE'

# else
# 	printf "\nERROR: PROBLEMS WITH GENOTYPECALLER"
# 	exit 1
# fi


# end=`date +%s`
# runtime=$((end-start))
# printf '\nExecuting time: '$runtime 




