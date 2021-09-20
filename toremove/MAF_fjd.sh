#!/bin/sh

###############################
### TASK: MOSTDEPTH calling ###
###############################


###############
## ARGUMENTS ##
###############

local=$1
run=$2
MDAP=$3
sample=$4
REF=$5
intervals=$6
panel=$7
cvcf=$8
removebam=$9




#####################
## MODULES AND DBs ##
#####################


if [ "$local" != "True" ]; then


	module load miniconda/3.6
	module load bedtools
	MAF="/home/proyectos/bioinfo/fjd/MAF_FJD"

	softwareFile="${MDAP}/software_${run}.txt"
	title="MAF FILE GENERATION"
	if [ ! -f $softwareFile ] || [ `grep -q $title $softwareFile` ] ; then 

		printf "MAF FILE GENERATION:\n" >> ${softwareFile}
		module list 2>> ${softwareFile}
	
	fi




else

	#¿PATHS?
	softwareFile="${MDAP}/software_${run}.txt"
	title="MAF FILE GENERATION"
	
	if [ ! -f $softwareFile ] || [ `grep -q $title $softwareFile` ] ; then 

		printf "VEP ANNOTATION:\n" >> ${softwareFile}

		printf "\nVEP VERSION\n" >> ${softwareFile}
		$VEP --help | grep "Versions:" -A 5 | tail -n4 >> ${softwareFile}

	fi




fi








###############
## VARIABLES ##
###############


# snv results folder
SNV="${MDAP}/snv_results"
GEN=${MDAP}/genotyping 
QC=${MDAP}/QC 


# Temporal Folder
TMP=$MDAP/${sample}_tmp
mkdir $TMP


VCF_RAW="${SNV}/${name}_raw.vcf"
BAMOUT="${GEN}/${sample}_bamout.bam"






printf "\n\n\n- MOSDEPTH FOR BAMOUT"
printf "\n-----------------------\n"


############### ADD FILE WITH SAMPLE, DATE, RUN, ETC.


info="#RUN:"${run}"_PANEL:"${panel}


mosdepth --quantize 10: -n -x  ${QC}/${sample}_cov ${BAMOUT}

if [ "$?" = "0" ]; then
	printf '\nEXIT STATUS: 0'
	printf '\nMOSDEPTH for '${sample}' DONE\n' 
	rm ${QC}/${sample}_cov.mosdepth.global.dist.txt
else
	printf "\nERROR: PROBLEMS WITH MOSDEPTH"
	exit 1
fi


# make padding file and intersect depth if intervals are given


if [ "$intervals" != "False" ]; then

	basenamePanel="$(basename -- $panel)"
	awk '{print $1"\t"$2-1000"\t"$3+1000}' $panel > ${QC}/${basenamePanel}_padding1000_tmp.bed
	bedtools merge -i ${QC}/${basenamePanel}_padding1000_tmp.bed > ${QC}/${basenamePanel}_padding1000_tmp_merged.bed
	bedtools intersect -a ${QC}/${basenamePanel}_padding1000_tmp_merged.bed -b ${QC}/${sample}_cov.quantized.bed.gz -wb | awk '{print $1"\t"$2"\t"$3"\t"$7}' >  ${QC}/${sample}_cov.quantized_Genotyped.bed


	if [ "$?" = "0" ]; then
		printf '\nEXIT STATUS: 0'
		printf '\nBEDTOOLS INTERSECT for '${sample}' DONE\n' 
		rm ${QC}/${basenamePanel}_padding1000_tmp.bed
		rm ${QC}/${basenamePanel}_padding1000_tmp_merged.bed
		rm ${QC}/${sample}_cov.quantized.bed.gz*
		sed -i 1i${info} ${QC}/${sample}_cov.quantized_Genotyped.bed
		cp   ${QC}/${sample}_cov.quantized_Genotyped.bed ${MAF}/coverage/.
		cp   ${VCF_RAW} ${MAF}/variants/.


	else
		printf "\nERROR: PROBLEMS WITH BEDTOOLS INTERSECT"
		exit 1
	fi


else

	sed  -i 1i${info} filename ${QC}/${sample}_cov.quantized.bed.gz
	cp ${QC}/${sample}_cov.quantized.bed.gz ${MAF}/coverage/${sample}_cov.quantized_Genotyped.bed.gz
	cp ${VCF_RAW} ${MAF}/variants/.


