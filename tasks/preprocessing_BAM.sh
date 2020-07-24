#!/bin/sh

####################################
### TASK: Preprocessing BAM file ###
####################################

# Mark Duplicates
# Bam sorting and indexing
# BQSR calibration



###############
## ARGUMENTS ##
###############

local=$1
run=$2
MDAP=$3
sample=$4
duplicates=$5
REF=$6



#####################
## MODULES AND DBs ##
#####################


if [ "$local" != "True" ]; then

	module load samtools/1.9
	module load picard/2.18.9
	module load gatk/4.1.5.0
	module load bedtools/2.27.0
	module load R
	alias picard='java -jar /usr/local/bioinfo/picard-tools/2.18.9/picard.jar'
	alias gatk='java -jar /usr/local/bioinfo/gatk/4.1.5.0/gatk-package-4.1.5.0-local.jar'	


	softwareFile="${MDAP}/software_${run}.txt"
	title="PREPROCESSING"
	if [ ! -f $softwareFile ] || [ `grep -q $title $softwareFile` ] ; then 

		printf "PREPROCESSING ALIGNMENT FILE:\n" >> ${softwareFile}
		module list 2>> ${softwareFile}
	
	fi




else

	export SFT=/mnt/genetica3/marius/pipeline_practicas_marius/software
	alias samtools='$SFT/samtools/samtools'
	alias picard='java -Xmx10g -jar $SFT/picard/build/libs/picard.jar'
	alias gatk='java  -Xmx10g -jar $SFT/gatk/build/libs/gatk-package-4.0.6.0-22-g9d9484f-SNAPSHOT-local.jar'

	softwareFile="${MDAP}/software_${run}.txt"
	title="PREPROCESSING"
	
	if [ ! -f $softwareFile ] || [ `grep -q $title $softwareFile` ] ; then 

		printf "PREPROCESSING ALIGNMENT FILE:\n" >> ${softwareFile}
		
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

# Alignment fastq folder
AF=$MDAP/alignmentfastq 
mkdir $AF

# Pre-processing folder
PPB=$MDAP/preprocessingbam 
mkdir $PPB


# Metrics folder
METRICS=$MDAP/metrics
mkdir $METRICS


# final BAM folder
FB=$MDAP/bams
mkdir $FB



# Temporal Folder
TMP=$MDAP/${sample}_tmp
mkdir $TMP


# Files
alignment=$AF/${sample}_mapped_merged.bam 
dedupped=$PPB/${sample}_dedupped.bam
bqsr=$FB/${sample}_alignment.bam







##############
## COMMANDS ##
##############


if [ "$duplicates" != "True" ]; then


	printf "\n\n\n- MARKING DUPLICATES (PICARD)"
	printf "\n--------------------------------\n"

	# Mark duplicate reads to avoid counting non-independent observations

	# Task is assuming query-sorted input so that the Secondary and Supplementary reads get marked correctly.
	# This works because the output of BWA is query-grouped and therefore, so is the output of MergeBamAlignment.
	# While query-grouped isn't actually query-sorted, it's good enough for MarkDuplicates with ASSUME_SORT_ORDER="queryname"
	
	start=`date +%s`

	bammerged_unsorted=$DD/dedupped_${sample}.bam

	printf '\nStart picard MarkDuplicates '${sample}''
	
	# lacking TMP folder 
	gatk MarkDuplicates --TMP_DIR=$TMP\
	-I $alignment \
	-O $dedupped \
	-M $METRICS/marked_dup_metrics_${sample}.txt \
	--REMOVE_DUPLICATES false --ASSUME_SORT_ORDER queryname \
	--OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 --VALIDATION_STRINGENCY SILENT



	if [ "$?" = "0" ]; then
		printf '\nEXIT STATUS: 0'
		printf  '\nPICARD MarkDuplicates '${sample}' DONE' 
		#rm $alignment

	else
		printf "\nERROR: PROBLEMS WITH MARKING DUPLICATES"
		exit 1
	fi


	end=`date +%s`
	runtime=$((end-start))
	printf '\nExecuting time: '$runtime 


else

	dedupped=$alignment

fi





printf "\n\n\n- SORTING BAM (PICARD) "
printf "\n------------------------\n"


# Sort BAM file by coordinate order and fix tag values for NM and UQ

printf '\nRun picard SortSam and SetNmMdAndUqTags for '${sample}'\n'
start=`date +%s`


deduppedsorted=${dedupped%.bam}_sorted.bam



gatk SortSam --TMP_DIR=$TMP\
  --INPUT ${dedupped} \
  --OUTPUT /dev/stdout \
  --SORT_ORDER "coordinate" \
  --CREATE_INDEX false \
  --CREATE_MD5_FILE false \
| \
gatk  SetNmMdAndUqTags --TMP_DIR=$TMP \
  --INPUT /dev/stdin \
  --OUTPUT ${deduppedsorted} \
  --CREATE_INDEX true \
  --CREATE_MD5_FILE true \
  --REFERENCE_SEQUENCE $fasta



if [ "$?" = "0" ]; then
	printf '\nEXIT STATUS: 0'
	printf '\nPicard SortSam for '${sample}'  DONE' 
	#rm ${dedupped} 


else
	printf "\nERROR: PROBLEMS WITH SORTING"
	exit 1
fi



end=`date +%s`
runtime=$((end-start))
printf '\nExecuting time: '$runtime 











printf "\n\n\n- BQSR RECALIBRATION TABLE"
printf "\n-------------------------------\n"


#Recalibrating  the reads using base quality score reads.
#GATK BaseRecalibration first table
#BaseRecalibration + table
printf '\nStarts GATK '${sample}' Recalibrator'
start=`date +%s`


#--known-sites $REF/dbsnp_138.hg19.vcf


# Known sites
myknwonsitesVCF=""
REFfolder=$(dirname "${REF}")
var="--known-sites"
for i in $REFfolder/*.vcf; do myknwonsitesVCF=`echo $var $i $myknwonsitesVCF`; done




# Can be done individually for each group read.  Generate the recalibration model by interval

gatk BaseRecalibrator --tmp-dir=$TMP \
-R $fasta \
-I $deduppedsorted \
--use-original-qualities \
$myknwonsitesVCF \
-O $METRICS/before_recalibrated_bqsr_data_${sample}.recal.table


if [ "$?" = "0" ]; then
	printf '\nEXIT STATUS: 0'
	printf  '\nGATK BaseRecalibrator '${sample}' DONE' 

else
	printf "\nERROR: PROBLEMS WITH BQSR RECALIBRATION TABLE"
	exit 1
fi


end=`date +%s`
runtime=$((end-start))
printf '\nExecuting time: '$runtime 










printf "\n\n\n- APPLYING BQSR  (GATK) "
printf "\n-------------------------\n"


#Applying the recalibration table to the bam file to continue the analysis.

printf '\nStarts '${sample}' BQSR'
start=`date +%s`


gatk ApplyBQSR --tmp-dir=$TMP \
-R $fasta \
-I $deduppedsorted  \
--bqsr $METRICS/before_recalibrated_bqsr_data_${sample}.recal.table \
-O $bqsr \
--static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
--add-output-sam-program-record \
--create-output-bam-md5 \
--use-original-qualities \
--create-output-bam-index


if [ "$?" = "0" ]; then
	printf '\nEXIT STATUS: 0'
	printf  '\nGATK ApplyBQSR '${sample}' DONE'
	#rm ${dedupped%.bam}_sorted.bam
	#rm ${dedupped%.bam}_sorted.bai
	#printf 'rm '$METRICS'/before_recalibrated_bqsr_data_'${sample}'.recal.table'

else
	printf "\nERROR: PROBLEMS WITH BQSR"
	exit 1
fi


end=`date +%s`
runtime=$((end-start))
printf '\nExecuting time: '$runtime 





# Removing temporal forder


rm -r $TMP