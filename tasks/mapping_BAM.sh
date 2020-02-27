#!/bin/sh

######################################
### Task: Mapping and unmapped BAM ###
######################################

###############
## ARGUMENTS ##
###############

local=$1
run=$2
MDAP=$3
sample=$4
forward=$5
reverse=$6
REF=$7
threads=$8


#####################
## MODULES AND DBs ##
#####################


if [ "$local" != "True" ]; then

	module load bwa/0.7.17
	module load gatk/4.1.2.0
	module load samtools/1.9
	alias gatk='java -jar /usr/local/bioinfo/gatk/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar'	
	alias picard='java -jar /usr/local/bioinfo/picard-tools/2.18.9/picard.jar'


	softwareFile="${MDAP}/software_${run}.txt"
	title="MAPPING"
	if [ ! -f $softwareFile ] || [ `grep -q $title $softwareFile` ] ; then 

		printf "MAPPING:\n" >> ${softwareFile}
		module list 2>> ${softwareFile}
	
	fi




else

	export SFT=/mnt/genetica3/marius/pipeline_practicas_marius/software
	alias bwa='$SFT/bwa/bwa'
	alias gatk='java  -Xmx10g -jar $SFT/gatk/build/libs/gatk-package-4.0.6.0-22-g9d9484f-SNAPSHOT-local.jar'

	softwareFile="${MDAP}/software_${run}.txt"
	title="MAPPING"
	
	if [ ! -f $softwareFile ] || [ `grep -q $title $softwareFile` ] ; then 

		printf "MAPPING:\n\n" >> ${softwareFile}
		
		printf "BWA VERSION\n" >> ${softwareFile}

		bwa 2>&1 | head -n4 | tail -n2 | head -n1 >> ${softwareFile}
		
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


# Temporal Folder
TMP=$MDAP/${sample}_tmp
mkdir $TMP








#############
## COMMAND ##
#############



printf "\n\n\n- INDEXING REFERENCE FILES (BWA) "
printf "\n----------------------------------\n"

# BWA index
if [ ! -f ${REF}.sa ]; then
	printf '\nRUN BWA INDEX!!!'
	exit 1
	#printf 'Starts BWA INDEX'
	#bwa index $HG19/ucsc.hg19.fasta
	#printf sampleFile '\nBWA INDEX DONE' 
else
	printf '\nBWA INDEX already existing. Reference indexing for BWA skipped'
fi




# Genome index
if [ ! -f ${REF}.fai ]; then
	printf '\nRUN REFERNCE INDEX!!!'
	exit 1		
	# printf '\nCreate .FAI file, using samtools faidx'
	# samtools faidx $HG19/ucsc.hg19.fasta 
	# printf sampleFile '\nucsc.hg19.fasta.FAI DONE'
else
	printf '\nREFERENCE INDEX already existing. Fasta indexing skipped'


fi



# Genome dict
REFext="${REF%.*}"
if [ ! -f ${REFext}.dict ]; then
	printf '\nCREATE REFERENCE DICT!!!'
	exit 1	
	# Creating .DICT in HG19.
	# printf 'Create .DICT file, using picardtools CreateSequnceDictionary'
	# picard CreateSequenceDictionary R=$HG19/ucsc.hg19.fasta O=$HG19/ucsc.hg19.2.dict
	# printf sampleFile 'ucsc.hg19.DICT DONE'
else
	printf '\nREFERENCE DICT already existing. Dictionary generation skipped'

fi









printf "\n\n\n- MAPPING READS (BWA) "
printf "\n----------------------------\n"

start=`date +%s`

#mapping fastq files to reference_genome after BWA INDEX

printf '\nStart '${sample}' mapping...\n'


bwa mem -v 3 -t $threads -Y \
$fasta \
$forward \
$reverse |  samtools view -1 > $AF/${sample}_mapped.bam


if [ "$?" = "0" ]; then
	printf '\nEXIT STATUS: 0'
	printf  '\nBWA MEM '${sample}'  DONE'
else
	printf "\nERROR: PROBLEMS WITH MAPPING"
	exit 1
fi

end=`date +%s`
runtime=$((end-start))
printf '\nExecuting time: '$runtime 











printf "\n\n\n- FASTQ TO uBWA "
printf "\n---------------------\n"


start=`date +%s`


gatk FastqToSam -TMP_DIR=$TMP \
--FASTQ ${forward} \
--FASTQ2 ${reverse} \
--OUTPUT $AF/${sample}_unmapped.bam \
--READ_GROUP_NAME ${sample} \
--SAMPLE_NAME ${sample} \
--PLATFORM illumina 


if [ "$?" = "0" ]; then
	printf '\nEXIT STATUS: 0'
	printf  '\nuBAM generation '${sample}'  DONE'
	if [ "$cat" = "True" ] || [ "$basespace" = "True" ]; then
		rm $forward
		rm $reverse
	fi
else
	printf "\nERROR: PROBLEMS WITH FASTQ TO SAM"
	exit 1
fi

end=`date +%s`
runtime=$((end-start))
printf '\nExecuting time: '$runtime 












printf "\n\n\n- Merge uBAM and BWA-aligned BAM "
printf "\n------------------------------------\n"


# Merge original input uBAM file with BWA-aligned BAM file

start=`date +%s`

echo $fasta

gatk MergeBamAlignment -TMP_DIR=$TMP \
  --VALIDATION_STRINGENCY SILENT \
  --EXPECTED_ORIENTATIONS FR \
  --ATTRIBUTES_TO_RETAIN X0 \
  --ALIGNED_BAM $AF/${sample}_mapped.bam \
  --UNMAPPED_BAM $AF/${sample}_unmapped.bam  \
  --OUTPUT $AF/${sample}_mapped_merged.bam \
  --REFERENCE_SEQUENCE $fasta \
  --PAIRED_RUN true \
  --SORT_ORDER "unsorted" \
  --IS_BISULFITE_SEQUENCE false \
  --ALIGNED_READS_ONLY false \
  --CLIP_ADAPTERS false \
  --ADD_MATE_CIGAR true \
  --MAX_INSERTIONS_OR_DELETIONS -1 \
  --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
  --UNMAPPED_READ_STRATEGY COPY_TO_TAG \
  --ALIGNER_PROPER_PAIR_FLAGS true \
  --UNMAP_CONTAMINANT_READS true

 # removed!
 # PROGRAM_GROUP_VERSION "${bwa_version}" \
 # PROGRAM_GROUP_COMMAND_LINE "${bwa_commandline}" \


if [ "$?" = "0" ]; then
	printf '\nEXIT STATUS: 0'
	printf  '\nuBAM and BWA-aligned BAM merging '${sample}'  DONE'
	rm $AF/${sample}_mapped.bam
	rm $AF/${sample}_unmapped.bam

else
	printf "\nERROR: PROBLEMS WITH MERGE BAM ALIGNMENT"
	exit 1
fi

end=`date +%s`
runtime=$((end-start))
printf '\nExecuting time: '$runtime 




