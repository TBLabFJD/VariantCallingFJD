#!/bin/sh

############################################################
### FJD pipeline - Combined Genotyping and Trio analysis ###
############################################################






if [ "$local" != "True" ]; then

	module load samtools/1.9
	module load picard/2.18.9
	module load gatk/4.1.5.0
	module load bedtools/2.27.0
	module load R
	alias picard='java -jar /usr/local/bioinfo/picard-tools/2.18.9/picard.jar'
	alias gatk='java -jar /usr/local/bioinfo/gatk/4.1.5.0/gatk-package-4.1.5.0-local.jar'	


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





###########
## INPUT ##
###########

# options 


MDAP=$1
filetocombine=$2
sample=$3

# files

fasta=$4 # Ref Fasta



# genotyping folder

GEN=$MDAP/genotyping 
mkdir $GEN


# Temporal Folder

TMP=$MDAP/${sample}_tmp
mkdir $TMP




### RUN:

# Merge GVCFs generated per-interval for the same sample


echo -e "\n\n\n- merge GVCFS (GATK) "
echo -e "-----------------------------\n"

echo -e 'mkdir combined_gvcf_data'
mkdir $CGVCFD
echo -e '\nUsing GATK COMBINEGVCFs for merging GVCFs'
start=`date +%s`

gatk MergeVcfs  \
-I $filetocombine \
-O $GEN/$sample.g.vcf


if [ "$?" = "0" ]; then
	echo -e  '\nEXIT STATUS: 0'
	echo -e '\nGATK COMBINEGVCFs DONE'

else
	echo -e  "ERROR: PROBLEMS WITH COMBINEGVCFs"
	exit 1
fi

end=`date +%s`
runtime=$((end-start))
echo -e  '\nExecuting time: '$runtime 


