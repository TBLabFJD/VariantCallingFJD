#!/bin/sh

####################################
### TASK: Lost of Heterocigosity ###
####################################


###############
## ARGUMENTS ##
###############

local=$1
run=$2
MDAP=$3
sample=$4
cvcf=$5



#####################
## MODULES AND DBs ##
#####################


if [ "$local" != "True" ]; then

	alias plink='/usr/local/bioinfo/plink/plink'

	softwareFile="${MDAP}/software_${run}.txt"
	title="LOH"

	if [ ! -f $softwareFile ] || [ `grep -q $title $softwareFile` ] ; then 

		printf "LOH:\n" >> ${softwareFile}
		/usr/local/bioinfo/plink/plink --version >> ${softwareFile}
	
	fi




else

	alias plink='/usr/local/bioinfo/plink/plink'

	softwareFile="${MDAP}/software_${run}.txt"
	title="LOH"
	
	if [ ! -f $softwareFile ] || [ `grep -q $title $softwareFile` ] ; then 

		printf "LOH:\n" >> ${softwareFile}
		
		printf "\nSamtools VERSION\n" >> ${softwareFile}
		samtools --version 2>&1 | head -n2 >> ${softwareFile}


	fi




fi





###############
## VARIABLES ##
###############


# Ref Fasta
fasta=$REF/ucsc.hg19.fasta


# snv results folder
SNV="${MDAP}/snv_results"



# LOH folder
PLINK=$MDAP/plink 
mkdir $PLINK


# Temporal Folder
TMP=$MDAP/${sample}_tmp
mkdir $TMP





##############
## COMMANDS ##
##############



start=`date +%s`

params="--allow-extra-chr --homozyg --homozyg-window-het 1 --vcf-filter"

plink $params --vcf $SNV/${sample}_raw.vcf --out $PLINK/${sample} 1>&2



if [ "$?" = "0"  ]; then
	echo -e  '\nEXIT STATUS: 0'
	echo -e   '\nPLINK HOMOZYGOSITY for '${run}' DONE'
else
	echo -e  "ERROR: PROBLEMS WITH PLINK"
	exit 1
fi 


end=`date +%s`
runtime=$((end-start))
echo -e  '\nExecuting time: '$runtime 





if [ "$cvcf" = "True" ]; then 
	printf '\n' $PLINK/${sample}.hom >> $PLINK/LOH_to_combine_$run.list

fi

