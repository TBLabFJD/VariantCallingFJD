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
name=$4
cvcf=$5
fasta=$6



#####################
## MODULES AND DBs ##
#####################


		
module load gatk/4.2.0
module load bcftools
module unload python/2.7.15

alias gatk='java -jar /usr/local/bioinfo/gatk/4.2.0/gatk-package-4.2.0.0-local.jar'	
alias plink='/usr/local/bioinfo/plink/plink'

softwareFile="${MDAP}/software_${run}.txt"
title="LOH"

if [ ! -f $softwareFile ] || ! grep -q $title $softwareFile  ; then
	printf "LOH:\n" >> ${softwareFile}
	/usr/local/bioinfo/plink/plink --version >> ${softwareFile}

fi










###############
## VARIABLES ##
###############




# snv results folder
SNV="${MDAP}/snvs"



# LOH folder
PLINK=$MDAP/plink 
mkdir $PLINK




##############
## COMMANDS ##
##############



start=`date +%s`
params="--allow-extra-chr --homozyg --homozyg-window-het 1 --vcf-filter"



if [ "$cvcf" = "True" ]; then 
	

	# if cvcf we need to split vcf into individual samples
	printf '\n\nSplitting VCF into individual samples to compute LOH...'


	for sample in `bcftools query -l $SNV/${name}.gatkLabeled.vcf`
	do
		gatk SelectVariants --exclude-non-variants -R $fasta \
		-V $SNV/${name}.gatkLabeled.vcf \
		-O $SNV/${sample}.gatkLabeled.vcf -sn $sample --exclude-non-variants

		if [ "$?" = "0"  ]; then
			echo -e  '\nEXIT STATUS: 0'
			echo -e  '\nVCF splitting for '${sample}' DONE'
		else
			echo -e  "ERROR: PROBLEMS WITH VCF SPLITTING FOR SAMPLE '${sample}'"
			exit 1
		fi 


		plink $params --vcf $SNV/${sample}.gatkLabeled.vcf  --out $PLINK/${sample} 1>&2

		if [ "$?" = "0"  ]; then
			echo -e  '\nEXIT STATUS: 0'
			echo -e   '\nPLINK HOMOZYGOSITY for '${sample}' DONE'
			sed 1d $PLINK/${sample}.hom >> $PLINK/${run}.hom
			#rm $PLINK/${sample}.*

		else
			echo -e  "ERROR: PROBLEMS WITH PLINK FOR SAMPLE '${sample}'"
			exit 1
		fi 


		if [ "$single" != "True" ]; then
			rm $SNV/${sample}.gatkLabeled.vcf*
		fi

	done 



	sed -i '1 i\header' $PLINK/${run}.hom



else 


	plink $params --vcf $SNV/${name}.gatkLabeled.vcf --out $PLINK/${name} 1>&2



	if [ "$?" = "0"  ]; then
		echo -e  '\nEXIT STATUS: 0'
		echo -e   '\nPLINK HOMOZYGOSITY for '${name}' DONE'
	else
		echo -e  "ERROR: PROBLEMS WITH PLINK"
		exit 1
	fi 


fi



end=`date +%s`
runtime=$((end-start))
echo -e '\nExecuting time: '$runtime 

#printf '\n' $PLINK/${sample}.hom >> $PLINK/LOH_to_combine_$run.list















