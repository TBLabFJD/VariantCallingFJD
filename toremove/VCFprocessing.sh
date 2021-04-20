#!/bin/sh

############################
### TASK: VCF processing ###
############################

# Given an annotated VCF, transforms into txt file and annotates spanish frequency and LOH.


###############
## ARGUMENTS ##
###############


local=$1
run=$2
MDAP=$3
name=$4 # sample name or run name
pathology=$5
genefilter=$6
tasksPath=$7





#####################
## MODULES AND DBs ##
#####################


if [ "$local" != "True" ]; then


	softwareFile="${MDAP}/software_${run}.txt"
	title="MAPPING"
	if [ ! -f $softwareFile ] || [ `grep -q $title $softwareFile` ] ; then 

		printf "VEP ANNOTATION:\n" >> ${softwareFile}
		module list 2>> ${softwareFile}
	
	fi




else


	softwareFile="${MDAP}/software_${run}.txt"
	title="MAPPING"
	
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
SNV="${MDAP}/snvs"
PLINK="${MDAP}/plink"

# Temporal Folder
TMP=$MDAP/${name}_tmp
mkdir $TMP



# Input
VCF="${SNV}/${name}.annotated.MAFfiltered.vcf" 





##############
## COMMANDS ##
##############





start=`date +%s`

printf '\nFrom VEP annotated VCF to txt file...\n'
start=`date +%s`


python $taskPath/vep2tsv_woFreq.py $VCF \
-o ${VCF%.*}.txt -t -g -f GT,AD,DP

if [ "$?" = "0" ]; then
	printf '\nEXIT STATUS: 0'
	printf '\nVCF PROCESSING for '${name}' DONE'

else
	printf "\nERROR: PROBLEMS WITH VCF PREPROCESSING"
	exit 1
fi


end=`date +%s`
runtime=$((end-start))
printf '\nExecuting time: '$runtime 







printf '\nAnnotating extra features as LOH...\n'

start=`date +%s`

python $taskPath/PVM_Cluster.py ${VCF%.*}.txt -l $local -k $PLINK/${name}.hom -o ${VCF%.*}.pvm.txt -P $pathology


if [ "$?" = "0" ]; then
	printf '\nEXIT STATUS: 0'
	printf '\nPOST-VEP MODIFICATIONS for '${name}' DONE'
	rm $VCF_FINAL

else
	printf "\nERROR: PROBLEMS WITH POST-VEP MODIFICATIONS"
	exit 1

fi

end=`date +%s`
runtime=$((end-start))
printf '\nExecuting time: '$runtime 




# Filter PVM files based on gene list

if [ "$genefilter" != "False" ]; then

	printf '\nFiltering VCF by gene list...\n'

	start=`date +%s`
	
	python $taskPath/filtering_geneList.py -i ${VCF%.*}.pvm.txt -f ${genefilter} -o ${VCF%.*}.pvm.genelist.txt


	if [ "$?" = "0" ]; then
		printf '\nEXIT STATUS: 0'
		printf '\nVARIANT FILTERING BY GENE LIST for '${name}' DONE'

	else
		printf "\nERROR: PROBLEMS WITH GENE FILTERING"
		exit 1

	fi

	end=`date +%s`
	runtime=$((end-start))
	printf '\nExecuting time: '$runtime 

fi



# Remove temporary files.

printf '\n'
rm -r $TMP
