#!/bin/sh

####################################################################
#### Pipeline  2: Preprocessing - SNV calling - VEP annotation ####
####################################################################


### PIPELINE ARGUMENTS 

INPUT=$1
MDAP=$2
sample=$3
threads=$4
run=$5
panel=$6
basespace=$7
cat=$8
fastqFolder=$9
analysis=${10}
cvcf=${11}
skipmapping=${12}    
REF=${13}
local=${14}
pathology=${15}
intervals=${16}
duplicates=${17}
removebam=${18}
genefilter=${19}
user=${20}
softwarePath=${21}
interval_n=${22}
tasksPath=${softwarePath}/tasks


if [ "$skipmapping" != "True" ]; then


	printf "\n.......................\n"
	printf "  PRE-PROCESSING $sample \n"
	printf ".........................\n"

	$tasksPath/preprocessing_BAM.sh $local $run $MDAP $sample $duplicates $REF

	FB=$MDAP/bams
	bamfile=$FB/${sample}_alignment.bam

else
	
	bamfile=${INPUT}/${sample}*.bam

fi



# printf "\n..........................\n"
# printf "  CALCULATE COVERAGE $sample \n"
# printf "............................\n"


# mosdepth --by $panel ${QC}/${sample} ${bamfile}

# cov=$(awk '{if($1=="total" && $3=="0.51"){print $2}}' ${QC}/${sample}.mosdepth.region.dist.txt)
# echo $sample'\t0.51\t'$cov'\t'$run'\t'$bamfile'\n' >> ${QC}/minCovFilterResults.txt

# if [ "$?" = "0" ]; then
# 	printf '\nEXIT STATUS: 0'
# 	printf '\nMOSDEPTH for '${sample}' DONE\n' 
# else
# 	printf "\nERROR: PROBLEMS WITH MOSDEPTH"
# 	exit 1
# fi





if [ "$analysis" = "mapping" ]; then

	exit 0

fi




printf "\n\n\n\n.......................\n"
printf "  SNV CALLING $sample \n"
printf ".........................\n"


$tasksPath/SNVcalling.sh $local $run $MDAP $sample $REF $bamfile $intervals $panel $cvcf $removebam $interval_n




if [ "$cvcf" = "True" ]; then

	exit 0

fi






printf "\n\n..........................\n"
printf "  VARIANT FILTERING $sample \n"
printf ".............................\n"


# CNN / HF

ftype="HF"
#ftype="CNN"

$tasksPath/SNVfiltering.sh $local $run $MDAP $sample $REF $ftype $cvcf








printf "\n\n\n\n.......................\n"
printf "  HOMOCIGOSITY $sample \n"
printf ".........................\n"


$tasksPath/LOH.sh $local $run $MDAP $sample $cvcf









if [ "$cvcf" = "True" ]; then

	rm -r $TMP
	exit 0

fi









printf "\n\n.......................\n"
printf "\n\n\n- VARIANT ANNOTATION (VEP - ENSEMBL) "
printf "\n---------------------------------------\n"


$tasksPath/VEP_annotation.sh $local $run $MDAP $sample $threads












printf "\n\n.......................\n"
printf "\n\n\n- OUTPUT PROCESSING "
printf "\n-------------------------\n"




## VCF TO TXT


start=`date +%s`

printf '\nFrom VEP annotated VCF to txt file...\n'
start=`date +%s`


python $utilitiesPath/vep2tsv_woFreq.py $VCF_FILTERED_3 \
-o $VCF_FINAL -t -g -f GT,AD,DP

if [ "$?" = "0" ]; then
	printf '\nEXIT STATUS: 0'
	printf '\nVCF PROCESSING for '${sample}' DONE'

else
	printf "\nERROR: PROBLEMS WITH VCF PREPROCESSING"
	exit 1
fi


end=`date +%s`
runtime=$((end-start))
printf '\nExecuting time: '$runtime 









## FEATURE ANNOTATION as LOH and SPANISH FREQUENCY



printf '\nAnnotating extra features...\n'

start=`date +%s`

#python $utilitiesPath/LOHmerge.py $PLINK/${sample}.hom $VCF_FINAL ${VCF_FINAL_PVM}

echo python $utilitiesPath/PVM_Cluster.py $VCF_FINAL -l $local -k $PLINK/${sample}.hom -o ${VCF_FINAL_PVM} -P $pathology

python $utilitiesPath/PVM_Cluster.py $VCF_FINAL -l $local -k $PLINK/${sample}.hom -o ${VCF_FINAL_PVM} -P $pathology


if [ "$?" = "0" ]; then
	printf '\nEXIT STATUS: 0'
	printf '\nPOST-VEP MODIFICATIONS for '${sample}' DONE'
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

	
	printf '\nFiltering variants by gene list...\n'


	start=`date +%s`

	VCF_FINAL_PVM_FILTER="${VEPVCFAD}/${sample}_filteredAnnotated_pvm_GENELIST.txt"
	
	python $utilitiesPath/filtering_geneList.py -i ${VCF_FINAL_PVM} -f ${genefilter} -o ${VCF_FINAL_PVM_FILTER}


	if [ "$?" = "0" ]; then
		printf '\nEXIT STATUS: 0'
		printf '\nVARIANT FILTERING BY GENE LIST for '${sample}' DONE'

	else
		printf "\nERROR: PROBLEMS WITH GENE FILTERING"
		exit 1

	fi

	end=`date +%s`
	runtime=$((end-start))
	printf '\nExecuting time: '$runtime 

fi









printf "\n\n..........................\n"
printf "\n\n\n- MAF FILE PROCESSING "
printf "\n------------------------------\n"


$tasksPath/MAF_fjd.sh $local $run $MDAP $sample $REF $intervals $panel $cvcf $removebam








# REMOVE TEMPORARY DIRECTORY

printf '\n'
rm -r $TMP