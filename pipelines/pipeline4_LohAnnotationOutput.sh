#!/bin/sh

###########################################################################
#### Pipeline  4: LOH - Annotation - MAF Filtering - Output processing ####
###########################################################################

### PIPELINE ARGUMENTS 

MDAP=$1
name=$2 
threads=$3
run=$4
panel=$5
cvcf=$6 
REF=$7
local=$8 
pathology=$9
genefilter=${10}
single=${11}
softwarePath=${12}
tasksPath=${softwarePath}/tasks
mafincorporation=${13}






printf "\n.............................\n"
printf "  VARIANT ANNOTATION \n"
printf ".............................\n"


vcf="${MDAP}/snvs/${name}.gatkLabeled.vcf"
softwareFile="${MDAP}/software_${run}.txt"


$tasksPath/VEPannotation.sh $threads $vcf $softwareFile
if [ "$?" != "0" ]  ; then exit 1; fi

vepoutput=${MDAP}/snvs/${name}.annotated.MAFfiltered.vcf






printf "\n...............\n"
printf "  LOH \n"
printf ".................\n"

$tasksPath/LOH.sh $local $run $MDAP $name $cvcf $REF


if [ "$?" = "0"  ]; then
	echo -e  '\nEXIT STATUS: 0'
	echo -e   '\nPLINK HOMOZYGOSITY for '${run}' DONE'
else
	echo -e  "ERROR: PROBLEMS WITH PLINK"
	exit 1
fi 






printf "\n.............................\n"
printf " VEP VCF TO TXT FILE \n"
printf ".............................\n"


printf '\nFrom VEP annotated VCF to txt file...\n'

input=$vepoutput
output=${input%.*}.txt
python $tasksPath/vep2tsv_woFreq.py $input -o $output -t -g -f GT,AD,DP

if [ "$?" = "0" ]; then
	printf '\nEXIT STATUS: 0'
	printf '\nVCF PROCESSING for '${sample}' DONE'

else
	printf "\nERROR: PROBLEMS WITH VCF PREPROCESSING"
	exit 1
fi





printf "\n.............................\n"
printf " POST-VEP ANNOTATIONS  \n"
printf ".............................\n"


printf '\nAnnotating extra features as Spanish Frequency and LOH...\n'

PLINK=$MDAP/plink

input=$output
output=${output%.*}.pvm.txt
python $tasksPath/PVM_Cluster.py ${input} -l $local -o ${output} -P "healthy" -k $PLINK/${name}.hom

if [ "$?" = "0" ]; then
	printf '\nEXIT STATUS: 0'
	printf '\nPOST-VEP MODIFICATIONS for '${sample}' DONE'
	rm $VCF_FINAL

else
 	printf "\nERROR: PROBLEMS WITH POST-VEP MODIFICATIONS"
 	exit 1

 fi




printf "\n.............................\n"
printf " FILTER VCF BY GENE LIST \n"
printf ".............................\n"



# Filter PVM files based on gene list

if [ "$genefilter" != "False" ]; then

	printf '\nFiltering VPM by gene list...\n'
	
 	input=$output
 	output=${output%.*}.genelist.txt
 	python $tasksPath/filtering_geneList.py -i $input -f ${genefilter} -o $output
 	if [ "$?" = "0" ]; then
 		printf '\nEXIT STATUS: 0'
 		printf '\nVARIANT FILTERING BY GENE LIST for '${sample}' DONE'

 	else
 		printf "\nERROR: PROBLEMS WITH GENE FILTERING"
 		exit 1

 	fi

 fi




if [ "${mafincorporation}" != "False" ]; then


	printf "\n.............................\n"
	printf " MOVING FILES TO MAF FOLDER\n"
	printf ".............................\n"


	
	# Moving vcf and coverage files to MAF folder


	printf '\nMoving files to MAF folder...\n'

	source ../pipeline.config

	maf_coverage_path=${db_coverage_path}
	maf_vcf_path=${db_vcf_path}

	#maf_coverage_path="/home/proyectos/bioinfo/fjd/MAF_FJD_v3.0/coverage/new_bed/"
	#maf_vcf_path="/home/proyectos/bioinfo/fjd/MAF_FJD_v3.0/individual_vcf/new_vcf/"


	if [ "$cvcf" != "False" ]; then

		module load bcftools

		vcf_cvcf="${MDAP}/snvs/${name}.final.vcf"

		for sample in `bcftools query -l ${vcf_cvcf}`
		do

			vcf=${MDAP}/snvs/${sample}.final.vcf 
			coverage="${MDAP}/qc/${sample}_padding.quantized.bed"
			
			cp $vcf $maf_vcf_path
			if [ "$?" != "0" ]  ; then exit 1; printf "\nERROR: PROBLEMS COPYING VCF FILE TO MAF"; fi

			cp $coverage $maf_coverage_path
			if [ "$?" != "0" ]  ; then exit 1; printf "\nERROR: PROBLEMS COPYING COVERAGE FILE TO MAF"; fi

			if [ "$single" != "True" ]; then rm ${vcf}*; fi

		done

	else

		vcf="${MDAP}/snvs/${name}.final.vcf"
		coverage="${MDAP}/qc/${name}_padding.quantized.bed"

		cp $vcf $maf_vcf_path
		if [ "$?" != "0" ]  ; then exit 1; printf "\nERROR: PROBLEMS COPYING VCF FILE TO MAF"; fi

		cp $coverage $maf_coverage_path
		if [ "$?" != "0" ]  ; then exit 1; printf "\nERROR: PROBLEMS COPYING COVERAGE FILE TO MAF"; fi

	fi

fi











