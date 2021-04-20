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
softwarePath=${11}
tasksPath=${softwarePath}/tasks






printf "\n...............\n"
printf "  LOH \n"
printf ".................\n"

$tasksPath/LOH.sh $local $run $MDAP $name $cvcf $REF
if [ "$?" != "0" ]  ; then exit 1; fi






printf "\n.............................\n"
printf "  VARIANT ANNOTATION \n"
printf ".............................\n"


vcf="${MDAP}/snvs/${name}.gatkLabeled.vcf"
softwareFile="${MDAP}/software_${run}.txt"


$tasksPath/VEPannotation.sh $local $threads $vcf $softwareFile
if [ "$?" != "0" ]  ; then exit 1; fi







printf "\n.............................\n"
printf " VEP VCF TO TXT FILE \n"
printf ".............................\n"


start=`date +%s`

printf '\nFrom VEP annotated VCF to txt file...\n'
start=`date +%s`

input=${MDAP}/snvs/${name}.annotated.MAFfiltered.txt
outputtxt=${MDAP}/snvs/${name}.annotated.MAFfiltered.txt

python $tasksPath/vep2tsv_woFreq.py $input \
-o ${outputtxt} -t -g -f GT,AD,DP
if [ "$?" != "0" ]  ; then exit 1; fi





printf "\n.............................\n"
printf " POST-VEP ANNOTATIONS  \n"
printf ".............................\n"


start=`date +%s`

printf '\nAnnotating extra features as Spanish Frequency...\n'
start=`date +%s`

outputpvm=${MDAP}/snvs/${name}.annotated.MAFfiltered.pvm.txt
python $taskPath/PVM_Cluster.py ${outputtxt} -l $local -o ${outputpvm} -P "healthy"
if [ "$?" != "0" ]  ; then exit 1; fi







printf "\n.............................\n"
printf " FILTER VCF BY GENE LIST \n"
printf ".............................\n"



# Filter PVM files based on gene list

if [ "$genefilter" != "False" ]; then

	printf '\nFiltering VCF by gene list...\n'

	start=`date +%s`
	
	python $tasksPath/filtering_geneList.py -i ${output%.*}.txt -f ${genefilter} -o ${output%.*}.genelist.txt
	if [ "$?" != "0" ]  ; then exit 1; fi


fi









