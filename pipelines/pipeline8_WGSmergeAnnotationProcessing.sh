#!/bin/sh

################################################################
#### Pipeline  7: MERGE ANNOTATED VCFS (WGS) and PROCESSING ####
################################################################

### PIPELINE ARGUMENTS 

MDAP=$1
name=$2 
threads=$3
run=$4
local=$5
softwarePath=${6}
genefilter=${7}
tasksPath=${softwarePath}/tasks




printf "\n.................................\n"
printf "MERGE CHR-BASED ANNOTATED VCFS  \n"
printf ".................................\n"



filestocombine=${MDAP}/snvs/${name}.list
printf '%s\n' "${MDAP}/snvs/${name}"*".final.vcf" > $filestocombine
output=${MDAP}/snvs/${name}.final.vcf
$tasksPath/mergeGVCF.sh $filestocombine $output $run $MDAP
if [ "$?" != "0" ]  ; then exit 1; fi


filestocombine=${MDAP}/snvs/${name}.list
printf '%s\n' "${MDAP}/snvs/${name}"*"annotated.vcf" > $filestocombine
output=${MDAP}/snvs/${name}.annotated.vcf
$tasksPath/mergeGVCF.sh $filestocombine $output $run $MDAP
if [ "$?" != "0" ]  ; then exit 1; fi



filestocombine=${MDAP}/snvs/${name}.list
printf '%s\n' "${MDAP}/snvs/${name}"*"annotated.MAFfiltered.vcf" > $filestocombine
output=${MDAP}/snvs/${name}.annotated.MAFfiltered.vcf
$tasksPath/mergeGVCF.sh $filestocombine $output $run $MDAP
if [ "$?" != "0" ]  ; then exit 1; fi






printf "\n.............................\n"
printf " VEP VCF TO TXT FILE \n"
printf ".............................\n"


start=`date +%s`

printf '\nFrom VEP annotated VCF to txt file...\n'
start=`date +%s`

outputtxt=${MDAP}/snvs/${name}.annotated.MAFfiltered.txt

python $tasksPath/vep2tsv_woFreq.py $output \
-o ${outputtxt} -t -g -f GT,AD,DP
if [ "$?" != "0" ]  ; then exit 1; fi





# printf "\n.............................\n"
# printf " POST-VEP ANNOTATIONS  \n"
# printf ".............................\n"


# start=`date +%s`

# printf '\nAnnotating extra features as Spanish Frequency...\n'
# start=`date +%s`

# outputpvm=${MDAP}/snvs/${name}.annotated.MAFfiltered.pvm.txt
# python $tasksPath/PVM_Cluster.py ${outputtxt} -l $local -o ${outputpvm} -P "healthy"
# if [ "$?" != "0" ]  ; then exit 1; fi







# printf "\n.............................\n"
# printf " FILTER VCF BY GENE LIST \n"
# printf ".............................\n"



# # Filter PVM files based on gene list

# if [ "$genefilter" != "False" ]; then

# 	printf '\nFiltering VCF by gene list...\n'

# 	start=`date +%s`
	
# 	python $tasksPath/filtering_geneList.py -i ${output%.*}.txt -f ${genefilter} -o ${output%.*}.genelist.txt
# 	if [ "$?" != "0" ]  ; then exit 1; fi


# fi

