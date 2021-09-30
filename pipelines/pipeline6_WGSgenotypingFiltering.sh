#!/bin/sh

####################################################################
#### Pipeline  7: MERGE GVCFS (WGS), SNV CALLING AND fILTERING ####
####################################################################

### PIPELINE ARGUMENTS 

MDAP=$1
name=$2 
REF=$3
run=$4
local=$5
cvcf=$6
softwarePath=${7}
tasksPath=${softwarePath}/tasks




printf "\n.............................\n"
printf "MERGE INTERVAL GVCFS  \n"
printf "............................\n"


GEN=$MDAP/genotyping
mkdir $GEN 
filestocombine=$GEN"/"${name}.list
printf '%s\n' "$MDAP/interval_"*"/genotyping/${name}.g.vcf" > $filestocombine
output=$GEN/${name}.g.vcf


$tasksPath/mergeGVCF.sh $filestocombine $output $run $MDAP
if [ "$?" != "0" ]  ; then exit 1; else rm $filestocombine; rm $MDAP/interval_*/genotyping/${name}.g.vcf*; rm $MDAP/interval_*/genotyping/${name}_bamout*; rm $MDAP/interval_*/genotyping/my_list_of_gvcfs_files*.list; rm $MDAP/interval_*/software*txt; rm $MDAP/interval_*/qc/${name}*; rm $MDAP/interval_*/qc/minCovFilterResults_snv*; fi


 



printf "\n.............................\n"
printf " GENOTYPING \n"
printf ".............................\n"

input=$output
output=$GEN/${name}.vcf
 # pedrigree is null, no option yet for WGS
$tasksPath/genotyping.sh $local $run $MDAP $name $input $output $REF "null"
if [ "$?" != "0" ]  ; then exit 1; else rm ${input}*; fi





printf "\n\n......................\n"
printf "VARIANT FILTERING \n"
printf "........................\n"


# CNN / HF

ftype="HF"
#ftype="CNN"

$tasksPath/SNVfiltering.sh $local $run $MDAP $name $REF $ftype $cvcf $softwarePath
if [ "$?" != "0" ]  ; then exit 1; fi



