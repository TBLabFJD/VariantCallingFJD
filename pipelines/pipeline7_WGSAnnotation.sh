#!/bin/sh

#####################################
#### Pipeline  4: Annotation WGS ####
#####################################

### PIPELINE ARGUMENTS 

MDAP=$1
name=$2 
threads=$3
run=$4
local=$5
softwarePath=${6}
interval=${7}
tasksPath=${softwarePath}/tasks




printf "\n.............................\n"
printf "  VARIANT ANNOTATION \n"
printf ".............................\n"


vcf="${MDAP}/snvs/${name}.gatkLabeled.vcf"
softwareFile="${MDAP}/software_${run}.txt"


$tasksPath/VEPannotation.sh $threads $vcf $softwareFile $interval
if [ "$?" != "0" ]  ; then exit 1; fi


