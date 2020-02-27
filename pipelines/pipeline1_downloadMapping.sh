#!/bin/sh

#####################################################
### Pipeline 1: Basemount, Mapping and unmapped BAM ###
#####################################################


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

tasksPath=${softwarePath}/tasks





printf "\n.........................\n"
printf "  DEFINING INPUT $sample \n"
printf ".........................\n"




if [ "$cat" = "True" ] || [ "$basespace" = "True" ]; then

	python $utilitiesPath/bscpCat_sample.py $basespace $fastqFolder $INPUT $sample $user

	if [ "$?" != "0" ]; then
		printf "\nERROR: PROBLEMS WITH BASESPACE DATA DOWNLOADING/CONCATENATION"
		exit 1
	else
		printf '\nEXIT STATUS: 0'
		printf '\nDOWNLOAD/CAT FOR '${sample}'  DONE'
	fi

	INPUT=$fastqFolder

	if [[ "$basespace" = "True" ]]; then
		sample=${sample:0:7}
	fi
fi

forward="${INPUT}/${sample}*_R1*.f*q.gz"
reverse="${INPUT}/${sample}*_R2*.f*q.gz"












printf "\n.....................\n"
printf "  MAPPING $sample \n"
printf ".....................\n"



$tasksPath/mapping_BAM.sh $local $run $MDAP $sample $forward $reverse $REF $threads


if [ "$cat" = "True" ] || [ "$basespace" = "True" ]; then
	rm $foward
	rm $reverse
fi


