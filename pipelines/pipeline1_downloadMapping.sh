#!/bin/sh

########################################################
### Pipeline 1: Basemount, Mapping and unmapped BAM ###
########################################################


### PIPELINE ARGUMENTS 

INPUT=$1
MDAP=$2
sample=$3
threads=$4
run=$5
basespace=$6
cat=$7
fastqFolder=$8
REF=${9}
local=${10}
user=${11}
softwarePath=${12}
tasksPath=${softwarePath}/tasks





printf "\n.........................\n"
printf "DEFINING INPUT $sample \n"
printf ".........................\n"




if [ "$cat" = "True" ] || [ "$basespace" = "True" ]; then

	python $tasksPath/fastqDownloadCat.py $basespace $fastqFolder $INPUT $sample $user

	if [ "$?" != "0" ]; then
		printf "\nERROR: PROBLEMS WITH BASESPACE DATA DOWNLOADING/CONCATENATION"
		exit 1
	else
		printf '\nEXIT STATUS: 0'
		printf '\nDOWNLOAD/CAT FOR '${sample}'  DONE\n'
	fi

	INPUT=$fastqFolder

	# if [[ "$basespace" = "True" ]]; then
	# 	sample=${sample:0:7}
	# fi
fi

forward="${INPUT}/${sample}*_R1*.f*q.gz"
reverse="${INPUT}/${sample}*_R2*.f*q.gz"






printf "\n.....................\n"
printf "  MAPPING $sample \n"
printf ".....................\n"



$tasksPath/mapping.sh $run $MDAP $sample $forward $reverse $REF $threads
s1="$?"


if [ "$cat" = "True" ] || [ "$basespace" = "True" ]; then
	rm $forward
	rm $reverse
fi


if [ "$s1" != "0" ]  ; then exit 1; fi

