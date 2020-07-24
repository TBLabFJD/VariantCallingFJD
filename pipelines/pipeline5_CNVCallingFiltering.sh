#!/bin/sh

####################################################################
#### Pipeline  2: Preprocessing - SNV calling - VEP annotation ####
####################################################################


### PIPELINE ARGUMENTS 

# Arguments
MDAP=$1
HCGVCFD=$2
CNV="${MDAP}/copy_number_variation_data"
SAMPLEFILE=$3
run=$4
threads=$5
panel=$6
window=$7
utilitiesPath=$8
local=$9
methods=${10}
depth=${11}
genefilter=${12}
genome=${13}
softwarePath=${14}
tasksPath=${softwarePath}/tasks



if echo "$methods" | grep -q "QC"; then


	printf "\n.......................\n"
	printf "  QUALITY-CONTROL  \n"
	printf ".........................\n"


	# sort, remove small regions and create new window size
	printf "\n\nSorting and removing small targets in bed file\n"

	panel_10bp=${MDAP}/$(basename "$panel" .bed)_${run}_10bp.bed

	awk '{if(($3-$2)>10 && $1!="chrX" && $1!="chrY" && $1!="Y" && $1!="X"){print $0}}' $panel | sort -V -k1,1 -k2,2 > $panel_10bp
	#awk '{if(($3-$2)>10){print $0}}' $panel | sort -V -k1,1 -k2,2 > $panel_10bp
	
	if [ "$window" != "False" ]; then
		printf "\n\nGenerating sliding windows across the BED file"
		panel_out=${MDAP}/$(basename "$panel_10bp" .bed)_125window.bed
		python $utilitiesPath/CNV_windowSize.py $panel_10bp 125 ${panel_out}_unsorted 
		sort -V -k1,1 -k2,2 ${panel_out}_unsorted | uniq > $panel_out
		rm ${panel_out}_unsorted 
	fi
	
fi










