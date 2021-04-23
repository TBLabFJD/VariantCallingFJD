#!/bin/sh

############################################################
#### Pipeline  2: BAM QC  - SNV calling - SNV filtering ####
############################################################


### PIPELINE ARGUMENTS 

INPUT=$1
MDAP=$2
sample=$3
mem=$4
run=$5
panel=$6
cvcf=$7
skipmapping=$8
REF=$9
local=${10}
intervals=${11}
removebam=${12}
padding=${13}
analysis=${14}
softwarePath=${15}
tasksPath=${softwarePath}/tasks


# Definition of bam file 

if [ "$skipmapping" != "True" ]; then bamF=$MDAP/bams; else bamF=${INPUT}; fi
bamfile=${bamF}/${sample}.bam




if [ "$analysis" != "WGS" ] 

then 


	printf "\n..........................\n"
	printf "CALCULATE COVERAGE $sample \n"
	printf "..........................\n"


	$tasksPath/mosdepth.sh $local $run $MDAP $sample $bamfile $panel $padding snv
	s1="$?"
	if [ "$s1" = "2" ] ; then exit 0;  elif [ "$s1" = "1" ]; then exit 1; fi


fi


printf "\n\n.......................\n"
printf "SNV CALLING $sample \n"
printf ".........................\n"


$tasksPath/SNVcalling.sh $local $run $MDAP $sample $REF $bamfile $intervals $panel $padding $cvcf $removebam $mem
s1="$?"
if [ "$s1" != "0" ]  ; then exit 1; fi


if [ "$cvcf" = "True" ]; then

	exit 0

fi






printf "\n\n..........................\n"
printf "VARIANT FILTERING $sample \n"
printf ".............................\n"


# CNN / HF

ftype="HF"
#ftype="CNN"

$tasksPath/SNVfiltering.sh $local $run $MDAP $sample $REF $ftype $cvcf
s1="$?"
if [ "$s1" != "0" ]  ; then exit 1; fi






