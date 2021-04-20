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
cvcf=$7
skipmapping=$8
REF=$9
local=${10}
intervals=${11}
removebam=${12}
padding=${13}
softwarePath=${14}
tasksPath=${softwarePath}/tasks


if [ "$skipmapping" != "True" ]; then

	FB=$MDAP/bams
	bamfile=$FB/${sample}_alignment.bam

else
	
	bamfile=${INPUT}/${sample}*.bam

fi



printf "\n..........................\n"
printf "  CALCULATE COVERAGE $sample \n"
printf "............................\n"

$tasksPath/mosdepth.sh $local $run $MDAP $sample $bamfile $panel $padding snv
s1="$?"
if [ "$s1" != "0" ]  ; then exit 1; fi



printf "\n\n\n\n.......................\n"
printf "  SNV CALLING $sample \n"
printf ".........................\n"


$tasksPath/SNVcalling.sh $local $run $MDAP $sample $REF $bamfile $intervals $panel $padding $cvcf $removebam
s1="$?"
if [ "$s1" != "0" ]  ; then exit 1; fi


if [ "$cvcf" = "True" ]; then

	rm -r $TMP
	exit 0

fi






printf "\n\n..........................\n"
printf "  VARIANT FILTERING $sample \n"
printf ".............................\n"


# CNN / HF

ftype="HF"
#ftype="CNN"

$tasksPath/SNVfiltering.sh $local $run $MDAP $sample $REF $ftype $cvcf
s1="$?"
if [ "$s1" != "0" ]  ; then exit 1; fi






