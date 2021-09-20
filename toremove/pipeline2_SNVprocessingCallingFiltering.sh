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
padding=${20}
user=${21}
mem=${22}
softwarePath=${23}
tasksPath=${softwarePath}/tasks


if [ "$skipmapping" != "True" ]; then


	printf "\n.......................\n"
	printf "  PRE-PROCESSING $sample \n"
	printf ".........................\n"

	$tasksPath/preprocessing_BAM.sh $local $run $MDAP $sample $duplicates $REF

	bamF=$MDAP/bams

else
	
	bamF=${INPUT}

fi

bamfile=${bamF}/${sample}.bam





if [ "$analysis" = "mapping" ]; then

	exit 0

fi




printf "\n\n\n\n.......................\n"
printf "  SNV CALLING $sample \n"
printf ".........................\n"


$tasksPath/SNVcalling.sh $local $run $MDAP $sample $REF $bamfile $intervals $panel $padding $cvcf $removebam $mem



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







# # REMOVE TEMPORARY DIRECTORY

# printf '\n'
# rm -r $TMP