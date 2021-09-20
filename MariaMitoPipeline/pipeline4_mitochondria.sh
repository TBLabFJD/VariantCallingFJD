#!/bin/sh

##################################################
#### Pipeline  2: SNV calling in mitochondria ####
##################################################


### PIPELINE ARGUMENTS 

local=$1
run=$2
MDAP=$3
sample=$4
forward=$5
reverse=$6
threads=$8
duplicates=True
tasksPath="/home/proyectos/bioinfo/lodela/VariantCallingFJD/tasks/"





####### NOPMAL

REF=$7


printf "\n...................\n"
printf "  MAPPING $sample \n"
printf ".....................\n"



$tasksPath/mapping.sh $local $run $MDAP $sample $forward $reverse $REF $threads





printf "\n.......................\n"
printf "  PRE-PROCESSING $sample \n"
printf ".........................\n"


$tasksPath/preprocessing_BAM.sh $local $run $MDAP $sample $duplicates $REF

FB=$MDAP/bams
bamfile=$FB/${sample}_alignment.bam






# ####### SHIFED

# REF=$7



# printf "\n...................\n"
# printf "  MAPPING $sample \n"
# printf ".....................\n"



# $tasksPath/mapping_BAM.sh $local $run $MDAP $sample $forward $reverse $REF $threads





# printf "\n.......................\n"
# printf "  PRE-PROCESSING $sample \n"
# printf ".........................\n"


# $tasksPath/preprocessing_BAM.sh $local $run $MDAP $sample $duplicates $REF

# FB=$MDAP/bams
# bamfile=$FB/${sample}_alignment.bam

