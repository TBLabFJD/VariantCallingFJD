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
user=${20}
softwarePath=${21}

tasksPath=${softwarePath}/tasks








printf "\n.....................\n"
printf "  MERGE VCFs \n"
printf "......................\n"


# Merging filtered VCFs


bcftools merge -O z $/raw_${run}.vcf -l $GEN/my_list_of_vcfs_files_to_combine_$run.list





# Merge LOH results

for sampleFile in `$PLINKPATH/listof .hom files`
do
	sed 1d $sampleFile >> $PLINK/${run}_plink.hom
	rm $PLINK/${samplee}.*

done 







printf "\n..............................\n"
printf "  GENOTYPE REFINEMENT $sample \n"
printf "................................\n"


$tasksPath/preprocessing_BAM.sh $local $ftype $MDAP $sample $HG19















printf "\n\n\n- VARIANT ANNOTATION (VEP - ENSEMBL) "
printf "\n---------------------------------------\n"


$taskPath/VEP_annotation.sh $local $run $MDAP $sample $threads










printf "\n\n\n- OUTPUT PROCESSING "
printf "\n-------------------------\n"


$taskPath/PVM_processing.sh $local $run $MDAP $sample $threads





###############







