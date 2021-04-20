#!/bin/bash

# Running a somatic samples from bam

input=$1  # folder where bam files are stored
output=$2 # ooutput folder
file=$3

bedfile="/home/proyectos/bioinfo/lodela/BioinfoUnit/martaLymp/3110861_Covered.bed"

mutectscript="/home/proyectos/bioinfo/mrod/.....BioinfoUnit/martaLymp/newResults_Nov19/scripts/mutect2.sh"


echo $file
basename=`basename $file .bam`
sample=${basename%*_alignment}
echo $sample
sbatch  -o ${output}/mutect2_${sample}.out -e ${output}/mutect2_${sample}.err --account=bioinfo_serv --partition=fastbioinfo $mutectscript $sample $input $output $bedfile


