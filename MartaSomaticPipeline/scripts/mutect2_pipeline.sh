#!/bin/bash

# Running somatic samples from bam

input=$1
output=$2
bedfile="/home/proyectos/bioinfo/lodela/BioinfoUnit/martaLymp/3110861_Covered.bed"

vepscript="/home/proyectos/bioinfo/lodela/BioinfoUnit/martaLymp/newResults_Nov19/scripts/vep.sh"
mutectscript="/home/proyectos/bioinfo/lodela/BioinfoUnit/martaLymp/newResults_Nov19/scripts/mutect2.sh"


cd $input
echo $input

for file in  $input/*bam
do 
	echo $file
	basename=`basename $file .bam`
	sample=${basename%*_alignment}
	echo $sample
	mutectID=$(sbatch  -o ${output}/mutect2_${sample}.out -e ${output}/mutect2_${sample}.err --account=bioinfo_serv --partition=fastbioinfo $mutectscript $sample $input $output $bedfile)
	vcf_in="$output/${sample}_filtered.vcf"
	vcf_out="$output/${sample}_filteredAnnotated.vcf"
	job_id="$(cut -d " " -f4 <<< $mutectID)"
	echo $job_id
	sbatch --dependency=afterok:$job_id -o ${output}/vep_${sample}.out -e ${output}/vep_${sample}.err --account=bioinfo_serv --partition=bioinfo $vepscript $vcf_in $vcf_out vcf
	#sbatch -o ${output}/vep_${sample}.out -e ${output}/vep_${sample}.err --account=bioinfo_serv --partition=bioinfo $vepscript $vcf_in $vcf_out vcf
	#sbatch --dependency=afterok:$job_id -o ${output}/postvep_${sample}.out -e ${output}/postvep_${sample}.err --account=bioinfo_serv --partition=fastbioinfo $vepscript $vcf_in $vcf_out pvm


done



