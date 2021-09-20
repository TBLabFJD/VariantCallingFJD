#!/bin/sh

#####################################
#### SNV calling in mitochondria ####
#####################################

### PIPELINE ARGUMENTS 

threads=$1
sample=$2
input_type=$3
input_dir=$4
output_dir=$5
ref_fasta=$6
duplicates=$7
polymorphism_vcf=$8

tasksPath=""

### CREATE OUTPUT DIR FOR SAMPLE

output=${output_dir}/$sample
mkdir $output

mapping_dir=${output}/alignment
mkdir $mapping_dir

### ALIGNMENT

printf "\n...................\n"
printf "  MAPPING $sample \n"
printf ".....................\n"

case $input_type in
	fastq )
		forward=${input_dir}/${sample}*1*f*q
		reverse=${input_dir}/${sample}*2*f*q
		$tasksPath/fastqtobam.sh $sample $forward $reverse $mapping_dir $ref_fasta $threads $polymorphism_vcf $duplicates
		;;
	bam )
		bam=${input_dir}/${sample}*bam
		$tasksPath/bamtobam.sh $sample $bam $mapping_dir $ref_fasta $threads $polymorphism_vcf $duplicates
		;;
esac

alignment_bam="${mapping_dir}/${sample}_alignment.bam"



### GET CONTAMINATION
# Runs haplochecker to study major and minor haplogroups

$tasksPath/get_contamination.sh $sample $alignment_bam $output $ref_fasta

printf "\n.......................\n"
printf "  CALLING SNV (SAMPLE: $sample) \n"
printf ".........................\n"

$tasksPath/lofreq.sh $sample $mapping_dir $output $ref_fasta $threads

printf "\n.......................\n"
printf "  CALLING INDEL (SAMPLE: $sample) \n"
printf ".........................\n"

MS_output="$output/mutserve"
mkdir $MS_output

$tasksPath/mutserve.sh $sample $alignment_bam $MS_output $ref_fasta

#Eliminar snvs


#Unir ambos vcf

