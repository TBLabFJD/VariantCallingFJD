#!/bin/sh

sample=$1
input=$2
output=$3
ref_fasta=$4
threads=$5

export SFT=/mnt/genetica3/software/mito
alias lofreq='$SFT/lofreq_star-2.1.3.1/bin/lofreq'

bam="$input/${sample}_alignment.bam"
#bam_indelqual="$input/${sample}_alignment_indelqual.bam"
vcf="$output/${sample}_lofreq.vcf"


### when illumina libraries are available... (default parameter for illumina libraries)

#lofreq indelqual --dindel \
#  -f $ref_fasta \
#  -o $bam_indelqual \
#  $bam

#samtools index $bam_indelqual

lofreq call-parallel --pp-threads $threads \
 -f $ref_fasta \
 -o $vcf \
 $bam
# --call/indels


#rm ${bam_indelqual}*


