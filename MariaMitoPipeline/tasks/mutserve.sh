#!/bin/sh

sample=$1
bam=$2
output=$3
ref_fasta=$4

export SFT=/mnt/genetica3/software/mito
alias mutserve='java -jar $SFT/mutserve-1.3.0/mutserve-1.3.0.jar analyse-local'

mutserve --deletions True \
 --insertions True \
 --input $bam \
 --output $output/${sample}_mutserve.vcf \
 --reference ${ref_fasta} 

cp $output/${sample}_mutserve.vcf $output/..


# Filtering SNVs



#¿Eliminamos .txt?
