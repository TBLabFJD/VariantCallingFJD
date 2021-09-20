#!/bin/sh

sample=$1
bam=$2
output=$3
ref_fasta=$4

haplochecker_out="$output/haplochecker_out/"

export SFT=/mnt/genetica3/software/mito
alias mitolib='java -jar $SFT/mitolib-0.1.2/mitolib-0.1.2.jar'
alias parse_table='python3 $SFT/tasks/parse_table.py'

mitolib haplochecker --in $bam \
--ref $ref_fasta \
--out $haplochecker_out \
--MAPQ 30 \
--QUAL 20 \
--VAF 0.01

parse_table $sample $haplochecker_out 