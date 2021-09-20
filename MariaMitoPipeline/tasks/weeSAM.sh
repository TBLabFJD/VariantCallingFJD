#!/bin/sh

sample=$1
dir=$2

export SFT=/mnt/genetica3/software/
alias weeSAM="$SFT/weeSAM/weeSAM"

coverage_dir="${dir}/coverage"
mkdir $coverage_dir 

weeSAM --overwrite \
 --bam ${dir}/${sample}_alignment.bam \
 --out ${coverage_dir}/${sample}_coverage.txt \
 --html ${coverage_dir}/${sample}_coverage.html
