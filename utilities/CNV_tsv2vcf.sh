#!/bin/bash 
echo "hola"
echo $1

combined_file=$1
echo ${combined_file}
combined_file_tmp=${combined_file}"_tmp"
echo ${combined_file_tmp}

combined_file_noheader=${combined_file}"_noheader"
echo ${combined_file_noheader}

nprograms=$2


outputpos_tmp=${combined_file}"_pos"

fasta="/home/proyectos/bioinfo/references/hg19/ucsc.hg19.fasta"


awk -F '\t' '{if($1!~/^#/ && $14>'$nprograms'){print $5"\t"$6-1"\t"$6}}'  ${combined_file} |sed 1d | sed 's/^/chr/g'  > ${combined_file_tmp}
awk -F '\t' '{if($1!~/^#/ && $14>'$nprograms'){print $0}}'  ${combined_file} | sed 1d > ${combined_file_noheader}

bedtools getfasta -fi $fasta -bed $combined_file_tmp  -tab -fo $outputpos_tmp

awk -F '\t' 'NR==FNR { a[FNR]=$2;next } {if($4=="HET"){print $5"\t"$6"\t.\t"a[FNR]"\t<"$3">\t.\tPASS\tSVTYPE="$3";END="$7"\tGT\t0/1"}else{print $5"\t"$6"\t.\t"a[FNR]"\t<"$3">\t.\tPASS\tSVTYPE="$3";END="$7"\tGT\t1/1"}}' $outputpos_tmp $combined_file_noheader | sort | uniq > $combined_file".vcf"






	

