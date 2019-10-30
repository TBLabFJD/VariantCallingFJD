#!/bin/bash 

# CNV results format change: tabulated files to VCF

dir=$1
run=$2
fasta=$3/ucsc.hg19.fasta
nprograms=0

for file in $dir/copy_number_variation_data/$run*combined*; do

	combined_file=$file
	combined_file_tmp=${combined_file}"_tmp"
	combined_file_noheader=${combined_file}"_noheader"

	outputpos_tmp=${combined_file}"_pos"
	outputfile=`basename $file .txt`


	awk -F '\t' '{if($1!~/^#/ && $14>'$nprograms'){print $5"\t"$6-1"\t"$6}}'  ${combined_file} |sed 1d | sed 's/^/chr/g'  > ${combined_file_tmp}
	awk -F '\t' '{if($1!~/^#/ && $14>'$nprograms'){print $0}}'  ${combined_file} | sed 1d > ${combined_file_noheader}

	bedtools getfasta -fi $fasta -bed $combined_file_tmp  -tab -fo $outputpos_tmp

	awk -F '\t' 'NR==FNR { a[FNR]=$2;next } {if($4=="HET"){print $5"\t"$6"\t.\t"a[FNR]"\t<"$3">\t.\tPASS\tSVTYPE="$3";END="$7"\tGT\t0/1"}else{print $5"\t"$6"\t.\t"a[FNR]"\t<"$3">\t.\tPASS\tSVTYPE="$3";END="$7"\tGT\t1/1"}}' $outputpos_tmp $combined_file_noheader | sort | uniq > $dir/${outputfile}.vcf

	rm ${combined_file_tmp}
	rm ${combined_file_noheader}
	rm ${outputpos_tmp}

done



	

