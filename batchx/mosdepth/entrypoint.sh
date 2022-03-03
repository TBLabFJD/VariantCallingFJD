#!/usr/bin/env bash
inputFile=/batchx/input/input.json

# Input data
bamFile=$(cat $inputFile | jq -r .bamFile)
bamIndex=$(cat $inputFile | jq -r .bamIndex)
panelFile=$(cat $inputFile | jq -r .panelInfo.panelFile)
intervalPadding=$(cat $inputFile | jq -r .panelInfo.intervalPadding)
sample=$(cat $inputFile | jq -r .sample)
analysisType=$(cat $inputFile | jq -r .analysisType)
minCoverage=$(cat $inputFile | jq -r .minCoverage)
if [ "$minCoverage" = "null" ]; then minCoverage=10; fi

## Place bam and bai in the same folder
ln -s $bamFile /tmp/$sample.bam
ln -s $bamIndex /tmp/$sample.bam.bai
## Point to link
bamFile="/tmp/$sample.bam"

# snv results folder
QC="/batchx/output/qc"
mkdir -p $QC/
cd $QC
# Ratio of covered regions
covPerc=0.9

# 1. QC to check coverage for input bams. Writing results to File

cov=0
ret=0

if [ "$panelFile" != "null" ]; then
	# run mosdepth to check QC at panel regions
	mosdepth -x --no-per-base --by $panelFile ${sample}_qc_${analysisType} ${bamFile}
	s1="$?"
	cov=$(awk -v reads="$minCoverage" '{if($1=="total" && $2==reads){print $3}}' ${sample}_qc_${analysisType}.mosdepth.region.dist.txt)
	s2="$?"
else
	# run mosdepth to check QC (whole genome)
	mosdepth -x --no-per-base ${sample}_qc_${analysisType} ${bamFile}
	ls
	s1="$?"
	cov=$(awk -v reads="$minCoverage" '{if($1=="total" && $2==reads){print $3}}' ${sample}_qc_${analysisType}.mosdepth.global.dist.txt)
	s2="$?"
fi
if [ "$s1" = "0" ]  &&  [ "$s2" = "0" ]; then
	printf '\nEXIT STATUS: 0'
	printf '\nCOVERAGE QC for '${sample}' DONE\n'
	pass=$(awk 'BEGIN{if ('$cov'>'$covPerc') print 0}')
	if [ "$pass" = "0" ]; then
		echo -e $sample"\t"$minCoverage"\t"$cov"\tPASSED\t"$bamfile >> minCovFilterResults_${analysisType}.txt
		printf 'COVERAGE LEVELS for '${sample}' ARE OK.\n'
	else
		printf 'COVERAGE LEVELS for '${sample}' ARE LOW.\n'
		echo -e $sample"\t"$minCoverage"\t"$cov"\tFAILED\t"$bamfile >> minCovFilterResults_${analysisType}.txt
		printf '\nERROR ***************SKIPPING '${analysisType}' ANALYSIS FOR SAMPLE '${sample}'.\n'
		ret=2
	fi
else
	printf "\nERROR: PROBLEMS WITH COVERAGE QC"
	ret=1
fi

# 2. If coverage is OK for sample and in case of SNVs analysis we need to compute coverage for MAF ddbb
if [ $ret = 0 ] && [ $analysisTypeType != "cnv" ]; then
	mosdepth --quantize 10: -n -x  ${sample}_cov ${bamFile}
	s1="$?"
	# in case of WES or panels we generate bed file with padding to consider only those regions
	if [ "$panel" != "null" ]; then
		# create panel with padding regions: regions where SNVs are called
		bedpadding="$QC/$(basename -- $panel)_padding_${sample}.bed"
		awk  -v pad="$padding" '{print $1"\t"$2-pad"\t"$3+pad"\t"$4}' $panel | bedtools merge -i - > $bedpadding
		s2="$?"
		# we just keep coverage at padding regions
		bedtools intersect -a ${bedpadding} -b ${sample}_cov.quantized.bed.gz -wb | awk '{print $1"\t"$2"\t"$3"\t"$NF}'  > ${sample}_padding.quantized.bed
		s3="$?"
	fi
	if [ "$s1" = "0" ]  &&  [ "$s2" = "0" ]  &&  [ "$s3" = "0" ]; then
		printf '\nEXIT STATUS: 0'
		printf '\nTARGET COVERAGE FOR '${sample}' DONE\n'
	else
		printf "\nERROR: PROBLEMS WITH TARGET COVERAGE CALCULATION"
		ret=1
	fi
fi
echo "{\"report\":\"$QC\"}" >> /batchx/output/output.json
exit $ret