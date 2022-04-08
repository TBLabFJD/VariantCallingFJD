#!/usr/bin/env bash
set -e
inputFile=/batchx/input/input.json

# Input data
refFolder=$(cat $inputFile | jq -r .refFolder)
bamFile=$(cat $inputFile | jq -r .bamFile)
bamIndex=$(cat $inputFile | jq -r .bamIndex)
sample=$(cat $inputFile | jq -r .sample)
panelFile=$(cat $inputFile | jq -r .panelInfo.panelFile)
intervalPadding=$(cat $inputFile | jq -r .panelInfo.intervalPadding)
gvcfCalling=$(cat $inputFile | jq -r .panelInfo.gvcfCalling)
if [ "$intervalPadding" = "null" ]; then
  intervalPadding=0
fi
refFile=$(ls $refFolder/*.fa.gz 2>/dev/null)
if [ -z "$refFile" ]; then
  refFile=$(ls $refFolder/*.fa 2>/dev/null)
fi

mkdir -p /batchx/output/
mkdir -p /tmp/tmp/
ln -s $bamFile /tmp/$(basename -- "$bamFile")
ln -s $bamIndex /tmp/$(basename -- "$bamIndex")
bamFile=/tmp/$(basename -- "$bamFile")

printf "\n\n\n- HAPLOTYPECALLER (GATK)"
printf '\nGATK HaplotypeCallerGVCF for '${sample}' STARTS'
outputBam=/batchx/output/${sample}.bamout.bam
if [ "$gvcfCalling" = "true" ]; then
  if [ "$panelFile" = "null" ]; then
    outputVcf=/batchx/output/${sample}.g.vcf
    gatk HaplotypeCaller --tmp-dir /tmp/tmp \
      -R $refFile \
      -I $bamFile \
      -ERC GVCF \
      -bamout $outputBam \
      -O $outputVcf \
      -G StandardAnnotation \
      -G AS_StandardAnnotation \
      -G StandardHCAnnotation \
      -A FisherStrand -A StrandOddsRatio -A RMSMappingQuality -A MappingQualityRankSumTest -A ReadPosRankSumTest -A DepthPerSampleHC -A BaseQualityRankSumTest -A ExcessHet \
      --annotate-with-num-discovered-alleles true
  else
    gatk HaplotypeCaller --tmp-dir /tmp/tmp \
      -R $refFile \
      -I $bamFile \
      -ERC GVCF \
      -bamout $outputBam \
      -O $outputVcf \
      -G StandardAnnotation \
      -G AS_StandardAnnotation \
      -G StandardHCAnnotation \
      -A FisherStrand -A StrandOddsRatio -A RMSMappingQuality -A MappingQualityRankSumTest -A ReadPosRankSumTest -A DepthPerSampleHC -A BaseQualityRankSumTest -A ExcessHet \
      --annotate-with-num-discovered-alleles true \
      -L $panelFile -ip $intervalPadding
  fi
  printf '\nGATK HaplotypeCaller in GVCF mode for '${sample}' DONE\n'
else
  outputVcf=/batchx/output/${sample}.vcf
  if [ "$panelFile" = "null" ]; then
    gatk HaplotypeCaller --tmp-dir /tmp/tmp \
      -R $refFile \
      -I $bamFile \
      -bamout $outputBam \
      -O $outputVcf \
      --annotate-with-num-discovered-alleles true
  else
    gatk HaplotypeCaller --tmp-dir /tmp/tmp \
      -R $refFile \
      -I $bamFile \
      -bamout $outputBam \
      -O $outputVcf \
      --annotate-with-num-discovered-alleles true \
      -L $panelFile -ip $intervalPadding
  fi
  printf '\nGATK HaplotypeCaller in VCF mode for '${sample}' DONE\n'
fi
echo "{\"outputVcf\":\"$outputVcf\",\"outputBam\":\"$outputBam\"}" >>/batchx/output/output.json
