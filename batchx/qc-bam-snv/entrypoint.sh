#!/bin/bash
set -e
source /batchx/bx-start.sh
inputFile=/batchx/input/input.json

# Input data
refFolder=$(cat $inputFile | jq -r .refFolder)
bamFile=$(cat $inputFile | jq -r .bamFile)
bamIndex=$(cat $inputFile | jq -r .bamIndex)
gvcfCalling=$(cat $inputFile | jq -r .gvcfCalling)
panelInfo=$(cat $inputFile | jq -r .panelInfo)
sample=$(cat $inputFile | jq -r .sample)
minCoverage=$(cat $inputFile | jq -r .minCoverage)

printf "\n.......................\n"
printf "  PRE-PROCESSING $sample \n"
printf ".........................\n"

if [ "$panelInfo" != "null" ]; then
  panelInfoField=",\"panelInfo\": $panelInfo"
fi
if [ "$minCoverage" != "null" ]; then
  mosdepth_output=$(bx run -l -v=1 -m=4000 iis-fjd@bioinformatics/fjd/mosdepth:0.0.3 \
  "{
    \"analysisType\": \"$analysisType\",
    \"bamFile\": \"$bamFile\",
    \"bamIndex\": \"$bamIndex\",
    \"sample\": \"$sample\",
    \"minCoverage\": \"$minCoverage\"
    $panelInfoFiled
  }")
fi
printf "\n\n\n\n.......................\n"
printf "  SNV CALLING $sample \n"
printf ".........................\n"

if [ "$gvcfCalling" != "null" ]; then
  gvcfCallingField=",\"gvcfCalling\": $gvcfCalling"
fi
snvcalling_output=$(bx run -l -v=1 -m=4000 iis-fjd@bioinformatics/fjd/snv-calling:0.0.1 \
"{
  \"refFolder\": \"$refFolder\",
  \"bamFile\": \"$outputBamFolder/$sample.bam\",
  \"bamIndex\": \"$outputBamFolder/$sample.bai\",
  \"outputPrefix\": \"$sample\"
  $gvcfCallingField
  $panelInfoField
}")
vcf=$(echo $snvcalling_output | jq -r .outputVcf)

printf "\n\n..........................\n"
printf "  VARIANT FILTERING $sample \n"
printf ".............................\n"

snvfiltering_output=$(bx run -l -v=1 -m=4000 iis-fjd@bioinformatics/fjd/snv-filtering:0.0.2 \
"{
  \"outputPrefix\": \"$sample\",
  \"refFolder\": \"$refFolder\",
  \"vcfFile\": \"$vcf\"
}")
filtered_vcf=$(echo $snvfiltering_output | jq -r .filteredVcf)

echo "{\"vcf\":\"$filtered_vcf\"}" >>/batchx/output/output.json