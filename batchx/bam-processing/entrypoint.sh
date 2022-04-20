#!/bin/bash
set -e
source /batchx/bx-start.sh
inputFile=/batchx/input/input.json

# Input data
refFolder=$(cat $inputFile | jq -r .refFolder)
bamFile=$(cat $inputFile | jq -r .bamFile)
gvcfCalling=$(cat $inputFile | jq -r .gvcfCalling)
panelInfo=$(cat $inputFile | jq -r .panelInfo)
sample=$(cat $inputFile | jq -r .sample)
knownSites=$(cat $inputFile | jq -r .knownSites)

printf "\n.......................\n"
printf "  PRE-PROCESSING $sample \n"
printf ".........................\n"

if [ "$knownSites" != "null" ]; then
  knownSitesField=",\"knownSites\": $knownSites"
fi
preprocessing_output=$(bx run -l -v=1 -m=4000 iis-fjd@bioinformatics/fjd/bam-preprocessing:0.0.6 \
"{
  \"refFolder\": \"$refFolder\",
  \"bamFile\": \"$bamFile\",
  \"duplicates\": \"mark\",
  \"outputPrefix\": \"$sample\"
  $knownSitesField
}")
outputBamFolder=$(echo $preprocessing_output | jq -r .outputBam)

printf "\n\n\n\n.......................\n"
printf "  SNV CALLING $sample \n"
printf ".........................\n"

if [ "$gvcfCalling" != "null" ]; then
  gvcfCallingField=",\"gvcfCalling\": $gvcfCalling"
fi
if [ "$panelInfo" != "null" ]; then
  $panelInfoField=",\"panelInfo\": $panelInfo"
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