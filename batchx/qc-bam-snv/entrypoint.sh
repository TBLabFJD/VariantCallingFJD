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
  vcpus=$(cat $inputFile | jq -r .resources.mosdepth.vcpus)
  if [ "$vcpus" == "null" ]; then
    vcpus="1"
  fi
  memory=$(cat $inputFile | jq -r .resources.mosdepth.memory)
  if [ "$memory" == "null" ]; then
    memory="4000"
  fi
  mosdepth_output=$(bx run -l -v=$vcpus -m=$memory iis-fjd@mosdepth:0.0.1 \
  "{
    \"analysisType\": \"snv\",
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
vcpus=$(cat $inputFile | jq -r .resources.snvCalling.vcpus)
if [ "$vcpus" == "null" ]; then
  vcpus="1"
fi
memory=$(cat $inputFile | jq -r .resources.snvCalling.memory)
if [ "$memory" == "null" ]; then
  memory="4000"
fi
snvcalling_output=$(bx run -l -v=$vcpus -m=$memory iis-fjd@snv-calling:0.0.2 \
"{
  \"refFolder\": \"$refFolder\",
  \"bamFile\": \"$bamFile\",
  \"bamIndex\": \"$bamIndexXX\",
  \"sample\": \"$sample\"
  $gvcfCallingField
  $panelInfoField
}")
vcf=$(echo $snvcalling_output | jq -r .outputVcf)

printf "\n\n..........................\n"
printf "  VARIANT FILTERING $sample \n"
printf ".............................\n"
vcpus=$(cat $inputFile | jq -r .resources.snvFiltering.vcpus)
if [ "$vcpus" == "null" ]; then
  vcpus="1"
fi
memory=$(cat $inputFile | jq -r .resources.snvFiltering.memory)
if [ "$memory" == "null" ]; then
  memory="4000"
fi
snvfiltering_output=$(bx run -l -v=$vcpus -m=$memory iis-fjd@snv-filtering:0.0.1 \
"{
  \"sample\": \"$sample\",
  \"refFolder\": \"$refFolderXX\",
  \"vcfFile\": \"$vcf\"
}")
filtered_vcf=$(echo $snvfiltering_output | jq -r .filteredVcf)
echo "{\"vcf\":\"$filtered_vcf\"}" >>/batchx/output/output.json
