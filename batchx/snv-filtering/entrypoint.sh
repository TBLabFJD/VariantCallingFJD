#!/usr/bin/env bash
set -e
inputFile=/batchx/input/input.json

# Input data
vcfFile=$(cat $inputFile | jq -r .vcfFile)
refFolder=$(cat $inputFile | jq -r .refFolder)
sample=$(cat $inputFile | jq -r .sample)

refFile=$(ls $refFolder/*.fa.gz 2>/dev/null)
if [ -z "$refFile" ]; then
  refFile=$(ls $refFolder/*.fa 2>/dev/null)
fi

mkdir -p /batchx/output/
mkdir -p /tmp/tmp/

printf "\n\n\n- Hard Filtering (classical way GATK)"
printf "\n---------------------------------------\n"
#HARD FILTERING
#First step extacting the SNP's
#Second step extracting the INDEL's
### SNPs

#1.Extract the SNP's from the call set.

gatk SelectVariants --tmp-dir /tmp/tmp \
  -R $refFile \
  -V $vcfFile \
  --select-type-to-include SNP \
  -O /batchx/output/${sample}.snp.vcf

rm -rf /tmp/tmp/*
#2.Apply the filters to the SNP's callset.

gatk VariantFiltration --tmp-dir /tmp/tmp \
  -V /batchx/output/${sample}.snp.vcf \
  -filter "QD < 2.0" --filter-name "QD2" \
  -filter "QUAL < 30.0" --filter-name "QUAL30" \
  -filter "SOR > 3.0" --filter-name "SOR3" \
  -filter "FS > 60.0" --filter-name "FS60" \
  -filter "MQ < 40.0" --filter-name "MQ40" \
  -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
  -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
  -O /batchx/output/${sample}.filtered.snp.vcf

rm -rf /tmp/tmp/*

### INDELs

#3. Extract the INDELS from the ORIGINAL call set.
gatk SelectVariants --tmp-dir /tmp/tmp \
  -R $refFile \
  -V $vcfFile \
  --select-type-to-include INDEL \
  -O /batchx/output/${sample}.indel.vcf

rm -rf /tmp/tmp/*

#4.Apply the filters to the INDEL's callset.
gatk VariantFiltration --tmp-dir /tmp/tmp \
  -V /batchx/output/${sample}.indel.vcf \
  -filter "QD < 2.0" --filter-name "QD2" \
  -filter "QUAL < 30.0" --filter-name "QUAL30" \
  -filter "FS > 200.0" --filter-name "FS200" \
  -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
  -O /batchx/output/${sample}.filtered.indel.vcf

rm -rf /tmp/tmp/*
# Combine Variants after using SNPS and INDELS filtering into a single file and get it ready for annotation.

vcf=/batchx/output/${sample}.gatkLabeled.vcf
gatk MergeVcfs --TMP_DIR /tmp/tmp \
  -R $refFile \
  -I /batchx/output/${sample}.filtered.snp.vcf \
  -I /batchx/output/${sample}.filtered.indel.vcf \
  -O $vcf

echo "{\"filteredVcf\":\"$vcf\"}" >>/batchx/output/output.json
