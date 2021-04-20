#!/bin/sh

sample=$1
forward=$2
reverse=$3
output=$4
ref_fasta=$5
threads=$6
polymorphism_vcf=$7
duplicates=$8

export SFT=/mnt/genetica3/software
alias bwa='$SFT/bwa-0.7.17/bwa'
alias gatk='java -jar $SFT/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar'
alias samtools='$SFT/samtools-1.9/samtools' 


TMP=$output/tmp
mkdir $TMP

METRICS=$output/metrics
mkdir $METRICS

bwa mem -v 3 -t $threads -Y \
$ref_fasta \
${forward} \
${reverse} |  samtools view -1 > $output/${sample}_mapped.bam

gatk FastqToSam -TMP_DIR $TMP \
--FASTQ ${forward} \
--FASTQ2 ${reverse} \
--OUTPUT $output/${sample}_unmapped.bam \
--READ_GROUP_NAME ${sample} \
--SAMPLE_NAME ${sample} \
--PLATFORM illumina 

gatk MergeBamAlignment -TMP_DIR $TMP \
  --VALIDATION_STRINGENCY SILENT \
  --EXPECTED_ORIENTATIONS FR \
  --ATTRIBUTES_TO_RETAIN X0 \
  --ATTRIBUTES_TO_REMOVE NM \
  --ATTRIBUTES_TO_REMOVE MD \
  --ALIGNED_BAM $output/${sample}_mapped.bam \
  --UNMAPPED_BAM $output/${sample}_unmapped.bam  \
  --OUTPUT $output/${sample}_mapped_merged.bam \
  --REFERENCE_SEQUENCE $ref_fasta \
  --PAIRED_RUN true \
  --SORT_ORDER "unsorted" \
  --IS_BISULFITE_SEQUENCE false \
  --MAX_RECORDS_IN_RAM 2000000 \
  --ALIGNED_READS_ONLY false \
  --CLIP_ADAPTERS false \
  --ADD_MATE_CIGAR true \
  --MAX_INSERTIONS_OR_DELETIONS -1 \
  --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
  --UNMAPPED_READ_STRATEGY COPY_TO_TAG \
  --ALIGNER_PROPER_PAIR_FLAGS true \
  --UNMAP_CONTAMINANT_READS true \
  --ADD_PG_TAG_TO_READS false

if [ "$duplicates" != "True" ]; then
  dedupped=$output/${sample}_dedupped.bam
  gatk MarkDuplicates \
  -I $output/${sample}_mapped_merged.bam \
  -O $dedupped \
  -M $METRICS/${sample}_marked_dup_metrics.txt \
  --REMOVE_DUPLICATES false --ASSUME_SORT_ORDER queryname \
  --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 --VALIDATION_STRINGENCY SILENT
else

  dedupped=$output/${sample}_mapped_merged.bam
fi

gatk SortSam \
  --INPUT $dedupped \
  --OUTPUT /dev/stdout \
  --SORT_ORDER "coordinate" \
  --CREATE_INDEX false \
  --CREATE_MD5_FILE false \
| gatk  SetNmMdAndUqTags \
  --INPUT /dev/stdin \
  --OUTPUT $output/${sample}_sorted.bam \
  --CREATE_INDEX true \
  --REFERENCE_SEQUENCE $ref_fasta

gatk BaseRecalibrator --tmp-dir $TMP \
 -R $ref_fasta \
 -I  $output/${sample}_sorted.bam\
 --use-original-qualities \
 --known-sites $polymorphism_vcf \
 -O $METRICS/before_recalibrated_bqsr_data_${sample}.recal.table

gatk ApplyBQSR --tmp-dir $TMP \
 -R $ref_fasta \
 -I $output/${sample}_sorted.bam  \
 --bqsr $METRICS/before_recalibrated_bqsr_data_${sample}.recal.table \
 -O $output/${sample}_alignment.bam \
 --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
 --add-output-sam-program-record \
 --create-output-bam-md5 \
 --use-original-qualities \
 --create-output-bam-index

if [ "$?" = "0" ]; then
  rm $output/${sample}_dedupped.bam
  rm $output/${sample}_mapped.bam
  rm $output/${sample}_unmapped.bam 
  rm $output/${sample}_mapped_merged.bam
  rm $output/${sample}_sorted.bam
  rm -r $TMP
  rm -r $METRICS
fi


