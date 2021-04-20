#!/bin/sh

#######################################
# TASK: Mapping and preprocessing BAM #
#######################################

sample=$1
bam=$2
output=$3
ref_fasta=$4
threads=$5
polymorphism_vcf=$6
duplicates=$7

export SFT=/mnt/genetica3/software
alias bwa='$SFT/bwa-0.7.17/bwa'
alias gatk='java -jar $SFT/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar'
alias samtools='$SFT/samtools-1.9/samtools' 

TMP=$output/tmp
mkdir $TMP

METRICS=$output/metrics
mkdir $METRICS

id=$(samtools view $bam | head -n 1 | cut -f 1-4 -d':' | sed 's/@//' | sed 's/:/_/g')

gatk SamToFastq --INPUT $bam \
 --FASTQ /dev/stdout \
 --INTERLEAVE true \
 --NON_PF true | bwa mem -v 3 -t $threads -R "@RG\tID:${id}\tSM:${sample}\tLB:${id}\tPL:ILLUMINA" -p -Y $ref_fasta /dev/stdin |  \
 samtools view -1 > $output/${sample}_mapped.bam

gatk RevertSam --INPUT $bam \
 --OUTPUT $output/${sample}_unmapped.bam

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

  dedupped_bam=$output/${sample}_dedupped.bam

  gatk MarkDuplicates \
   --INPUT $output/${sample}_mapped.bam \
   --OUTPUT $dedupped_bam \
   --METRICS_FILE $METRICS/${sample}.metrics \
   --VALIDATION_STRINGENCY SILENT \
   --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
   --ASSUME_SORT_ORDER "queryname" \
   --CLEAR_DT "false" \
   --ADD_PG_TAG_TO_READS false
else
  dedupped_bam=$output/${sample}_mapped.bam
fi

gatk SortSam \
 --INPUT $dedupped_bam \
 --OUTPUT $output/${sample}_sorted.bam  \
 --SORT_ORDER "coordinate" \
 --CREATE_INDEX false \
 --CREATE_MD5_FILE false 

gatk BaseRecalibrator --tmp-dir=$TMP \
 -R $ref_fasta \
 -I $output/${sample}_sorted.bam \
 --use-original-qualities \
 --known-sites $polymorphism_vcf \
 -O $METRICS/before_recalibrated_bqsr_data_${sample}.recal.table

gatk ApplyBQSR --tmp-dir=$TMP \
 -R $ref_fasta \
 -I $output/${sample}_sorted.bam  \
 --bqsr $METRICS/before_recalibrated_bqsr_data_${sample}.recal.table \
 -O $output/${sample}_alignment.bam \
 --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
 --add-output-sam-program-record \
 --use-original-qualities \
 --create-output-bam-index

if [ "$?" = "0" ]; then
  rm $output/${sample}_mapped.bam
  rm $output/${sample}_unmapped.bam
  rm $output/${sample}_mapped_merged.bam
  rm $output/${sample}_dedupped.bam
  rm $output/${sample}_sorted.bam
  rm -r $TMP
  rm -r $METRICS
fi
