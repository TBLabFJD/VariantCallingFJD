#!/bin/sh

#######################################################
### FJD pipeline - Processing and Haplotype Calling ###
#######################################################

### FJD PIPELINE ARGUMENTS:

INPUT=$1
MDAP=$2
sample=$3
threads=$4
run=$5
panel=$6
basespace=$7
cat=$8
fastqFolder=$9
analysis=${10}
cvcf=${11}
skipmapping=${12}    
HG19=${13}
local=${14}
pathology=${15}
intervals=${16}
duplicates=${17}
removebam=${18}
genefilter=${19}
user=${20}
utilitiesPath=${21}

printf "\n......................................................\n"
printf "  MAPPING OR/AND SNV CALLING FOR SAMPLE $sample \n"
printf "......................................................\n"


# REQUIRED MODULES AND DATABASES:

if [ "$local" != "True" ]; then

#	module load miniconda/2.7
	module load bwa/0.7.15
	module load samtools/1.9
	module load picard/2.18.9
	module load gatk/4.1.2.0
	module load vep/release98
	module load bedtools/2.27.0
	module load R
	alias plink='/usr/local/bioinfo/plink/plink'
	alias picard='java -jar /usr/local/bioinfo/picard-tools/2.18.9/picard.jar'
	alias gatk='java -jar /usr/local/bioinfo/gatk/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar'

	VEP="/usr/local/bioinfo/vep/ensembl-vep/vep"
	FILTER_VEP='/usr/local/bioinfo/vep/ensembl-vep/filter_vep'
	VEP_CACHE='/usr/local/bioinfo/vep/ensembl-vep/t/testdata/cache/homo_sapiens'
	VEP_FASTA="/home/proyectos/bioinfo/references/VEPfasta/Homo_sapiens.GRCh37.dna.primary_assembly.fa"
	PLUGIN_DIR=/usr/local/bioinfo/vep/ensembl-vep/Plugins
	PLUGIN_DBS="/home/proyectos/bioinfo/references/VEPdbs"
	dbNSFP_DB="${PLUGIN_DBS}/dbNSFP3.5a_hg19.gz"
	CCS_DB="/home/proyectos/bioinfo/references/CCS/ccrs.autosomes.v2.20180420.bed.gz"


	softwareFile="${MDAP}/software_${run}.txt"
	title="MAPPING"
	if [ ! -f $softwareFile ] || [ `grep -q $title $softwareFile` ] ; then 

		printf "MAPPING / SINGLE NUCLEOTIDE VARIANT CALLING:\n" >> ${softwareFile}
		module list 2>> ${softwareFile}
	
	fi




else

	export SFT=/mnt/genetica3/marius/pipeline_practicas_marius/software
	alias plink='$SFT/plink/plink'
	alias bwa='$SFT/bwa/bwa'
	alias samtools='$SFT/samtools/samtools'
	alias picard='java -Xmx10g -jar $SFT/picard/build/libs/picard.jar'
	alias gatk='java  -Xmx10g -jar $SFT/gatk/build/libs/gatk-package-4.0.6.0-22-g9d9484f-SNAPSHOT-local.jar'
	VEP="${SFT}/variant_effect_predictor/ensembl-vep/vep"
	FILTER_VEP="${SFT}/variant_effect_predictor/ensembl-vep/filter_vep"
	VEP_CACHE='/mnt/genetica3/marius/pipeline_practicas_marius/software/variant_effect_predictor/.vep'
	VEP_FASTA="${VEP_CACHE}/homo_sapiens/93_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa"
	PLUGIN_DIR="${VEP_CACHE}/Plugins"
	PLUGIN_DBS="${VEP_CACHE}/dbs"
	dbNSFP_DB="${PLUGIN_DBS}/dbNSFP_hg19.gz"
	#CCS_DB="/home/proyectos/bioinfo/references/CCS"

	softwareFile="${MDAP}/software_${run}.txt"
	title="MAPPING"
	
	if [ ! -f $softwareFile ] || [ `grep -q $title $softwareFile` ] ; then 

		printf "MAPPING / SINGLE NUCLEOTIDE VARIANT CALLING:\n\n" >> ${softwareFile}
		
		printf "BWA VERSION\n" >> ${softwareFile}

		bwa 2>&1 | head -n4 | tail -n2 | head -n1 >> ${softwareFile}
		
		printf "\nSamtools VERSION\n" >> ${softwareFile}
		samtools --version 2>&1 | head -n2 >> ${softwareFile}
		
		printf "\nPicard VERSION\n" >> ${softwareFile}
		picard MarkDuplicates --version  2>&1 | head -n2 >> ${softwareFile}
		
		printf "\nGATK VERSION\n" >> ${softwareFile}
		gatk ApplyBQSR 2>&1 | head -n4 | tail -n1 >> ${softwareFile}

		printf "\nVEP VERSION\n" >> ${softwareFile}
		$VEP --help | grep "Versions:" -A 5 | tail -n4 >> ${softwareFile}

	fi




fi




### FJD PIPELINE VARIABLES:
#··························


#Mapped_data
MD="${MDAP}/mapped_data"

#Sorted_data
SD="${MDAP}/sorted_data"
	
#Dedupped_data
DD="${MDAP}/dedupped_data"

#Recalibrated_base_quality_scores_data: generates a recalibration table  based on specified covariates, read_group,reported_quiality_score,machine_cycle and nucleotide_context, after GATK.
RBQSRD="${MDAP}/recalibrated_bqsr_data"

#First using the APPLIED_RECALIBRATION_BQSR_DATA BAM file, running again BaseRecalibrator and generating the plots with AnalyzeCovariates.
PRD="${MDAP}/metrics"


#Haplotype_caller_gvcf_data:calling for SNPs and indels via local re-assembly of HAPLOTYPES using HAPLOTYPECALLER by GATK.
#Perform joint genotyping on one or more samples precalled with Haplotype_caller, if one sample,
#straight after haplotypecaller and if more than one use Combine_gvcfs.
HCGVCFD="${MDAP}/haplotype_caller_gvcf_data"


#Genotyped_vcf_data: single or single-multisample GVCF as input, output will be a VCF.
GVCFD="${MDAP}/genotyped_vcf_data"

#HARD_FILTERING, filters for SNPs and filters for INDELs.
#Variant_filtration_vcf_data:hard filtering process to select based on the INFO and FORMAT annotations(QD,MQO,FS,MQ,ReadPosrankSum)
VFVCFD="${MDAP}/snv_results"

# Hard filtering, vep filtering and vep annotation
VEPVCFAD="${MDAP}/snv_results"

# LOH regions
PLINK="${MDAP}/plink_results"

TMP=$MDAP/tmp_${sample}
mkdir $TMP





### INPUT DATA RECOVERY:
#·······················


if [ "$skipmapping" != "True" ]; then

	if [ "$cat" = "True" ] || [ "$basespace" = "True" ]; then

		python2 $utilitiesPath/bscpCat_sample.py $basespace $fastqFolder $INPUT $sample $user
		if [ "$?" != "0" ]; then
			printf "\nERROR: PROBLEMS WITH BASESPACE DATA DOWNLOADING/CONCATENATION"
			exit 1
		else
			printf '\nEXIT STATUS: 0'
			printf '\nDOWNLOAD/CAT FOR '${sample}'  DONE'
		fi

		INPUT=$fastqFolder

		if [[ "$basespace" = "True" ]]; then
			sample=${sample:0:7}
		fi
	fi

	foward="${INPUT}/${sample}*_R1*.f*q.gz"
	reverse="${INPUT}/${sample}*_R2*.f*q.gz"

	ABQSRD="${MDAP}/bams"
	bqsr_bamfile=${ABQSRD}/${sample}_alignment.bam
	

else
	
	bqsr_bamfile=${INPUT}/${sample}*.bam
	echo $bqsr_bamfile

fi





### ANALYSIS PIPELINE:
#····················..

if [ "$skipmapping" != "True" ]; then



	printf "\n\n\n- INDEXING REFERENCE FILES (BWA) "
	printf "\n----------------------------------\n"

	# BWA index
	if [ ! -f $HG19/ucsc.hg19.fasta.sa ]; then
		printf '\nRUN BWA INDEX!!!'
		exit 0
		#printf 'Starts BWA INDEX'
		#bwa index $HG19/ucsc.hg19.fasta
		#printf sampleFile '\nBWA INDEX DONE' 
	else
		printf '\nBWA INDEX already existing. Reference indexing for BWA skipped'
	fi




	# Genome index
	if [ ! -f $HG19/ucsc.hg19.fasta.fai ]; then
		printf '\nRUN REFERNCE INDEX!!!'
		exit 0		
		# printf '\nCreate .FAI file, using samtools faidx'
		# samtools faidx $HG19/ucsc.hg19.fasta 
		# printf sampleFile '\nucsc.hg19.fasta.FAI DONE'
	else
		printf '\nREFERENCE INDEX already existing. Fasta indexing skipped'


	fi



	# Genome dict
	if [ ! -f $HG19/ucsc.hg19.dict ]; then
		printf '\nCREATE REFERENCE DICT!!!'
		exit 0	
		#Creating .DICT in HG19.
		# printf 'Create .DICT file, using picardtools CreateSequnceDictionary'
		# picard CreateSequenceDictionary R=$HG19/ucsc.hg19.fasta O=$HG19/ucsc.hg19.2.dict
		# printf sampleFile 'ucsc.hg19.DICT DONE'
	else
		printf '\nREFERENCE DICT already existing. Dictionary generation skipped'

	fi







	printf "\n\n\n- MAPPING READS (BWA) "
	printf "\n----------------------------\n"

	#mapping fastq files to reference_genome after BWA INDEX
	#reference genome is: UCSC.HG19.FASTA

	mkdir $MD
	printf '\nmkdir mapped_data' 


	printf '\nBuilding the header for '${sample}' ongoing...\n'
	#build header
	start=`date +%s`

	header=$(zcat $foward | head -n 1)
	id=$(echo  $header | head -n 1 | cut -f 1-4 -d':' | sed 's/@//' | sed 's/:/_/g')
	#sm=$(echo  $header | head -n 1 | grep -Eo '[ATGCN]+$')
	echo  -e "\nThis is how the new header looks\n"
	echo  -e '@RG\tID:'$id'\tSM:'${sample}'\tPL:ILLUMINA'

	echo -e 'bwa mem -v 3 -t '$threads' -R @RG\tID:'$id'\tSM:'${sample}'\tPL:ILLUMINA '$HG19'/ucsc.hg19.fasta'$foward$reverse'>'$MD'/mapped_'${sample}'.sam'


	bwa mem -v 3 -t $threads -R '@RG\tID:'$id'\tSM:'${sample}'\tPL:ILLUMINA'\
	$HG19/ucsc.hg19.fasta \
	$foward \
	$reverse > $MD/mapped_${sample}.sam

	if [ "$?" = "0" ]; then
		printf '\nEXIT STATUS: 0'
		printf  '\nBWA MEM '${sample}'  DONE'
		if [ "$cat" = "True" ] || [ "$basespace" = "True" ]; then
			rm $foward
			rm $reverse
		fi
	else
		printf "\nERROR: PROBLEMS WITH MAPPING"
		exit 1
	fi

	end=`date +%s`
	runtime=$((end-start))
	printf '\nExecuting time: '$runtime 



	printf "\n\n\n- SORTING BAM (PICARD) "
	printf "\n------------------------\n"
	mkdir $SD
	printf '\nmkdir sorted_data'
	#Sorting mapped data.
	printf '\nRun picard SortSam for '${sample}''
	start=`date +%s`

	picard SortSam TMP_DIR=$TMP I=$MD/mapped_${sample}.sam  \
	         O=$SD/sorted${sample}.bam \
	         SO=coordinate


	if [ "$?" = "0" ]; then
		printf '\nEXIT STATUS: 0'
		printf '\nPicard SortSam for '${sample}'  DONE' 
		rm $MD/mapped_${sample}.sam

	else
		printf "\nERROR: PROBLEMS WITH SORTING"
		exit 1
	fi


	end=`date +%s`
	runtime=$((end-start))
	printf '\nExecuting time: '$runtime 


	mkdir $PRD
	printf '\nmkdir metrics' 



	if [ "$duplicates" != "True" ]; then


		printf "\n\n\n- MARKING DUPLICATES (PICARD)"
		printf "\n--------------------------------\n"

		#Selecting the duplicates reads from the mapped and sorted reads.
		mkdir $DD
		printf '\nmkdir dedupped_data'

		start=`date +%s`

		#Mark duplicates PICARD

		bambeforebqsr=$DD/dedupped_${sample}.bam
		baibeforebqsr=$DD/dedupped_${sample}.bai

		printf '\nStart picard MarkDuplicates '${sample}''
		picard MarkDuplicates TMP_DIR=$TMP \
		 I=$SD/sorted${sample}.bam \
		 O=$bambeforebqsr \
		 M=$PRD/marked_dup_metrics_${sample}.txt \
		 REMOVE_DUPLICATES=false \
		 AS=SortOrder



		if [ "$?" = "0" ]; then
			printf '\nEXIT STATUS: 0'
			printf  '\nPICARD MarkDuplicates '${sample}' DONE' 
			rm $SD/sorted${sample}.bam


		else
			printf "\nERROR: PROBLEMS WITH MARKING DUPLICATES"
			exit 1
		fi


		end=`date +%s`
		runtime=$((end-start))
		printf '\nExecuting time: '$runtime 


	else

		bambeforebqsr=$SD/sorted${sample}.bam 
		baibeforebqsr=$SD/sorted${sample}.bai


	fi




	printf "\n\n\n- Dedupped BAM INDEX (PICARD)"
	printf "\n-----------------------------------\n"

	#	Create a .BAI file to compare original vs removed duplicates reads.
	#	#Indexing the BAM files.
	#	#Generating the .BAI files from the DEDUPPED (markedDuplicates from the original SAM/BAM file).
	printf '\nIndexing '${sample}' BAM files'
	start=`date +%s`

	picard BuildBamIndex TMP_DIR=$TMP \
	 I=$bambeforebqsr \
	 O=$baibeforebqsr


	if [ "$?" = "0" ]; then
		printf '\nEXIT STATUS: 0'
		printf  '\nPICARD BuildBamIndex '${sample}' DONE'

	else
		printf "\nERROR: PROBLEMS WITH THE BUILDING OF THE BAM INDEX"
		exit 1
	fi


	end=`date +%s`
	runtime=$((end-start))
	printf '\nExecuting time: '$runtime 









	printf "\n\n\n- BQSR RECALIBRATION TABLE (GATK)"
	printf "\n------------------------------------\n"


	#Recalibrating  the reads using base quality score reads.
	mkdir $RBQSRD
	printf '\nmkdir recalibrated_bqsr_data' 

	#GATK BaseRecalibration first table
	#BaseRecalibration + table
	printf '\nStarts GATK '${sample}' Recalibrator'
	start=`date +%s`

	gatk BaseRecalibrator --tmp-dir=$TMP \
	-R $HG19/ucsc.hg19.fasta \
	-I $bambeforebqsr \
	--known-sites $HG19/dbsnp_138.hg19.vcf \
	--known-sites $HG19/1000G_phase1.indels.hg19.sites.vcf \
	--known-sites $HG19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
	-O $RBQSRD/before_recalibrated_bqsr_data_${sample}.recal.table
	#--bqsr 1st_racalibration.table


	if [ "$?" = "0" ]; then
		printf '\nEXIT STATUS: 0'
		printf  '\nGATK BaseRecalibrator '${sample}' DONE' 

	else
		printf "\nERROR: PROBLEMS WITH BQSR RECALIBRATION TABLE"
		exit 1
	fi


	end=`date +%s`
	runtime=$((end-start))
	printf '\nExecuting time: '$runtime 










	printf "\n\n\n- APPLYING BQSR  (GATK) "
	printf "\n-------------------------\n"


	#Applying the recalibration table to the bam file to continue the analysis.
	#ApplyBQSR

	printf '\nStarts picard  '${sample}' ApplyBQSR'
	start=`date +%s`
	printf '\nmkdir bams' 
	mkdir $ABQSRD

	gatk ApplyBQSR --tmp-dir=$TMP \
	-R $HG19/ucsc.hg19.fasta \
	-I $bambeforebqsr \
	--bqsr $RBQSRD/before_recalibrated_bqsr_data_${sample}.recal.table \
	-O $bqsr_bamfile

	if [ "$?" = "0" ]; then
		printf '\nEXIT STATUS: 0'
		printf  '\nGATK ApplyBQSR '${sample}' DONE'
		rm $bambeforebqsr
		rm $baibeforebqsr
	else
		printf "\nERROR: PROBLEMS WITH BQSR"
		exit 1
	fi

	end=`date +%s`
	runtime=$((end-start))
	printf '\nExecuting time: '$runtime 







	printf "\n\n\n- ANALYZE COVARIATES FOR PDF PLOTS. COMPARISON OF THE RECALIBRATION DATA (GATK)"
	printf "\n---------------------------------------------------------------------------------\n"


	printf '\nStarts GATK Second Recalibration\n'

	#GATK BaseRecalibration second table for next step AnalyzeCovariates.
	#Generates the second pass table.
	#instead of second table, we use the BAM created by ApplyBQSR to regenrate a new TABLE for plot.
	start=`date +%s`

	gatk BaseRecalibrator --tmp-dir=$TMP \
	-I $bqsr_bamfile \
	-R $HG19/ucsc.hg19.fasta \
	--known-sites $HG19/dbsnp_138.hg19.vcf \
	--known-sites $HG19/1000G_phase1.indels.hg19.sites.vcf \
	--known-sites $HG19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
	-O $PRD/after_recalibrated_bqsr_data_${sample}.recal.table
	printf '\nSecond recalibration GATK table --> DONE'

	#Full generating the Plots of the recalibration tables using AnalyzeCovariates and saving a csv copy.
	#Analyze  the tables.

	printf   '\nGenerating Plots, pdf and csv files...'
	gatk AnalyzeCovariates --tmp-dir=$TMP \
	-before $RBQSRD/before_recalibrated_bqsr_data_${sample}.recal.table \
	-after $PRD/after_recalibrated_bqsr_data_${sample}.recal.table \
	-csv $PRD/BQSR_${sample}.csv \
	-plots $PRD/AnalyzeCovariates_bqsr_${sample}.pdf

	#Obtaining an CSV and PDF file of the comparrisson between first and second pass of the recalibration applied to the bam

	if [ "$?" = "0" ]; then
		rm $RBQSRD/before_recalibrated_bqsr_data_${sample}.recal.table
		printf '\nEXIT STATUS: 0'
		printf '\nPlot and CSV file for '${sample}' is DONE!'
		printf '\n Plots files generated --> DONE'

	else
		printf "\nERROR: PROBLEMS ANALYZING COVARIATES AND GENERATING PLOTS"
		#exit 1
	fi


	end=`date +%s`
	runtime=$((end-start))
	printf '\nExecuting time: '$runtime 



fi



### If just mapping step, then exit



if [ "$analysis" = "mapping" ]; then

	rm -r $TMP
	exit 0

fi






printf "\n\n\n- HAPLOTYPECALLER (GATK)"
printf "\n--------------------------\n"


#Ready to call for Variants.
mkdir $HCGVCFD
printf '\nmkdir haplotype_caller_gvcf_data' 

#HaplotypeCaller for each sample for later joint genotyping.
printf '\nGATK HaplotypeCallerGVCF for '${sample}' STARTS'
start=`date +%s`

if [ "$intervals" != "True" ]; then
	gatk HaplotypeCaller --tmp-dir=$TMP \
	-R $HG19/ucsc.hg19.fasta \
	-I $bqsr_bamfile \
	-ERC GVCF \
	-O $HCGVCFD/${sample}.g.vcf \
	-G StandardAnnotation \
	-G AS_StandardAnnotation \
	-G StandardHCAnnotation \
	-A FisherStrand -A StrandOddsRatio -A RMSMappingQuality -A MappingQualityRankSumTest -A ReadPosRankSumTest -A DepthPerSampleHC -A BaseQualityRankSumTest -A ExcessHet  \
	--annotate-with-num-discovered-alleles=true 
else
	gatk HaplotypeCaller --tmp-dir=$TMP \
	-R $HG19/ucsc.hg19.fasta \
	-I $bqsr_bamfile \
	-ERC GVCF \
	-O $HCGVCFD/${sample}.g.vcf \
	-G StandardAnnotation \
	-G AS_StandardAnnotation \
	-G StandardHCAnnotation \
	-A FisherStrand -A StrandOddsRatio -A RMSMappingQuality -A MappingQualityRankSumTest -A ReadPosRankSumTest -A DepthPerSampleHC -A BaseQualityRankSumTest -A ExcessHet \
	--annotate-with-num-discovered-alleles=true \
	-L $panel -ip 5000
fi


if [ "$?" = "0" ]; then
	printf '\nEXIT STATUS: 0'
	printf '\nGATK HaplotypeCallerGVCF ERC GVCF for '${sample}' DONE\n' 
	if [ "$cvcf" = "True" ]; then 
		printf '\n'$HCGVCFD/${sample}.g.vcf >> $HCGVCFD/my_list_of_gvcfs_files_to_combine_$run.list
	fi
	if [ "$removebam" = "single" ]; then
		basename=${bqsr_bamfile%.*}
		rm ${basename}.bam
		rm ${basename}.bai
	fi	
else
	printf "\nERROR: PROBLEMS WITH HAPLOTYPECALLER"
	exit 1
fi


end=`date +%s`
runtime=$((end-start))
printf '\nExecuting time: '$runtime 





# if joint vcf exit here


if [ "$cvcf" = "True" ]; then

 	rm -r $TMP
	exit 0

fi






printf "\n\n\n- GENOTYPECALLER (GATK) "
printf "\n--------------------------\n"


#GenotypeGVCFs into final VCF

mkdir $GVCFD
printf '\nmkdir genotyped_data_vcf'

#GenotypeGVCF into final VCF
printf  '\nUsing GATK GenotypeGVCFs for final VCF'
start=`date +%s`

gatk GenotypeGVCFs --tmp-dir=$TMP \
    -R $HG19/ucsc.hg19.fasta \
    -V $HCGVCFD/${sample}.g.vcf \
    -G StandardAnnotation \
    -O $GVCFD/genotyped_data_${sample}.vcf 

if [ "$?" = "0" ]; then
	printf '\nEXIT STATUS: 0'
	printf  '\nGATK GenotypeGVCFs for '${sample}' DONE'

else
	printf "\nERROR: PROBLEMS WITH GENOTYPECALLER"
	exit 1
fi


end=`date +%s`
runtime=$((end-start))
printf '\nExecuting time: '$runtime 










printf "\n\n\n- HARD FILTERING (classical way GATK)"
printf "\n---------------------------------------\n"


#HARD FILTERING
#First step extacting the SNP's
#Second step extracting the INDEL's


mkdir $VFVCFD
printf "mkdir snv_results\n"


### SNPs

#1.Extract the SNP's from the call set.
#printf "\nExtract the SNP's from the call set."
start=`date +%s`

gatk SelectVariants --tmp-dir=$TMP \
-R $HG19/ucsc.hg19.fasta \
-V $GVCFD/genotyped_data_${sample}.vcf \
--select-type-to-include SNP \
-O $VFVCFD/selected_raw_snp_${sample}.vcf
s1="$?"


#2.Apply the filters to the SNP's callset.

gatk VariantFiltration --tmp-dir=$TMP \
-R $HG19/ucsc.hg19.fasta \
-V $VFVCFD/selected_raw_snp_${sample}.vcf \
--filter-expression "QD < 2.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
--filter-name "my_SNP_filter" \
-O $VFVCFD/filtered_SNP_data_${sample}.vcf
s2="$?"




### INDELs

#3. Extract the INDELS from the ORIGINAL call set.
gatk SelectVariants --tmp-dir=$TMP \
-R $HG19/ucsc.hg19.fasta \
-V $GVCFD/genotyped_data_${sample}.vcf \
--select-type-to-include INDEL \
-O $VFVCFD/selected_raw_indels_${sample}.vcf
s3="$?"



#4.Apply the filters to the INDEL's callset.
gatk VariantFiltration --tmp-dir=$TMP \
-R $HG19/ucsc.hg19.fasta \
-V $VFVCFD/selected_raw_indels_${sample}.vcf \
--filter-expression "QD < 2.0  || ReadPosRankSum < -20.0" \
--filter-name "my_INDEL_filter" \
-O $VFVCFD/filtered_INDEL_data_${sample}.vcf
s4="$?"


# Combine Variants after using SNPS and INDELS filtering into a single file and get it ready for annotation.

gatk MergeVcfs --TMP_DIR=$TMP \
-R $HG19/ucsc.hg19.fasta \
-I $VFVCFD/filtered_SNP_data_${sample}.vcf \
-I $VFVCFD/filtered_INDEL_data_${sample}.vcf \
-O $VFVCFD/${sample}_raw.vcf
s5="$?"



if [ "$s1" = "0"  ] &&  [ "$s2" = "0" ] &&  [ "$s3" = "0" ] &&  [ "$s4" = "0" ] &&  [ "$s5" = "0" ]; then
	printf '\nEXIT STATUS: 0'
	printf  '\nHARD FILTERING for '${sample}' DONE'

	rm $GVCFD/genotyped_data_${sample}.vcf*
	rm $VFVCFD/selected_raw_snp_${sample}.vcf*
	rm $VFVCFD/filtered_SNP_data_${sample}.vcf*
	rm $VFVCFD/selected_raw_indels_${sample}.vcf*
	rm $VFVCFD/filtered_INDEL_data_${sample}.vcf*

else
	printf "\nERROR: PROBLEMS WITH HARD FILTERING"
	exit 1
fi


end=`date +%s`
runtime=$((end-start))
printf '\nExecuting time: '$runtime 




if [ "$cvcf" = "True" ]; then

	rm -r $TMP
	exit 0

fi



# if not combined analysis we remove all g.vcfs

rm $HCGVCFD/${sample}*




printf "\n\n\n- HOMOZYGOSITY (PLINK) "
printf "\n------------------------\n"


start=`date +%s`


params="--allow-extra-chr --homozyg --homozyg-window-het 1 --vcf-filter"
mkdir $PLINK

plink $params --vcf $VFVCFD/${sample}_raw.vcf --out $PLINK/${sample} 1>&2


if [ "$?" = "0"  ]; then
	echo -e  '\nEXIT STATUS: 0'
	echo -e   '\nPLINK HOMOZYGOSITY for '${run}' DONE'
else
	echo -e  "ERROR: PROBLEMS WITH PLINK"
	exit 1
fi 


end=`date +%s`
runtime=$((end-start))
echo -e  '\nExecuting time: '$runtime 












printf "\n\n\n- VARIANT ANNOTATION (VEP - ENSEMBL) "
printf "\n---------------------------------------\n"


VCF_IN="${VFVCFD}/${sample}_raw.vcf"
VCF_OUT="${VEPVCFAD}/${sample}_filtered.vcf"

VCF_FILTERED="${VEPVCFAD}/${sample}_annotated.vcf"
VCF_FILTERED_2="${VEPVCFAD}/${sample}_filtered_annotated_can_conseq.vcf"
VCF_FILTERED_3="${VEPVCFAD}/${sample}_filteredAnnotated.vcf" 
VCF_FINAL="${VEPVCFAD}/${sample}_filteredAnnotated.txt"
VCF_FINAL_PVM="${VEPVCFAD}/${sample}_filteredAnnotated_pvm.txt"


start=`date +%s`










printf '\n\nFiltering by quality and chromosome...'

perl $FILTER_VEP \
-i $VCF_IN -o $VCF_OUT \
--filter "(QUAL > 100) and (FILTER = PASS) and (CHROM in chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chr23,chrX,chrY)" --force_overwrite 




if [ "$?" = "0" ]; then
	printf '\nEXIT STATUS: 0'
	printf '\nVEP FILTERING quality and chromosome for '${sample}' DONE\n' 

else
	printf "\nERROR: PROBLEMS WITH VEP FILTERING FOR QUALITY AND CHR"
	exit 1
fi














printf '\n\nVEP annotation...'

perl $VEP \
--cache --offline --hgvs --refseq --dir $VEP_CACHE --dir_plugins $PLUGIN_DIR --v --fork $threads --assembly GRCh37 --fasta $VEP_FASTA --force_overwrite  --no_stats \
--biotype --regulatory --protein --symbol --allele_number --numbers --domains --uniprot --variant_class \
--canonical --vcf \
--sift p --polyphen p --af --max_af \
--format vcf \
--pubmed \
--plugin dbscSNV,$PLUGIN_DBS/dbscSNV1.1_GRCh37.txt.gz \
--plugin LoFtool,$PLUGIN_DIR/LoFtool_scores.txt \
--plugin ExACpLI,$PLUGIN_DIR/ExACpLI_values.txt \
--plugin dbNSFP,${dbNSFP_DB},gnomAD_exomes_AF,gnomAD_exomes_NFE_AF,1000Gp3_AF,1000Gp3_EUR_AF,ExAC_AF,ExAC_EAS_AF,ExAC_NFE_AF,ExAC_Adj_AF,rs_dbSNP150,phyloP20way_mammalian,phyloP20way_mammalian_rankscore,phastCons20way_mammalian,phastCons20way_mammalian_rankscore,GERP++_RS,GERP++_RS_rankscore,LRT_pred,MutationTaster_pred,MutationAssessor_pred,FATHMM_pred,PROVEAN_pred,MetaLR_pred,MetaSVM_pred,M-CAP_pred,Interpro_domain,GTEx_V6p_gene,GTEx_V6p_tissue \
--plugin MaxEntScan,$PLUGIN_DBS/maxEntScan \
--custom ${PLUGIN_DBS}/gnomad.genomes.r2.0.1.sites.noVEP.vcf.gz,gnomADg,vcf,exact,0,AF_NFE,AF_Male,AF_Female,Hom,POPMAX,AF_POPMAX \
--custom ${PLUGIN_DBS}/Kaviar-160204-Public/vcfs/Kaviar-160204-Public-hg19.vcf.gz,kaviar,vcf,exact,0,AF,AC,AN \
--custom ${CCS_DB},gnomAD_exomes_CCR,bed,overlap,0,ccr_pct \
--plugin CADD,${PLUGIN_DBS}/InDels.tsv.gz,${PLUGIN_DBS}/whole_genome_SNVs.tsv.gz \
-i $VCF_OUT -o $VCF_FILTERED



if [ "$?" = "0" ]; then
	printf '\nEXIT STATUS: 0'
	printf '\nVEP ANNOTATION for '${sample}' DONE\n' 
	rm $VCF_OUT

else
	printf "\nERROR: PROBLEMS WITH VEP ANNOTATION"
	exit 1
fi







printf '\n\nFiltering by population frequencies...\n'


perl $FILTER_VEP \
-i $VCF_FILTERED -o $VCF_FILTERED_3 \
--filter "(DP > 10) and (MAX_AF < 0.05 or not MAX_AF)" \
--force_overwrite


if [ "$?" = "0" ]; then
	printf '\nEXIT STATUS: 0'
	printf '\n VEP FREQUENCY FILTERING for '${sample}' DONE\n' 
	rm $VCF_FILTERED

else
	printf "\nERROR: PROBLEMS WITH VEP FREQUENCY FILTERING"
	exit 1
fi







printf '\nVEP annotation and filtering DONE.\n'


end=`date +%s`
runtime=$((end-start))
printf '\nExecuting time: '$runtime 






printf "\n\n\n- OUTPUT PROCESSING "
printf "\n------------------------\n"

start=`date +%s`

printf '\nFrom VEP annotated VCF to txt file...\n'
start=`date +%s`


python $utilitiesPath/vep2tsv_woFreq.py $VCF_FILTERED_3 \
-o $VCF_FINAL -t -g -f GT,AD,DP

if [ "$?" = "0" ]; then
	printf '\nEXIT STATUS: 0'
	printf '\nVCF PROCESSING for '${sample}' DONE'

else
	printf "\nERROR: PROBLEMS WITH VCF PREPROCESSING"
	exit 1
fi


end=`date +%s`
runtime=$((end-start))
printf '\nExecuting time: '$runtime 








printf "\n\n\n- ADDING LOH INFO"
printf "\n-----------------------\n"

printf '\nAnnotating extra features...\n'

start=`date +%s`

#python $utilitiesPath/LOHmerge.py $PLINK/${sample}.hom $VCF_FINAL ${VCF_FINAL_PVM}

echo python $utilitiesPath/PVM_Cluster.py $VCF_FINAL -l $local -k $PLINK/${sample}.hom -o ${VCF_FINAL_PVM} -P $pathology

python $utilitiesPath/PVM_Cluster.py $VCF_FINAL -l $local -k $PLINK/${sample}.hom -o ${VCF_FINAL_PVM} -P $pathology


if [ "$?" = "0" ]; then
	printf '\nEXIT STATUS: 0'
	printf '\nPOST-VEP MODIFICATIONS for '${sample}' DONE'
	rm $VCF_FINAL

else
	printf "\nERROR: PROBLEMS WITH POST-VEP MODIFICATIONS"
	exit 1

fi

end=`date +%s`
runtime=$((end-start))
printf '\nExecuting time: '$runtime 




# Filter PVM files based on gene list

if [ "$genefilter" != "False" ]; then


	start=`date +%s`

	VCF_FINAL_PVM_FILTER="${VEPVCFAD}/${sample}_filteredAnnotated_pvm_GENELIST.txt"
	
	python $utilitiesPath/filtering_geneList.py -i ${VCF_FINAL_PVM} -f ${genefilter} -o ${VCF_FINAL_PVM_FILTER}


	if [ "$?" = "0" ]; then
		printf '\nEXIT STATUS: 0'
		printf '\nVARIANT FILTERING BY GENE LIST for '${sample}' DONE'

	else
		printf "\nERROR: PROBLEMS WITH GENE FILTERING"
		exit 1

	fi

	end=`date +%s`
	runtime=$((end-start))
	printf '\nExecuting time: '$runtime 

fi



# Remove temporary files.

printf '\n'
rm -r $TMP







