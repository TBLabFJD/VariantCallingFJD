#!/bin/sh

#################################
### FJD pipeline -  Filtering ###
#################################




# ARGUMENTS:

local=$1
ftype=$2
MDAP=$3
#Haplotype Results
HCGVCFD="${MDAP}/haplotype_caller_gvcf_data"
#Genotyped Results
GVCFD="${MDAP}/genotyped_vcf_data"
# Hard filtering, vep filtering and vep annotation results
VEPVCFAD="${MDAP}/snv_results"
#Temporary directory
name=$4
TMP=$MDAP/tmp_${name}
HG19=$5




# REQUIRED MODULES AND DATABASES:

if [ "$local" != "True" ]; then


	CNN_REF=/home/proyectos/bioinfo/references/CNNgatk


else

	CNN_REF=/home/proyectos/bioinfo/references/CNNgatk


fi




# COMMANDS

# IMPORTANT: CNN method can be only applied to single samples since the default models were trained on single-sample VCFs. 


if [ "$ftype" != "HF" ]; then



	printf "\n\n\n- Convolutional Neural Network (CNN)"
	printf "\n---------------------------------------\n"



	# 1. Annotate a VCF with scores
	# 2. Filtering variants according to CNN score.


	gatk CNNScoreVariants \
	-I $HCGVCFD/${name}_bamout.bam \
	-V $GVCFD/genotyped_data_${name}.vcf  \
	-O $GVCFD/genotyped_data_cnn_scored_${name}.vcf  \
	-R $HG19/ucsc.hg19.fasta \
	--tensor-type read_tensor


	gatk FilterVariantTranches \
	-V $GVCFD/genotyped_data_cnn_scored_${name}.vcf \
	--resource $CNN_REF/b37_1000G_omni2.5.b37.vcf.gz \
	--resource $CNN_REF/b37_hapmap_3.3.b37.vcf.gz \
	--info-key CNN_2D \
	--snp-tranche 95.9 \
	--indel-tranche 95.0 \
	-O $GVCFD/genotyped_data_cnn_scored_filtered_${name}.vcf \
	--invalidate-previous-filters



	# # compare 2 vcf files

	# gatk Concordance \
	# -truth gs://gatk-tutorials/$WORKSHOP/2-germline/CNNScoreVariants/vcfs/hg001_na12878_b37_truth.vcf.gz \
	# -eval /home/jupyter-user/CNN/Output/my_2d_filtered.vcf \
	# -L 20:1000000-1432828 \
	# -S /home/jupyter-user/CNN/Output/2d_filtered_concordance.txt




else




	printf "\n\n\n- Hard Filtering (classical way GATK)"
	printf "\n---------------------------------------\n"


	#HARD FILTERING
	#First step extacting the SNP's
	#Second step extracting the INDEL's


	mkdir $VEPVCFAD
	printf "mkdir snv_results\n"


	### SNPs

	#1.Extract the SNP's from the call set.
	#printf "\nExtract the SNP's from the call set."
	start=`date +%s`

	gatk SelectVariants --tmp-dir=$TMP \
	-R $HG19/ucsc.hg19.fasta \
	-V $GVCFD/genotyped_data_${name}.vcf \
	--select-type-to-include SNP \
	-O $VEPVCFAD/selected_raw_snp_${name}.vcf
	s1="$?"


	#2.Apply the filters to the SNP's callset.

	gatk VariantFiltration --tmp-dir=$TMP \
	-R $HG19/ucsc.hg19.fasta \
	-V $VEPVCFAD/selected_raw_snp_${name}.vcf \
	--filter-expression "QD < 2.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
	--filter-name "my_SNP_filter" \
	-O $VEPVCFAD/filtered_SNP_data_${name}.vcf
	s2="$?"




	### INDELs

	#3. Extract the INDELS from the ORIGINAL call set.
	gatk SelectVariants --tmp-dir=$TMP \
	-R $HG19/ucsc.hg19.fasta \
	-V $GVCFD/genotyped_data_${name}.vcf \
	--select-type-to-include INDEL \
	-O $VEPVCFAD/selected_raw_indels_${name}.vcf
	s3="$?"



	#4.Apply the filters to the INDEL's callset.
	gatk VariantFiltration --tmp-dir=$TMP \
	-R $HG19/ucsc.hg19.fasta \
	-V $VEPVCFAD/selected_raw_indels_${name}.vcf \
	--filter-expression "QD < 2.0  || ReadPosRankSum < -20.0" \
	--filter-name "my_INDEL_filter" \
	-O $VEPVCFAD/filtered_INDEL_data_${name}.vcf
	s4="$?"


	# Combine Variants after using SNPS and INDELS filtering into a single file and get it ready for annotation.

	gatk MergeVcfs --TMP_DIR=$TMP \
	-R $HG19/ucsc.hg19.fasta \
	-I $VEPVCFAD/filtered_SNP_data_${name}.vcf \
	-I $VEPVCFAD/filtered_INDEL_data_${name}.vcf \
	-O $VEPVCFAD/${name}_raw.vcf
	s5="$?"



	if [ "$s1" = "0"  ] &&  [ "$s2" = "0" ] &&  [ "$s3" = "0" ] &&  [ "$s4" = "0" ] &&  [ "$s5" = "0" ]; then
		printf '\nEXIT STATUS: 0'
		printf  '\nHARD FILTERING for '${name}' DONE'

		rm $GVCFD/genotyped_data_${name}.vcf*
		rm $VEPVCFAD/selected_raw_snp_${name}.vcf*
		rm $VEPVCFAD/filtered_SNP_data_${name}.vcf*
		rm $VEPVCFAD/selected_raw_indels_${name}.vcf*
		rm $VEPVCFAD/filtered_INDEL_data_${name}.vcf*

	else
		printf "\nERROR: PROBLEMS WITH HARD FILTERING"
		exit 1
	fi


	end=`date +%s`
	runtime=$((end-start))
	printf '\nExecuting time: '$runtime 






fi


