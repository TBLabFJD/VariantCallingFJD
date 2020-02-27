#!/bin/sh

############################
### TASK: SNV filtering ###
###########################

# Two options: CNN or Hard Filtering


###############
## ARGUMENTS ##
###############

local=$1
run=$2
MDAP=$3
name=$4 # sample name
REF=$5
ftype=$6 # when multi VCF is used only "HF" can be used
cvcf=$7




#####################
## MODULES AND DBs ##
#####################


if [ "$local" != "True" ]; then

	module load miniconda/2.7
	module load samtools/1.9
	module load picard/2.18.9
	module load gatk/4.1.2.0
	module load bedtools/2.27.0
	module load R
	module load bcftools/1.3
	module load gcc/7.3.0
	source /usr/local/miniconda/python-2.7/etc/profile.d/conda.sh
	conda activate gatk

	alias gatk='java -jar /usr/local/bioinfo/gatk/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar'	
	CNN_REF=/home/proyectos/bioinfo/references/CNNgatk


	softwareFile="${MDAP}/software_${run}.txt"
	title="SNV FILTERING"
	if [ ! -f $softwareFile ] || [ `grep -q $title $softwareFile` ] ; then 

		printf "SNV FILTERING:\n" >> ${softwareFile}
		module list 2>> ${softwareFile}
	
	fi




else

	export SFT=/mnt/genetica3/marius/pipeline_practicas_marius/software
	alias gatk='java  -Xmx10g -jar $SFT/gatk/build/libs/gatk-package-4.0.6.0-22-g9d9484f-SNAPSHOT-local.jar'
	CNN_REF=/home/proyectos/bioinfo/references/CNNgatk

	softwareFile="${MDAP}/software_${run}.txt"
	echo $softwareFile
	title="SNV FILTERING"
	
	if [ ! -f $softwareFile ] || [ `grep -q $title $softwareFile` ] ; then 

		printf "SNV FILTERING:\n" >> ${softwareFile}
		
		printf "\nGATK VERSION\n" >> ${softwareFile}
		gatk ApplyBQSR 2>&1 | head -n4 | tail -n1 >> ${softwareFile}


	fi




fi






###############
## VARIABLES ##
###############

# Ref Fasta
fasta=$REF


# genotyping folder
GEN=$MDAP/genotyping 

# snv results folder
SNV="${MDAP}/snv_results"
mkdir $SNV


# Temporal Folder
TMP=$MDAP/${sample}_tmp
mkdir $TMP







##############
## COMMANDS ##
##############


# IMPORTANT: CNN method can be only applied to single samples since the default models were trained on single-sample VCFs. 


if [ "$ftype" != "HF" ]; then



	printf "\n\n\n- Convolutional Neural Network (CNN)"
	printf "\n---------------------------------------\n"



	# 1. Annotate a VCF with scores
	# 2. Filtering variants according to CNN score.


	gatk CNNScoreVariants \
	--disable-avx-check \
	-I $GEN/${name}_bamout.bam \
	-V $GEN/${name}.vcf  \
	-O $SNV/${name}_cnnScored.vcf  \
	-R $fasta \
	--tensor-type read_tensor

	echo $GEN/${name}_cnnScored.vcf
	echo ${CNN_REF}/b37_hapmap_3.3.b37.chr.vcf.gz

	gatk FilterVariantTranches \
	-V  $GEN/${name}_cnnScored.vcf \
	--resource ${CNN_REF}/b37_hapmap_3.3.b37.chr.vcf.gz \
	--resource ${CNN_REF}/b37_1000G_omni2.5.b37.chr.vcf.gz \
	--info-key CNN_2D \
	--snp-tranche 99.95 \
	--indel-tranche 99.4 \
	-O $SNV/${name}_gatkLabeled.vcf \
	--invalidate-previous-filters


	#bcftools view -i 'FILTER="."' $SNV/${name}_cnnLabeled.vcf  > $SNV/${name}_gatkFiltered.vcf
	#s6="$?"


	#--resource ${CNN_REF}/b37_1000G_omni2.5.b37.vcf.gz \


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


	### SNPs

	#1.Extract the SNP's from the call set.
	#printf "\nExtract the SNP's from the call set."
	start=`date +%s`

	gatk SelectVariants --tmp-dir=$TMP \
	-R $fasta \
	-V $GEN/${name}.vcf \
	--select-type-to-include SNP \
	-O $SNV/${name}_selected_raw_snp.vcf
	s1="$?"


	#2.Apply the filters to the SNP's callset.

	gatk VariantFiltration --tmp-dir=$TMP \
	-R $fasta \
	-V $SNV/${name}_selected_raw_snp.vcf \
	--filter-expression "QD < 2.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
	--filter-name "my_SNP_filter" \
	-O $SNV/${name}_filtered_snp.vcf
	s2="$?"




	### INDELs

	#3. Extract the INDELS from the ORIGINAL call set.
	gatk SelectVariants --tmp-dir=$TMP \
	-R $fasta \
	-V $GEN/${name}.vcf \
	--select-type-to-include INDEL \
	-O $SNV/${name}_selected_raw_indel.vcf
	s3="$?"



	#4.Apply the filters to the INDEL's callset.
	gatk VariantFiltration --tmp-dir=$TMP \
	-R $fasta \
	-V $SNV/${name}_selected_raw_indel.vcf \
	--filter-expression "QD < 2.0  || ReadPosRankSum < -20.0" \
	--filter-name "my_INDEL_filter" \
	-O $SNV/${name}_filtered_indel.vcf
	s4="$?"


	# Combine Variants after using SNPS and INDELS filtering into a single file and get it ready for annotation.

	gatk MergeVcfs --TMP_DIR=$TMP \
	-R $HG19/ucsc.hg19.fasta \
	-I $SNV/${name}_filtered_snp.vcf \
	-I $SNV/${name}_filtered_indel.vcf \
	-O $SNV/${name}_gatkLabeled.vcf
	s5="$?"

	#bcftools view -i 'FILTER="PASS"' $SNV/${name}_hLabeled.vcf  > $SNV/${name}_gatkFiltered.vcf
	#s6="$?"


	if [ "$s1" = "0"  ] &&  [ "$s2" = "0" ] &&  [ "$s3" = "0" ] &&  [ "$s4" = "0" ] &&  [ "$s5" = "0" ]; then
		printf '\nEXIT STATUS: 0'
		printf  '\nHARD FILTERING for '${name}' DONE'

		rm $GEN/${name}.vcf*
		rm $SNV/${name}_selected_raw_indel.vcf*
		rm $SNV/${name}_selected_raw_snp.vcf*
		rm $SNV/${name}_filtered_indel.vcf*
		rm $SNV/${name}_filtered_snp.vcf*
	else
		printf "\nERROR: PROBLEMS WITH HARD FILTERING"
		exit 1
	fi


	end=`date +%s`
	runtime=$((end-start))
	printf '\nExecuting time: '$runtime 






fi






if [ "$cvcf" = "True" ]; then 
	
	printf '\n' $SNV/${name}_raw.vcf >> $SNV/filteredVCFs_to_combine_$run.list

fi

