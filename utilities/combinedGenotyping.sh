#!/bin/sh

############################################################
### FJD pipeline - Combined Genotyping and Trio analysis ###
############################################################


### FJD PIPELINE ARGUMENTS:

MDAP=$1
run=$2
threads=$3
HG19=$4
ped=$5
local=$6
pathology=$7
single=$8
genefilter=$9
utilitiesPath=${10}





# REQUIRED MODULES AND DATABASES:

if [ "$local" != "True" ]; then

	module load python/gnu/4.4.7/2.7.3
	module load samtools/1.9
	module load picard/2.18.9
	module load gatk/4.1.2.0
	module load vep/release95
	module load bedtools/2.27.0
	module load R
	#module load plink
	alias picard='java -jar /usr/local/bio/picard/2.18.9/picard.jar'
	alias gatk='java -jar /usr/local/bio/GATK/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar'	#alias gatk3='java -jar /usr/local/bio/GATK/gatk-4.0.5.1/gatk-package-4.0.5.1-local.jar'
	alias bcftoolsl="/home/proyectos/bioinfo/software/bcftools-1.9/bcftools"
	VEP="/usr/local/bio/vep/vep"
 	FILTER_VEP='/usr/local/bio/vep/filter_vep'
	VEP_CACHE='/usr/local/bio/vep/t/testdata/cache/homo_sapiens'
	VEP_FASTA="/home/proyectos/bioinfo/references/VEPfasta/Homo_sapiens.GRCh37.dna.primary_assembly.fa"
	PLUGIN_DIR=/usr/local/bio/vep/plugins-95/VEP_plugins-release-95
	PLUGIN_DBS="/home/proyectos/bioinfo/references/VEPdbs"
	dbNSFP_DB="${PLUGIN_DBS}/dbNSFP3.5a_hg19.gz"
	softwareFile="${MDAP}/software_${run}.txt"

	title="GENOTYPING"
	if grep -q $title ${softwareFile}; then
		echo found
	else
		echo -e "COMBINED GENOTYPING:\n" >> ${softwareFile}
		module list 2>> ${softwareFile}
		module list 2>> ${softwareFile}

	fi



else
	export SFT=/mnt/genetica3/marius/pipeline_practicas_marius/software
	alias samtools='$SFT/samtools/samtools'
	alias picard='java -jar $SFT/picard/build/libs/picard.jar'
	alias gatk='java -jar $SFT/gatk/build/libs/gatk-package-4.0.6.0-22-g9d9484f-SNAPSHOT-local.jar'
	#alias gatk3='java -jar /mnt/genetica3/GeneticaPipeDB_updated/gatk-3.8/GenomeAnalysisTK.jar'
	alias plink='$SFT/plink/plink'
	VEP="${SFT}/variant_effect_predictor/ensembl-vep/vep"
	FILTER_VEP="${SFT}/variant_effect_predictor/ensembl-vep/filter_vep"
	VEP_CACHE='/mnt/genetica3/marius/pipeline_practicas_marius/software/variant_effect_predictor/.vep'
	VEP_FASTA="${VEP_CACHE}/homo_sapiens/93_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa"
	PLUGIN_DIR="${VEP_CACHE}/Plugins"
	PLUGIN_DBS="${VEP_CACHE}/dbs"
	dbNSFP_DB="${PLUGIN_DBS}/dbNSFP_hg19.gz"


	softwareFile="${MDAP}/software_${run}.txt"
	title="GENOTYPING"

	if grep -q $title ${softwareFile}; then
		echo found
	else
		echo -e "COMBINED GENOTYPING:\n" >> ${softwareFile}
		
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





### VARIABLES:

#Haplotype_caller_gvcf_data:calling for SNPs and indels via local re-assembly of HAPLOTYPES using HAPLOTYPECALLER by GATK.
HCGVCFD="${MDAP}/haplotype_caller_gvcf_data"

#If analyzing more than one sample or an family trio, use COMBINE_GVCFS to  combine them into a single GVCF.
CGVCFD="${MDAP}/combined_gvcf_data"

#Genotyped_vcf_data: single or single-multisample GVCF as input, output will be a VCF.
GVCFD="${MDAP}/genotyped_vcf_data"
	
#Variant_filtration_vcf_data:hard filtering process to select based on the INFO and FORMAT annotations(QD,MQO,FS,MQ,ReadPosrankSum)
VFVCFD="${MDAP}/snv_results"

#Vep_vcf_annotated_data: annotations added to the CSQ tag in INFO columnd to the VCF format.
VEPVCFAD="${MDAP}/snv_results"

# LOH regions
PLINK="${MDAP}/plink_results"

TMP=$MDAP/tmp_$run
mkdir $TMP








### RUN:



echo -e "\n\n\n- COMBINED GVCFS (GATK) "
echo -e "-----------------------------\n"

echo -e 'mkdir combined_gvcf_data'
mkdir $CGVCFD
echo -e '\nUsing GATK COMBINEGVCFs for merging GVCFs'
start=`date +%s`

gatk CombineGVCFs --tmp-dir=$TMP \
-R $HG19/ucsc.hg19.fasta \
--variant $HCGVCFD/my_list_of_gvcfs_files_to_combine_$run.list \
-O $CGVCFD/combinedg_$run.vcf


if [ "$?" = "0" ]; then
	echo -e  '\nEXIT STATUS: 0'
	echo -e '\nGATK COMBINEGVCFs DONE'

else
	echo -e  "ERROR: PROBLEMS WITH COMBINEGVCFs"
	exit 1
fi

end=`date +%s`
runtime=$((end-start))
echo -e  '\nExecuting time: '$runtime 







echo -e "\n\n\n- JOINT GENOTYPING (GATK) "
echo -e "-----------------------------\n"


mkdir $GVCFD
echo -e 'mkdir genotyped_data_vcf'
echo -e '\nUsing GATK GenotypeGVCFs for final VCF'

start=`date +%s`

if [ "$ped" != "null" ]; then

	gatk GenotypeGVCFs --tmp-dir=$TMP  \
	-R $HG19/ucsc.hg19.fasta \
	-V $CGVCFD/combinedg_$run.vcf \
	-O $GVCFD/genotyped_data_$run.vcf \
	-ped $ped

else

	gatk GenotypeGVCFs --tmp-dir=$TMP  \
	-R $HG19/ucsc.hg19.fasta \
	-V $CGVCFD/combinedg_$run.vcf \
	-O $GVCFD/genotyped_data_$run.vcf

fi


if [ "$?" = "0" ]; then
	echo -e  '\nEXIT STATUS: 0'
	echo -e   '\nGATK GenotypeGVCFs for '${run}' DONE'
	rm $CGVCFD/combinedg_$run.vcf*

else
	echo -e  "ERROR: PROBLEMS WITH GENOTYPECALLER"
	exit 1
fi


end=`date +%s`
runtime=$((end-start))
echo -e  '\nExecuting time: '$runtime 















echo -e "\n\n\n- HARD FILTERING (GATK) "
echo -e "-----------------------------\n"

#HARD FILTERING
#First step extacting the SNP's
#Second step extracting the INDEL's


mkdir $VFVCFD
echo -e  "mkdir snv_results"

start=`date +%s`


### SNPs
echo -e "Extract and filter SNP's from the call set."

	#1.Extract the SNP's from the call set.
gatk SelectVariants --tmp-dir=$TMP \
-R $HG19/ucsc.hg19.fasta \
-V $GVCFD/genotyped_data_$run.vcf \
--select-type-to-include SNP \
-O $VFVCFD/selected_raw_snp_$run.vcf
s1="$?"

	#2.Apply the filters to the SNP's callset.
gatk VariantFiltration --tmp-dir=$TMP \
-R $HG19/ucsc.hg19.fasta \
-V $VFVCFD/selected_raw_snp_$run.vcf \
--filter-expression "QD < 2.0  || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
--filter-name "my_SNP_filter" \
-O $VFVCFD/filtered_SNP_data_$run.vcf
s2="$?"

### INDELs
echo -e "Extract and filter INDEL's from the call set."

	#3. Extract the INDELS from the ORIGINAL call set.
gatk SelectVariants --tmp-dir=$TMP \
-R $HG19/ucsc.hg19.fasta \
-V $GVCFD/genotyped_data_$run.vcf \
--select-type-to-include INDEL \
-O $VFVCFD/selected_raw_indels_$run.vcf
s3="$?"

	#4.Apply the filters to the INDEL's callset.
gatk VariantFiltration --tmp-dir=$TMP \
-R $HG19/ucsc.hg19.fasta \
-V $VFVCFD/selected_raw_indels_$run.vcf \
--filter-expression "QD < 2.0  || ReadPosRankSum < -20.0" \
--filter-name "my_INDEL_filter" \
-O $VFVCFD/filtered_INDEL_data_$run.vcf
s4="$?"

# 5. Combine Variants after using SNPS and INDELS filtering into a single file and get it ready for annotation.
gatk MergeVcfs --TMP_DIR=$TMP \
-R $HG19/ucsc.hg19.fasta \
-I $VFVCFD/filtered_SNP_data_$run.vcf \
-I $VFVCFD/filtered_INDEL_data_$run.vcf \
-O $VFVCFD/raw_$run.vcf
s5="$?"


if [ "$s1" = "0"  ] &&  [ "$s2" = "0" ] &&  [ "$s3" = "0" ] &&  [ "$s4" = "0" ] &&  [ "$s5" = "0" ]; then
	echo -e  '\nEXIT STATUS: 0'
	echo -e   '\nHARD FILTERING for '${run}' DONE'


	for f in $(cat $HCGVCFD/my_list_of_gvcfs_files_to_combine_$run.list) ; do 
		y=${f%%.*}
		rm $y*
	done

	rm $HCGVCFD/my_list_of_gvcfs_files_to_combine_$run.list
	rm $GVCFD/genotyped_data_$run.vcf*
	rm $VFVCFD/selected_raw_snp_${run}.vcf*
	rm $VFVCFD/filtered_SNP_data_${run}.vcf*
	rm $VFVCFD/selected_raw_indels_${run}.vcf*
	rm $VFVCFD/filtered_INDEL_data_${run}.vcf*


else
	echo -e  "ERROR: PROBLEMS WITH HARD FILTERING"
	exit 1
fi


end=`date +%s`
runtime=$((end-start))
echo -e  '\nExecuting time: '$runtime 










echo -e  "\n\n\n- Split VCF into single-sample VCFs "
echo -e  "---------------------------------------------\n"


for samplee in `bcftoolsl query -l $VFVCFD/raw_${run}.vcf`
do
	gatk SelectVariants --exclude-non-variants -R $HG19/ucsc.hg19.fasta \
	-V $VFVCFD/raw_${run}.vcf \
	-O $VFVCFD/${samplee}_raw.vcf -sn $samplee --exclude-non-variants

done 


if [ "$?" = "0"  ]; then
	echo -e  '\nEXIT STATUS: 0'
	echo -e  '\nVCF splitting for '${run}' DONE'
else
	echo -e  "ERROR: PROBLEMS WITH VCF SPLITTING"
	exit 1
fi 









echo -e  "\n\n\n- HOMOZYGOSITY (PLINK) "
echo -e  "-------------------------\n"

start=`date +%s`

if [ "$local" != "True" ]; then
	module load plink/1.90beta
	alias plink="/usr/local/bio/Plink/plink"
fi 


params="--allow-extra-chr --homozyg --homozyg-window-het 1 --vcf-filter"
mkdir $PLINK


for samplee in `bcftoolsl query -l $VFVCFD/raw_${run}.vcf`
do
	plink $params --vcf $VFVCFD/${samplee}_raw.vcf --out $PLINK/${samplee} 1>&2

	if [ "$?" != "0"  ]; then
		echo -e  "ERROR: PROBLEMS WITH PLINK"
		exit 1
	fi 

	sed 1d $PLINK/${samplee}.hom >> $PLINK/${run}_plink.hom
	rm $PLINK/${samplee}.*

	if [ "$single" != "True" ]; then
		rm $VFVCFD/${samplee}_raw.vcf
	fi

done 


if [ "$?" = "0"  ]; then
	echo -e  '\nEXIT STATUS: 0'
	echo -e   '\nPLINK HOMOZYGOSITY for '${run}' DONE'
else
	echo -e  "ERROR: PROBLEMS WITH PLINK FILE PROCESSING"
	exit 1
fi 

end=`date +%s`
runtime=$((end-start))
echo -e  '\nExecuting time: '$runtime 



if [ "$local" != "True" ]; then
	module unload plink/1.90beta
fi











echo -e "\n\n\n- GENOTYPE POSTERIOR RECALCULATION BASED ON FAMILY PRIOR INFORMATION (GATK) "
echo -e "-------------------------------------------------------------------------------\n"


if [ "$ped" != "null" ]; then

	start=`date +%s`

	# Recalculate genotype posterior probabilities based on family prior information.
	# Compute the most likely genotype combination of trios and parent/child pairs given their genotype likelihoods and a mutation prior;

	gatk CalculateGenotypePosteriors \
	   -R $HG19/ucsc.hg19.fasta \
	   -V $VFVCFD/raw_$run.vcf \
	   -ped $ped \
	   -O $VFVCFD/filtered_INDEL_SNP_data_GP_$run.vcf \
	   --skip-population-priors  
	s1="$?"


	# Filter variants based on GQ: genotypes with GQ < 20 based on the posteriors are flagged for posterior filtering. 
	gatk VariantFiltration \
	-R $HG19/ucsc.hg19.fasta \
	-V $VFVCFD/filtered_INDEL_SNP_data_GP_$run.vcf \
	--filter-expression "GQ < 20" \
	--filter-name "lowGQ" \
	-O $VFVCFD/filtered_INDEL_SNP_data_GP_GQfiltered_$run.vcf 
	s2="$?"


	# Annotation of High and Low Confidence De Novo mutations: 

	gatk VariantAnnotator \
		-R $HG19/ucsc.hg19.fasta \
		-V $VFVCFD/filtered_INDEL_SNP_data_GP_GQfiltered_$run.vcf  \
		-O $VFVCFD/filtered_INDEL_SNP_data_GP_GQfiltered_DeNovoAnnot_$run.vcf \
		-A PossibleDeNovo  \
		-A StrandBiasBySample  \
		-A AS_FisherStrand \
		-A DepthPerAlleleBySample \
		-A DepthPerSampleHC \
		-ped $ped 
	s3="$?"


	VCF_IN="${VFVCFD}/filtered_INDEL_SNP_data_GP_GQfiltered_DeNovoAnnot_$run.vcf"


	if [ "$s1" = "0"  ] &&  [ "$s2" = "0" ] &&  [ "$s3" = "0" ] ; then
		echo -e  '\nEXIT STATUS: 0'
		echo -e '\nPosterior genotyping and de novo mutations annotation DONE\n'
		rm  $VFVCFD/filtered_INDEL_SNP_data_GP_$run.vcf*
		rm  $VFVCFD/filtered_INDEL_SNP_data_GP_GQfiltered_$run.vcf* 

	else
		echo -e  "ERROR: PROBLEMS WITH posterior genotyping and de novo mutations annotation"
		exit 1
	fi


	end=`date +%s`
	runtime=$((end-start))
	echo -e  '\nExecuting time: '$runtime 



else

	echo -e "\n NOT INPUT PED FILE: genotype posterior calculation skipped"

	VCF_IN="${VFVCFD}/raw_${run}.vcf"


fi








echo -e "\n\n\n- VARIANT ANNOTATION (VEP - ENSEMBL) "
echo -e "----------------------------------------\n"


VCF_OUT="${VEPVCFAD}/filtered_${run}.vcf"
VCF_FILTERED="${VEPVCFAD}/annotated_${run}.vcf"
VCF_FILTERED_2="${VEPVCFAD}/annotated_can_conseq_${run}.vcf"
VCF_FILTERED_3="${VEPVCFAD}/filteredAnnotated_${run}.vcf" 
VCF_FINAL="${VEPVCFAD}/filteredAnnotated_${run}.txt"
VCF_FINAL_PVM="${VEPVCFAD}/${run}_filteredAnnotated_pvm.txt"


start=`date +%s`









echo -e '\nFiltering by quality and chromosome...'

perl $FILTER_VEP \
-i $VCF_IN -o $VCF_OUT \
--filter "QUAL > 100 and FILTER = PASS" \
--filter "CHROM in chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chr23,chrX,chrY" --force_overwrite 

if [ "$?" = "0" ]; then
	printf 'EXIT STATUS: 0'
	printf '\nVEP FILTERING quality and chromosome for '${run}' DONE' 

else
	printf "ERROR: PROBLEMS WITH VEP FILTERING FOR QUALITY AND CHR"
	exit 1
fi











echo -e '\n\nVEP annotation...'

perl $VEP \
--cache --offline --hgvs --refseq --dir $VEP_CACHE --dir_plugins $PLUGIN_DIR --v --fork 16 --assembly GRCh37 --fasta $VEP_FASTA --force_overwrite \
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
--plugin CADD,${PLUGIN_DBS}/InDels.tsv.gz,${PLUGIN_DBS}/whole_genome_SNVs.tsv.gz \
-i $VCF_OUT -o $VCF_FILTERED




if [ "$?" = "0" ]; then
	printf '\nEXIT STATUS: 0'
	printf '\nVEP ANNOTATION for '${run}' DONE\n' 

else
	printf "ERROR: PROBLEMS WITH VEP ANNOTATION"
	exit 1
fi



echo -e '\n\nFiltering by population frequencies...'

perl $FILTER_VEP \
-i $VCF_FILTERED -o $VCF_FILTERED_2 \
--filter "DP > 10" \
--filter "CANONICAL = YES" \
--filter "Consequence is 3_prime_UTR_variant or Consequence is 5_prime_UTR_variant or Consequence is intron_variant or Consequence is splice_donor_variant or Consequence is splice_acceptor_variant or Consequence is splice_region_variant or Consequence is synonymous_variant or Consequence is missense_variant or Consequence is inframe_deletion or Consequence is inframe_insertion or Consequence is stop_gained or Consequence is frameshift_variant or Consequence is coding_sequence_variant" --force_overwrite 



if [ "$?" = "0" ]; then
	printf '\nEXIT STATUS: 0'
	printf '\nVEP FREQUENCY FILTERING 1 for '${run}' DONE\n' 

else
	printf "ERROR: PROBLEMS WITH VEP FREQUENCY FILTERING 1"
	exit 1
fi






perl $FILTER_VEP \
-i $VCF_FILTERED_2 -o $VCF_FILTERED_3 \
--filter "AF < 0.01 or not AF" \
--filter "(ExAC_EAS_AF < 0.01 or not ExAC_EAS_AF) and (1000Gp3_AF < 0.01 or not 1000Gp3_AF) and (1000Gp3_EUR_AF < 0.01 or not 1000Gp3_EUR_AF ) and (gnomAD_exomes_AF < 0.01 or not gnomAD_exomes_AF) and (ExAC_AF < 0.01 or not ExAC_AF) and (ExAC_Adj_AF < 0.01 or not ExAC_Adj_AF) and (ExAC_NFE_AF < 0.01 or not ExAC_NFE_AF) and (gnomAD_exomes_NFE_AF < 0.01 or not gnomAD_exomes_NFE_AF) and (gnomAD_genomes_AF < 0.01 or not gnomAD_genomes_AF)" --force_overwrite 



if [ "$?" = "0" ]; then
	printf '\nEXIT STATUS: 0'
	printf '\nVEP FREQUENCY FILTERING 2 for '${run}' DONE\n' 
	rm $VCF_OUT
	rm $VCF_FILTERED
	rm $VCF_FILTERED_2

else
	printf "ERROR: PROBLEMS WITH VEP FREQUENCY FILTERING 2"
	exit 1
fi


printf '\n\nVEP annotation and filtering DONE\n'


end=`date +%s`
runtime=$((end-start))
echo -e  '\nExecuting time: '$runtime 










echo -e  "\n\n\n- OUTPUT PROCESSING "
echo -e  "-------------------------\n"

start=`date +%s`

echo -e  '\nFrom VEP annotated VCF to txt file...\n'
start=`date +%s`


if [ "$ped" != "null" ]; then

	python $utilitiesPath/vep2tsv_woFreq.py $VCF_FILTERED_3 \
-o $VCF_FINAL -t -g -f GT,AD,DP -j -y

else

	python $utilitiesPath/vep2tsv_woFreq.py $VCF_FILTERED_3 \
-o $VCF_FINAL -t -g -f GT,AD,DP -j 

fi




if [ "$?" = "0" ]; then
	echo -e  '\nEXIT STATUS: 0'
	echo -e  '\nVCF PROCESSING for '${run}' DONE'

else
	echo -e  "ERROR: PROBLEMS WITH VCF PREPROCESSING"
	exit 1
fi


end=`date +%s`
runtime=$((end-start))
echo -e  '\nExecuting time: '$runtime 







echo -e  "\n\n\n- ADDING LOH INFO"
echo -e  "------------------------\n"

start=`date +%s`

#python $utilitiesPath/LOHmerge.py $PLINK/${run}_plink.hom $VCF_FINAL ${VCF_FINAL}_LOH.txt
python $utilitiesPath/PVM_Cluster.py $VCF_FINAL -l $local -k $PLINK/${run}_plink.hom -o ${VCF_FINAL_PVM} -P $pathology


if [ "$?" = "0" ]; then
	echo -e  '\nEXIT STATUS: 0'
	echo -e  '\nPOST-VEP MODIFICATIONS for '${run}' DONE'
else
	echo -e  "ERROR: PROBLEMS WITH POST-VEP MODIFICATIONS"
	exit 1
fi

end=`date +%s`
runtime=$((end-start))
echo -e  '\nExecuting time: '$runtime 




# Filter PVM files based on gene list

if [ "$genefilter" != "False" ]; then


	start=`date +%s`

	VCF_FINAL_PVM_FILTER="${VEPVCFAD}/${run}_filteredAnnotated_pvm_GENELIST.txt"
	
	python $utilitiesPath/filtering_geneList.py -i ${VCF_FINAL_PVM} -f ${genefilter} -o ${VCF_FINAL_PVM_FILTER}


	if [ "$?" = "0" ]; then
		printf '\nEXIT STATUS: 0'
		printf '\nVARIANT FILTERING BY GENE LIST for '${run}' DONE'

	else
		printf "\nERROR: PROBLEMS WITH GENE FILTERING"
		exit 1
	fi

	end=`date +%s`
	runtime=$((end-start))
	printf '\nExecuting time: '$runtime 

fi







# Remove files 

echo -e '\n'
rm -r $TMP











#### VCF PROCESSING

# if [ "$ped" != "null" ]; then


# 	gatk VariantsToTable \
# 	-R $HG19/ucsc.hg19.fasta \
# 	-V $VCF_FILTERED_3 \
# 	-F CHROM -F POS -F TYPE -F REF -F ALT -F HET -F HOM-VAR  -F HOM-REF -GF GT -GF GQ -F hiConfDeNovo -F lowConfDeNovo \
# 	--split-multi-allelic -O $MDAP/SNV_filtAnnot.txt --add-output-vcf-command-line

# else


# 	gatk VariantsToTable \
# 	-R $HG19/ucsc.hg19.fasta \
# 	-V $VCF_FILTERED_3 \
# 	-F CHROM -F POS -F TYPE -F REF -F ALT -F HET -F HOM-VAR  -F HOM-REF -GF GT -GF GQ -F hiConfDeNovo -F lowConfDeNovo \
# 	--split-multi-allelic -O $MDAP/SNV_filtAnnot.txt --add-output-vcf-command-line

# fi






# Trio phasing

#  This tool performs two functions:

#     Compute the most likely genotype combination of trios and parent/child pairs given their genotype likelihoods and a mutation prior;

# The tool ultimately reports the genotype combination (and hence phasing) probability.
# Ambiguous sites are:

#     Sites where all individuals are heterozygous
#     Sites where there is a Mendelian violation

# Missing genotypes are handled as follows:

#     In parent/child pairs: If an individual genotype is missing at one site, the other one is phased if it is homozygous. No phasing probability is emitted.
#     In trios: If the child is missing, parents are treated as separate individuals and phased if homozygous. No phasing probability is emitted.
#     In trios: If one of the parents is missing, it is handled like a parent/child pair. Phasing is done unless both the parent and child are heterozygous and a phasing probability is emitted.
#     In trios: If two individuals are missing, the remaining individual is phased if it is homozygous. No phasing probability is emitted.




