#!/bin/sh

################################
### FJD pipeline - Haplotype ###
################################

### FJD PIPELINE ARGUMENTS:



MDAP=$1  # output folder
sample=$2
local=${3}
pathology=${4}
utilitiesPath=${5}
threads=${6}

printf "\n......................................................\n"
printf "  MAPPING OR/AND SNV CALLING FOR SAMPLE $sample \n"
printf "......................................................\n"


# REQUIRED MODULES AND DATABASES:

if [ "$local" != "True" ]; then

	module load perl
	module load miniconda/2.7
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










































printf '\n\nVEP annotation...'

perl $VEP \
--cache --offline --hgvs --refseq --dir $VEP_CACHE --dir_plugins $PLUGIN_DIR --v --fork $threads --assembly GRCh37 --fasta $VEP_FASTA --force_overwrite \
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
-i $VCF_FILTERED -o $VCF_FILTERED_2 \
--filter "DP > 10" \
--filter "CANONICAL is YES" \
--filter "Consequence is 3_prime_UTR_variant or Consequence is 5_prime_UTR_variant or Consequence is intron_variant or Consequence is splice_donor_variant or Consequence is splice_acceptor_variant or Consequence is splice_region_variant or Consequence is synonymous_variant or Consequence is missense_variant or Consequence is inframe_deletion or Consequence is inframe_insertion or Consequence is stop_gained or Consequence is frameshift_variant or Consequence is coding_sequence_variant" --force_overwrite 


if [ "$?" = "0" ]; then
	printf '\nEXIT STATUS: 0'
	printf '\n VEP FREQUENCY FILTERING 1 for '${sample}' DONE\n' 
	rm $VCF_FILTERED

else
	printf "\nERROR: PROBLEMS WITH VEP FREQUENCY FILTERING 1"
	exit 1
fi







# perl $FILTER_VEP \
# -i $VCF_FILTERED_2 -o $VCF_FILTERED_3 \
# --filter "AF < 0.01 or not AF" \
# --filter "(ExAC_EAS_AF < 0.01 or not ExAC_EAS_AF) and (1000Gp3_AF < 0.01 or not 1000Gp3_AF) and (1000Gp3_EUR_AF < 0.01 or not 1000Gp3_EUR_AF ) and (gnomAD_exomes_AF < 0.01 or not gnomAD_exomes_AF) and (ExAC_AF < 0.01 or not ExAC_AF) and (ExAC_Adj_AF < 0.01 or not ExAC_Adj_AF) and (ExAC_NFE_AF < 0.01 or not ExAC_NFE_AF) and (gnomAD_exomes_NFE_AF < 0.01 or not gnomAD_exomes_NFE_AF) and (gnomAD_genomes_AF < 0.01 or not gnomAD_genomes_AF)" --force_overwrite 


# if [ "$?" = "0" ]; then
# 	printf '\nEXIT STATUS: 0'
# 	printf '\nVEP FREQUENCY FILTERING 2 for '${sample}' DONE\n' 
# 	rm $VCF_FILTERED_2

# else
# 	printf "\nERROR: PROBLEMS WITH VEP FREQUENCY FILTERING 2"
# 	exit 1
# fi






# printf '\nVEP annotation and filtering DONE.\n'


# end=`date +%s`
# runtime=$((end-start))
# printf '\nExecuting time: '$runtime 






# printf "\n\n\n- OUTPUT PROCESSING "
# printf "\n------------------------\n"

# start=`date +%s`

# printf '\nFrom VEP annotated VCF to txt file...\n'
# start=`date +%s`


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
fi

end=`date +%s`
runtime=$((end-start))
printf '\nExecuting time: '$runtime 





# Remove temporary files.

printf '\n'
rm -r $TMP








# echo -e  "\n\n\n- HOMOZYGOSITY (PLINK) "
# echo -e  "-------------------------\n"


# 	plink --vcf -recode -out mycancerplink2
# 	../plink --file hapmap1 --homozyg  --homozyg-snp 50 --homozyg-window-het 0
# 	vcftools --vcf genotyped_data_allsamples.vcf --out mycancerplink2 --plink
# 	../plink --bfile hola-temporary --homozyg  --homozyg-snp 50  --homozyg-window-het 0

# 	# Plink default: 
# 	--homozyg-kb 1,000
# 	--homozyg-snp 100

# 	# Plink improvemements (minimizes the trade-off between the exclusion of non-autozygous ROHs and the cost of missing shorter autozygous ROHs)
# 	# - 1. remove SNPs with MAF < 0.05
# 	# - 2. LD prunning in two levels: ligth to moderate outperforms ROH detection -> --indep 
# 	# Parameters of light prunning: r2>0.9 in a 50 SNP window (VIF > 10)
# 	# Parameters of moderate prunning: 50 SNP window that had r2>0.5 (VIF > 2)

# 	# - type1 and type2 error rates are computed at SNP level.

# 	# You may also want to try 'bcftools roh', which uses a HMM-based detection method. 
# 	#(We'll include a basic port of that command in PLINK 2.0 if there is sufficient interest.)

# 	#(For more accurate detection of smaller segments, one might consider approaches that also take
# 	# population parameters such as allele frequency and recombination rate into account, 
# 	# in a HMM approach for example: but for now, PLINK only supports this basic detection of long, homozygous segments). 

# 	#   The originality of this algorithm consists in its capabil-ity to incorporate distances between consecutive SNPs to discriminate between the homozygosity and the hetero-zygosity states. This feature of  H3M2  offers advantage in the detection of small-sized (<500 kb) and medium-sized (between 500 and 1,500 kb) ROH  [18] . 
# 	# Differently,  H3M2and PLINK show similar performances in the WES-based detection of long ROH (>1,500 kb)

# fi


