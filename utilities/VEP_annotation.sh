#!/bin/sh

##################################
### FJD pipeline -  ANNOTATION ###
##################################

# ARGUMENTS


local=$1
MDAP=$2
VEPVCFAD=${MDAP}/snv_results
name=$3 # sample name or run name
run=$4
threads=$5



# REQUIRED MODULES AND DATABASES:

if [ "$local" != "True" ]; then


	module load vep/release95

	VEP="/usr/local/bioinfo/vep/ensembl-vep/vep"
	FILTER_VEP='/usr/local/bioinfo/vep/ensembl-vep/filter_vep'
	VEP_CACHE='/usr/local/bioinfo/vep/ensembl-vep/t/testdata/cache/homo_sapiens'
	VEP_FASTA="/home/proyectos/bioinfo/references/VEPfasta/Homo_sapiens.GRCh37.dna.primary_assembly.fa"
	PLUGIN_DIR=/usr/local/bioinfo/vep/ensembl-vep/plugins-95/VEP_plugins-release-95
	PLUGIN_DBS="/home/proyectos/bioinfo/references/VEPdbs"
	dbNSFP_DB="${PLUGIN_DBS}/dbNSFP3.5a_hg19.gz"
	CCS_DB="/home/proyectos/bioinfo/references/CCS/ccrs.autosomes.v2.20180420.bed.gz"

	softwareFile="${MDAP}/software_${run}.txt"
	title="MAPPING"
	if [ ! -f $softwareFile ] || [ `grep -q $title $softwareFile` ] ; then 

		printf "VEP ANNOTATION:\n" >> ${softwareFile}
		module list 2>> ${softwareFile}
	
	fi




else

	export SFT=/mnt/genetica3/marius/pipeline_practicas_marius/software
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

		printf "VEP ANNOTATION:\n" >> ${softwareFile}

		printf "\nVEP VERSION\n" >> ${softwareFile}
		$VEP --help | grep "Versions:" -A 5 | tail -n4 >> ${softwareFile}

	fi




fi




VCF_IN="${VEPVCFAD}/${name}_raw.vcf"
VCF_OUT="${VEPVCFAD}/${name}_filtered.vcf"
VCF_FILTERED="${VEPVCFAD}/${name}_annotated.vcf"
VCF_FILTERED_2="${VEPVCFAD}/${name}_filtered_annotated_can_conseq.vcf"
VCF_FILTERED_3="${VEPVCFAD}/${name}_filteredAnnotated.vcf" 





# COMMANDS







printf '\n\nFiltering by quality and chromosome...'



start=`date +%s`

perl $FILTER_VEP \
-i $VCF_IN -o $VCF_OUT \
--filter "QUAL > 100 and FILTER = PASS" \
--filter "CHROM in chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chr23,chrX,chrY" --force_overwrite 




if [ "$?" = "0" ]; then
	printf '\nEXIT STATUS: 0'
	printf '\nVEP FILTERING quality and chromosome for '${name}' DONE\n' 

else
	printf "\nERROR: PROBLEMS WITH VEP FILTERING FOR QUALITY AND CHR"
	exit 1
fi














printf '\n\nVEP annotation...'

perl $VEP \
--cache --offline --hgvs --refseq --dir $VEP_CACHE --dir_plugins $PLUGIN_DIR --v --fork $threads --assembly GRCh37 --fasta $VEP_FASTA --force_overwrite \
--biotype --regulatory --protein --symbol --allele_number --numbers --domains --uniprot --variant_class --no_stats \
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
	printf '\nVEP ANNOTATION for '${name}' DONE\n' 
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
	printf '\n VEP FREQUENCY FILTERING 1 for '${name}' DONE\n' 
	rm $VCF_FILTERED

else
	printf "\nERROR: PROBLEMS WITH VEP FREQUENCY FILTERING 1"
	exit 1
fi







perl $FILTER_VEP \
-i $VCF_FILTERED_2 -o $VCF_FILTERED_3 \
--filter "AF < 0.01 or not AF" \
--filter "(ExAC_EAS_AF < 0.01 or not ExAC_EAS_AF) and (1000Gp3_AF < 0.01 or not 1000Gp3_AF) and (1000Gp3_EUR_AF < 0.01 or not 1000Gp3_EUR_AF ) and (gnomAD_exomes_AF < 0.01 or not gnomAD_exomes_AF) and (ExAC_AF < 0.01 or not ExAC_AF) and (ExAC_Adj_AF < 0.01 or not ExAC_Adj_AF) and (ExAC_NFE_AF < 0.01 or not ExAC_NFE_AF) and (gnomAD_exomes_NFE_AF < 0.01 or not gnomAD_exomes_NFE_AF) and (gnomAD_genomes_AF < 0.01 or not gnomAD_genomes_AF)" --force_overwrite 


if [ "$?" = "0" ]; then
	printf '\nEXIT STATUS: 0'
	printf '\nVEP FREQUENCY FILTERING 2 for '${name}' DONE\n' 
	rm $VCF_FILTERED_2

else
	printf "\nERROR: PROBLEMS WITH VEP FREQUENCY FILTERING 2"
	exit 1
fi






printf '\nVEP annotation and filtering DONE.\n'


end=`date +%s`
runtime=$((end-start))
printf '\nExecuting time: '$runtime 
















