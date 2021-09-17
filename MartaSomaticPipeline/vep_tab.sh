#"${VCF_IN}_annotated.tab""${VCF_IN}_annotated.tab"#!/bin/bash
#SBATCH --account=bioinfo_serv
#SBATCH --partition=bioinfo
#SBATCH --job-name=vep   #job name
#SBATCH --mail-type=FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ldelafuente.lorena@gmail.com # Where to send mail        
#SBATCH --mem-per-cpu=5gb # Per processor memory
#SBATCH --cpus-per-task=3
#SBATCH -t 15:00:00     # Walltime
#SBATCH -o tophat2_%j.out # Name output file 
##SBATCH --error=
##SBATCH --file=
##SBATCH --initaldir=


VCF_IN="B19-303205_unfiltered.vcf"
VCF_OUT="${VCF_IN}_annotated.tab"


module load python/2.7.15
module load perl
source ~/.Renviron
module load perl 
#module load miniconda/2.7
module load bwa/0.7.15
module load samtools/1.9
module load picard/2.18.9
module load gatk/4.1.2.0
module load vep/release98
module load bedtools/2.27.0
module load R
alias picard='java -jar /usr/local/bioinfo/picard-tools/2.18.9/picard.jar'
alias gatk='java -jar /usr/local/bioinfo/gatk/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar'

VEP="/usr/local/bioinfo/vep/ensembl-vep/vep"
FILTER_VEP='/usr/local/bioinfo/vep/ensembl-vep/filter_vep'
VEP_CACHE='/usr/local/bioinfo/vep/ensembl-vep/t/testdata/cache/homo_sapiens'
VEP_FASTA="/home/proyectos/bioinfo/references/VEPfasta/Homo_sapiens.GRCh37.dna.primary_assembly.fa"
PLUGIN_DIR='/usr/local/bioinfo/vep/ensembl-vep/Plugins'
PLUGIN_DBS="/home/proyectos/bioinfo/references/VEPdbs"
dbNSFP_DB="${PLUGIN_DBS}/dbNSFP3.5a_hg19.gz"
CCS_DB="/home/proyectos/bioinfo/references/CCS/ccrs.autosomes.v2.20180420.bed.gz"
utilitiesPath="/home/proyectos/bioinfo/fjd/VariantCallingFJD/utilities"


printf '\n\nVEP annotation...'


perl $VEP \
--cache --offline --hgvs --refseq --dir $VEP_CACHE --dir_plugins $PLUGIN_DIR --v --fork 3 --assembly GRCh37 --fasta $VEP_FASTA --force_overwrite \
--biotype --regulatory --protein --symbol --allele_number --numbers --domains --uniprot --variant_class \
--canonical --tab \
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
-i $VCF_IN -o $VCF_OUT
