#!/bin/bash
#SBATCH --account=bioinfo_serv
#SBATCH --partition=bioinfo
#SBATCH --job-name=mutect2   #job name
#SBATCH --mail-type=FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ldelafuente.lorena@gmail.com # Where to send mail        
#SBATCH --mem-per-cpu=2gb # Per processor memory
#SBATCH --cpus-per-task=5
#SBATCH -t 15:00:00     # Walltime
#SBATCH -o mutect2_%j.out # Name output file 
##SBATCH --error=
##SBATCH --file=
##SBATCH --initaldir=


sample=$1
input=$2
output=$3
bedfile=$4

echo $sample


module load miniconda/2.7
module load bwa/0.7.15
module load samtools/1.9
module load picard/2.18.9
module load gatk/4.1.2.0
module load vep/release98
module load bedtools/2.27.0
module load R
#module load plink
alias picard='java -jar /usr/local/bioinfo/picard-tools/2.18.9/picard.jar'
alias gatk='java -jar /usr/local/bioinfo/gatk/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar'

# mutect2

gnomad="/home/proyectos/bioinfo/references/VEPdbs/af-only-gnomad.raw.sites.hg19.vcf.gz"
reffasta="/home/proyectos/bioinfo/references/hg19/ucsc.hg19.fasta"
tumorbam="$input/${sample}_alignment.bam"
echo $tumorbam

unfilteredvcf="$output/${sample}_unfiltered.vcf"
assembly="$output/${sample}_assemblyIGV.txt"
profile="$output/${sample}_profile.txt"
f1r2="$output/${sample}_f1r2.tar.gz"
readmodel="$output/${sample}_read-orientation-model.tar.gz"
filteredvcf="$output/${sample}_filtered.vcf"
bamout="$output/${sample}_bamout.bam"
pileupsummaries="$output/${sample}_getpileupsummaries.table"



### Call Somatic SNV & INDELS with MUTECT2

gatk Mutect2 \
-R $reffasta \
-I $tumorbam \
-O $unfilteredvcf \
-bamout $bamout \
--activity-profile-out $profile -A AlleleFraction -A AS_FisherStrand -A StrandBiasBySample -A PossibleDeNovo -A DepthPerAlleleBySample \
--interval-padding 2000 --intervals $bedfile --assembly-region-out $assembly \
-germline-resource $gnomad \
--f1r2-tar-gz $f1r2



gatk LearnReadOrientationModel -I $f1r2 -O $readmodel

gatk FilterMutectCalls -R $reffasta -V $unfilteredvcf -O $filteredvcf --stats ${unfilteredvcf}.stats --ob-priors $readmodel









## RESOURCE: GATK4

# optimization!!!!
# for chromosome in {1..22}; do
# gatk Mutect2 -R ref.fasta -I tumor.bam -L $chromosome --f1r2-tar-gz ${chromosome}-f1r2.tar.gz -O ${chromosome}-unfiltered.vcf``
# done
# all_f1r2_input=`for chromosome in {1..22}; do printf -- "-I ${chromosome}-f1r2.tar.gz "; done`
# LearnReadOrientationModel $all_f1_r2_input -O read-orientation-model.tar.gz



 # Furthermore, we never liked how germline variants, which are handled in a more principled way with our germline filter, ended up as hard-filtered pon sites, so the panel of normals workflow now optionally excludes germline events from its output, keeping only technical artifacts.
