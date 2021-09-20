#!/bin/bash
#SBATCH --account=bioinfo_serv
#SBATCH --partition=bioinfo
#SBATCH --job-name=mutect2   #job name
#SBATCH --mail-type=FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=martarodmo@gmail.com # Where to send mail        
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
module load pird/2.18.9
module load gatk/4.1.5.0
module load vep/release98
module load bedtools/2.27.0
module load R


alias picard='java -jar /usr/local/bioinfo/picard-tools/2.18.9/picard.jar'
alias gatk='java -jar /usr/local/bioinfo/gatk/gatk-4.1.5.0/gatk-package-4.1.5.0-local.jar'

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
calculatecontaminationtable="$output/${sample}_calculatecontamination.table"



### Call Somatic SNV & INDELS with MUTECT2

gatk Mutect2 \
-R $reffasta \
-I $tumorbam \
-O $unfilteredvcf \
-bamout $bamout \
--intervals $bedfile --assembly-region-out $assembly \
-germline-resource $gnomad \
-pon /home/proyectos/bioinfo/mrod/newResults_Oct20/pon_results/pon_ganglios_db.vcf.gz  \
--f1r2-tar-gz $f1r2



### LearnReadOrientationModel

gatk LearnReadOrientationModel -I $f1r2 -O $readmodel   #### EL OUTPUT NO HEMOS REPETIDO SIEMPRE QUE DEBE TENER EL NOMBRE DE LA MUESTRA PORQUE SI LANZAS EN PARALELO SE SOBREESCRIBEN??!



### Run GetPileupSummaries to summarize read support for a set number of known variant sites



# # create intervals before running this script

# java -jar /usr/local/bioinfo/picard-tools/2.18.9/picard.jar BedToIntervalList \
#       I=/home/proyectos/bioinfo/lodela/BioinfoUnit/martaLymp/3110861_Covered.bed \
#       O=/home/proyectos/bioinfo/mrod/scriptsToMarta/list.interval_list \
#       SD=/home/proyectos/bioinfo/mrod/newResults_Oct20/pon_results/pon_ganglios_db.vcf.gz




gatk GetPileupSummaries -I $tumorbam -V $gnomad -L /home/proyectos/bioinfo/mrod/scriptsToMarta/list.interval_list -O $pileupsummaries   ############# LO MISMO CON EL OUTPUT, LO HE CAMBIADO YO PERO DEBERÍAS HABERLO CAMBIADO




### CalculateContamination

gatk CalculateContamination -I $pileupsummaries --tumor-segmentation segments.table -O $calculatecontaminationtable



### FilterMutectCallswith the -ob-priors argument
gatk FilterMutectCalls -R $reffasta -V $unfilteredvcf --tumor-segmentation segments.table --contamination-table $calculatecontaminationtable --ob-priors $readmodel -O $filtered.vcf




















## RESOURCE: GATK4

# optimization!!!!

#for chromosome in {1..22}; do
#gatk Mutect2 -R ref.fasta -I tumor.bam -L $chromosome --f1r2-tar-gz ${chromosome}-f1r2.tar.gz -O ${chromosome}-unfiltered.vcf``
#done
#all_f1r2_input=`for chromosome in {1..22}; do printf -- "-I ${chromosome}-f1r2.tar.gz "; done`
#LearnReadOrientationModel $all_f1_r2_input -O read-orientation-model.tar.gz




 # Furthermore, we never liked how germline variants, which are handled in a more principled way with our germline filter, ended up as hard-filtered pon sites, so the panel of normals workflow now optionally excludes germline events from its output, keeping only technical artifacts.
