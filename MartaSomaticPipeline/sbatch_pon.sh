#!/bin/bash
#SBATCH --account=bioinfo_serv
#SBATCH --partition=bioinfo
#SBATCH --job-name=PON   #job name
#SBATCH --mail-type=FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=martarodmo@gmail.com # Where to send mail
#SBATCH --mem-per-cpu=2gb # Per processor memory
#SBATCH --cpus-per-task=5
#SBATCH -t 15:00:00     # Walltime
#SBATCH -o mutect2_%j.out # Name output file
##SBATCH --error=
##SBATCH --file=
##SBATCH --initaldir=

module load miniconda/2.7
module load bwa/0.7.15
module load samtools/1.9
#module load pird/2.18.9
module load gatk/4.1.5.0
module load vep/release98
module load bedtools/2.27.0
module load R

reference="/home/proyectos/bioinfo/references/hg19/ucsc.hg19.fasta"
intervals="/home/proyectos/bioinfo/lodela/BioinfoUnit/martaLymp/3110861_Covered.bed"
gnomad="/home/proyectos/bioinfo/references/VEPdbs/af-only-gnomad.raw.sites.hg19.vcf.gz"

gatk Mutect2 -R $reference -I /home/proyectos/bioinfo/mrod/newResults_Oct20/nextseq_4/bams/B14-28162-2_alignment.bam -max-mnp-distance 0 -L $intervals \
	-O /home/proyectos/bioinfo/mrod/newResults_Oct20/pon_results/mutect2/B14-28162-2_unfiltered.vcf

gatk Mutect2 -R $reference -I /home/proyectos/bioinfo/mrod/newResults_Oct20/nextseq_4/bams/B14-34524-2_alignment.bam -max-mnp-distance 0 -L $intervals \
	-O /home/proyectos/bioinfo/mrod/newResults_Oct20/pon_results/mutect2/B14-34524-2_unfiltered.vcf

gatk Mutect2 -R $reference -I /home/proyectos/bioinfo/mrod/newResults_Oct20/nextseq_4/bams/B14-28669-5_alignment.bam -max-mnp-distance 0 -L $intervals \
	-O /home/proyectos/bioinfo/mrod/newResults_Oct20/pon_results/mutect2/B14-28669-5_unfiltered.vcf




# Create Somatic PanelOfNormals in DB

gatk GenomicsDBImport -R $reference -L $intervals \
       --genomicsdb-workspace-path pon_db \
       -V /home/proyectos/bioinfo/mrod/newResults_Oct20/pon_results/mutect2/B14-28162-2_unfiltered.vcf \
       -V /home/proyectos/bioinfo/mrod/newResults_Oct20/pon_results/mutect2/B14-34524-2_unfiltered.vcf \
       -V /home/proyectos/bioinfo/mrod/newResults_Oct20/pon_results/mutect2/B14-28669-5_unfiltered.vcf

gatk CreateSomaticPanelOfNormals -R $reference -V gendb://pon_db -O pon_ganglios_db.vcf.gz --germline-resource $gnomad 
