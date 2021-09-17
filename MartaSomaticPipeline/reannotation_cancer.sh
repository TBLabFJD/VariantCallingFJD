#!/bin/bash

sophia_dir=/home/proyectos/bioinfo/NOBACKUP/mrod/mutect2_20210623

ls ${sophia_dir}/*_unfiltered.vcf | while read f ; do


			
	sample=$(echo $f | cut -d"/" -f8 )
	dna=${sample%_u*}
	echo $dna

	sbatch --account=bioinfo_serv --partition=bioinfo --job-name=filter --mem-per-cpu=5gb --cpus-per-task=12  ./vepAnnotation_raquel.sh ${sophia_dir} ${dna} False healthy $sophia_dir 12


done



# sbatch --account=bioinfo_serv --partition=bioinfo --job-name=vep --mem-per-cpu=10gb --cpus-per-task=4 


# MDAP=$1  # output folder
# sample=$2
# local=${3}
# pathology=${4}
# utilitiesPath=${5}
# threads=${6}
