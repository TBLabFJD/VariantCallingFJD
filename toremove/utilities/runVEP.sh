#!/bin/bash

#SBATCH --account=bioinfo_serv
#SBATCH --partition=bioinfo
#SBATCH --job-name=vep   #job name
#SBATCH --mail-type=FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ldelafuente.lorena@gmail.com # Where to send mail        
#SBATCH --mem-per-cpu=15gb # Per processor memory
#SBATCH --cpus-per-task=4
#SBATCH -t 15:00:00     # Walltime
#SBATCH -o vep_%j.out # Name output file 
#SBATCH --error=vep_%j.err 




MDAP=$1  # output folder
sample=$2
local=${3}
pathology=${4}
utilitiesPath=${5}

bash $utilitiesPath/vepAnnotation.sh $MDAP $sample $local $pathology $utilitiesPath 4
