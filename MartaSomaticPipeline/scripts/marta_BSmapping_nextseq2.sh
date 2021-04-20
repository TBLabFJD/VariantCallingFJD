#!/bin/bash

module load miniconda/2.7
module load perl
source ~/.Renviron

project="/scratch/lodela/DS-Linfomas-FJD-3_2020_06_29_11_31_26/B20314403"
#project="/scratch/lodela/Linfomas_NextSeq500_23-09-2019"
#project="/scratch/lodela/Linfomas_NextSeq500_200420_demultiplexed_newManifest/Linfomas_NextSeq500_20-01-2020"
output="/home/proyectos/bioinfo/lodela/BioinfoUnit/martaLymp/newResults_Nov19/nextseq_round4"
pipeline="/home/proyectos/bioinfo/fjd/VariantCallingFJD/pipelineFJD19.py"
#samples="/home/proyectos/bioinfo/lodela/martaLymp/results_nextseq500_2/samples.txt"
samples="/home/proyectos/bioinfo/lodela/BioinfoUnit/martaLymp/newResults_Nov19/nextseq_round4/samples.txt"

python $pipeline -q -u martaLymp -i $project -o ${output} -s $samples -n Nextseq -a mapping -t 5 -M 2 -m ldelafuente.lorena@gmail.com

#python /home/proyectos/bioinfo/pipelineFJD2019v3/pipelineFJD19_bsNwait.py -i $project -o $output -a mapping -t 4 -b -u martaLymp

