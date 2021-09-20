#!/bin/bash

module load python/2.7.15
module load perl
source ~/.Renviron
module load R

Rscript /home/proyectos/bioinfo/lodela/VariantCallingFJD/utilities/codex_normalization.R


