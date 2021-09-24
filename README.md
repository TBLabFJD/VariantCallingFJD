# variantCallingFJD





The objective of the development of this pipeline was to automate the customized genomics analysis that we carry out in the **Bioinformatics Unit** for the **Department of Genetics and Genomics** of the [**Institituto de Investigación Sanitaria Fundación Jiménez Díaz (IIS-FJD)**](https://www.fjd.es/iis-fjd). This pipeline is designed to be run in a [Slurm Workload Manager](https://slurm.schedmd.com/documentation.html) system wich is an "open source, fault-tolerant, and highly scalable cluster management and job scheduling system for large and small Linux clusters".

This pipeline has been developed by [**Translational Bioinformatics Lab**](https://www.translationalbioinformaticslab.es/tblab-home-page) at the [**IIS-FJD**](https://www.fjd.es/iis-fjd). 

## License
VariantCallingFJD source code is provided under the [**Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)**](https://creativecommons.org/licenses/by-nc-sa/4.0/). VariantCallingFJD includes several third party packages provided under other open source licenses, please see COPYRIGHT.txt for additional details.

[![Licencia de Creative Commons](https://i.creativecommons.org/l/by-nc-sa/4.0/88x31.png)](http://creativecommons.org/licenses/by-nc-sa/4.0/)

## Requirements
- Slurm Workload Manager
- python/2.7.15
- perl/5.28.0
- bwa/0.7.17
- samtools/1.9
- picard/2.18.9
- R/3.5.0
- miniconda/3.6
- gatk/4.2.0
- samtools/1.9
- bedtools/2.27.0
- R/R
- annotsv/2.2
- PLINK v1.90b6.9 64-bit (4 Mar 2019)

R libraries
- CODEX2_1.3.0
- panelcn.mops_1.4.0
- cn.mops_1.28.0
- ExomeDepth_1.1.15

Python libraries


## Instalation
1. Clone this repository using:
```sh
git clone https://github.com/TBLabFJD/VariantCallingFJD.git
```
2. There is a CONFIG_FILE.txt that needs to be filled with the apropiate paths for:
   - Reference genome (in FASTA format)
   - Binaries
   - Data bases and tables for annotation
   - Slurm credentials
We recommed installing mosdepth with conda

## Getting Started

 
 




