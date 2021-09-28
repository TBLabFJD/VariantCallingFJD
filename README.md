# VariantCallingFJD


The objective of the development of this pipeline was to automate the customized genomics analysis that we carry out in the **Bioinformatics Unit** for the **Department of Genetics and Genomics** of the [**Institituto de Investigación Sanitaria Fundación Jiménez Díaz (IIS-FJD)**](https://www.fjd.es/iis-fjd). This pipeline is designed to be run in a [Slurm Workload Manager](https://slurm.schedmd.com/documentation.html) system wich is an "open source, fault-tolerant, and highly scalable cluster management and job scheduling system for large and small Linux clusters".

This pipeline has been developed by [**Translational Bioinformatics Lab**](https://www.translationalbioinformaticslab.es/tblab-home-page) at the [**IIS-FJD**](https://www.fjd.es/iis-fjd). 

## License
VariantCallingFJD source code is provided under the [**Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)**](https://creativecommons.org/licenses/by-nc-sa/4.0/). VariantCallingFJD includes several third party packages provided under other open source licenses, please see COPYRIGHT.txt for additional details.


[![Licencia de Creative Commons](https://i.creativecommons.org/l/by-nc-sa/4.0/88x31.png)](http://creativecommons.org/licenses/by-nc-sa/4.0/)

## Dependencies
**Job scheduler:**
- **Slurm Workload Manager:** As said in the introduction, this pipeline has been developped to be run in a Slurm Workload Manager system. We use the the function `sbatch` to submit the jobs to the queue. We also have installed the diferent programs in different modules to prevent incompatibilities.

**Programming languages:**
- **Python v2.7.15:** All python scrips were developped using python v2.7.15. In task `"Mosdepth.sh"` python v3.6.12 is used by loading the module miniconda/3.6 to use `mosdepth`.
- **Perl v5.28.0**
- **R v3.5.0**

**Bioinformatic tools:**
- **bwa v0.7.17**
- **samtools v1.9**
- **picard v2.18.9**
- **gatk v4.2.0**
- **mosdepth 0.2.5**
- **bedtools v2.27.0** 
- **bcftools v1.3**
- **annotsv v2.2**
- **PLINK v1.90b6.9 64-bit (4 Mar 2019)**
- **BaseSpaceCLI v1.0.0**
- **bscp v0.6.1.337**

**R libraries:**
- **dplyr v0.8.3**
- **optparser v1.6.6**
- **stringr v1.4.0**
- **CODEX2 v1.3.0**
- **panelcn.mops v1.4.0**
- **cn.mops v1.28.0**
- **ExomeDepth v1.1.15**
- **GenomicRanges v1.34.0**


**Python libraries**
- **csv v1.0**
- **argparse v1.1**
- **os**
- **subprocess**
- **sys**
- **re v2.2.1**
- **glob**
- **numpy v1.16.6**
- **pandas v0.24.2**
- **datetime**
- **shutil**
- **ConfigParser**
- **string**
- **json v2.0.9**
- **itertools**
- **operator**
- **collections**
- **time**

## Instalation
1. Install all the dependencies.
    1.1. We recommend installing all the above bioinformatics programs and programming languages in modules as follow:
    - Python v2.7.15 --> python/2.7.15
    - Perl v5.28.0 --> perl/5.28.0
    - R v3.5.0 --> R/R
    - bwa v0.7.17 --> bwa/0.7.17
    - samtools v1.9 --> samtools/1.9
    - picard v2.18.9 --> picard/2.18.9
    - gatk v4.2.0 --> gatk/4.2.0
    - mosdepth 0.2.5 --> miniconda/3.6
    - bedtools v2.27.0 --> bedtools/2.27.0
    - bcftools v1.3 --> bcftools/1.3
    - annotsv v2.2 --> annotsv/2.2
    - vep release 103 --> vep/release103
    
    1.2. Install all Python and R packages
    1.3 Download the following binaries and executables
    - The BaseSpace Sequence Hub CLI tool suite binary
    - BaseSpace-copy binary
    - GATK executable jar file
    - Picard executable jar file
    - PLINK binary

2. Clone this repository using:
```sh
git clone https://github.com/TBLabFJD/VariantCallingFJD.git
```



## Getting Started
0. (Optional) **BaseSpace credentials setup.** To be able to use the -b/--basespace option to automatically download samples from BaseSpace, run this command:
```sh 
BaseSpace_Sequence_Hub_CLI_tool_suite_binary_path/bs authenticate
```
 * This command generates a URL to copy into a web browser and login into BaseSpace. After login, it creates into your directory a `.basespace` file with the credentials so that when you run the pipeline it autamatically acces the data.

1. **Configuration file.** There is a [`pipeline.conf`](https://github.com/TBLabFJD/VariantCallingFJD/blob/master/pipeline.conf) file that needs to be filled with the apropiate information (some information is optional or required depending on the selected parameters):
   - Slurm credentials
     - account
     - partition
   - Files
     - Reference genome (in FASTA format)
     - Validated VCF with known sites of common variation. We use a VCF from The HapMap project. Available at: https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz
     - Validated VCF with known sites of common variation. We use a VCF from The 1000 Genomes Project. Available at: https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz
     - Index of the reference (.fasta.fai) genome to be used by CoNVading (CNV calling algorithm). This file should be a two column file with the name of the chromosomes (without the "chr") and their size (See CoNVaDING documentation for more information).
   - Directories
     - Temporal directory used for merging FASTQ files.
     - (Optional) Dictory where to move a copy of the VCF file with the SNVs and INDELs for data base creation.
     - (Optional) Dictory where to move a copy of the BED file with the covered regions for data base creation.
   - Binaries & Executables
     - The BaseSpace Sequence Hub CLI tool suite binary
     - BaseSpace-copy binary
     - GATK executable (.jar)
     - PICARD executable (.jar)
     - conda binary (same binary as the une installed for the miniconda)
     - PLINK binary
     - ANNOTSV executable (.tcl)
 2. Run the pipeline as follow
 ```sh
module load python/2.7.15
module load perl/5.28.0

pipeline="path_to_the_downloaded_git_repository/variantDiscoveryFJD_panelWES.py"
input_path="input_directory_with_fastq_files"
output_path="output_directory"
bed_path="path_to_bed_file"

python ${pipeline} -i ${input_path} -o ${output_path} -p ${bed_path} -a all -A -t 20 -M 20 -P healthy
```
 