This pipeline performs quality control, variant calling and variant filtering in a secuencial way as follows:

[mosdepth](https://platform.batchx.io/iis-fjd/tools/mosdepth) --> [snv-calling](https://platform.batchx.io/iis-fjd/tools/snv-calling) --> [snv-filtering](https://platform.batchx.io/iis-fjd/tools/snv-filtering)

# Inputs
## Required inputs
This tool has the following **required** inputs:
1. `refFolder`
Reference folder containing the reference genome in FASTA format, a dictionary file ending in `.dict` and an index file ending in `.fai`. The dictionary can be created using the command `gatk-launch CreateSequenceDictionary -R ref.fasta` where `ref.fasta` is the reference sequence. The index file can be created using the command `samtools faidx ref.fasta`. These two files are used by gatk HaplotypeCaller (and other tools) to efficiently access records of the genome. The reference must be the same as the one used to create de BAM file.
2. `bamFile`
File containing the mapped reads with the extension `.bam`. A BAM file (*.bam) is the compressed binary version of a SAM file that is used to represent aligned sequences. 
3. `bamIndex`
Index file of the bamfile with the extension `.bai` or `.bam.bai`. It stores the genomic positions to efficiently access records of the BAM file. The index file can be created using the command `samtools faidx sample.bam` where `sample.bam` is the BAM file from which we want to ceate the index.
4. `sample`
Sample ID.


## Optional inputs
1. `minCoverage` 
mosdepth minimum read depth.
1. `gvcfCalling` 
Used to perform snv calling in `GVCF` mode.
2. `panelFile`
BED file containing the captured regions of the genome based on the library design. It is a TSV file without a header. It has a minimum of 3 columns storing the chromosome, start and stop position for each target. For CNV calling, it needs a 4th column containing the gene name of the target. In this step, the `panelFile` is used to calculate the coverage distribution among the theoretical covered regions. If not provided, the coverage is calculated for all the genome.
3. `resources` 
Memory and CPUs to assign to each pipeline step.

# Outputs
1. `vcf`
- VCF (Variant Call Format) file containing all called SNVs and INDELs. (default)
- GVCF (Genomic Variant Call Format) file containing  all called SNVs and INDELs and all covered regions. (If `gvcfCalling` parameter is used)
   


