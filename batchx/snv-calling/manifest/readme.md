# Context
SNV and INDEL calling uses [gatk HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/4418062719899-HaplotypeCaller) in single-sample mode or in `GVCF` mode. A bed file can be provided to limit the variant calls inside the targeted regions plus padding (extension of x bases to each side of a target region). 
When run in single sample mode it generates a Variant Calling Format (VCF) file and when run in GVCF mode it generates a Genomic VCF (GVCF) file. 

While VCF file contains just variant records, GVCF contains records for all sites (variant and non variant sites). Single-sample VCF is used for analyzing one sample. GVCF can be used to create a multi-sample vcf to analyze multiple samples from a cohort. The extra regions from the GVCF are used during the merging step to differentiate between covered-non-variant sites and non-covered sites.

# Inputs
## Required inputs
This tool has the following **required** inputs:
1. `refFolder`
Reference folder containing the reference genome in FASTA format, a dictionary file ending in `.dict` and an index file ending in `.fai`. The dictionary can be created using the command `gatk-launch CreateSequenceDictionary -R ref.fasta` where `ref.fasta` is the reference sequence. The index file can be created using the command `samtools faidx ref.fasta`. These two files are used by gatk HaplotypeCaller (and other tools) to efficiently access records of the genome. The reference must be the same as the one used to create de BAM file.
2. `bamFile`
File containing the mapped reads with the extension `.bam`. A BAM file (*.bam) is the compressed binary version of a SAM file that is used to represent aligned sequences. 
3. `bamIndex`
Index file of the bamfile with the extension `.bai` or `.bam.bai`. It stores the genomic positions to efficiently access records of the BAM file. The index file can be created using the command `samtools faidx sample.bam` where `sample.bam` is the BAM file from which we want to ceate the index.
4. `outputPrefix`
Sample ID.

## Optional inputs
1. `gvcfCalling`
Argument used to run gatk HaplotypeCaller in GVCF mode.
2. `panelFile`
BED file containing the captured regions of the genome based on the library design. It is a TSV file without a header. It has a minimum of 3 columns storing the chromosome, start and stop position for each target. For CNV calling, it needs a 4th column containing the gene name of the target. In this step, the `panelFile` is used to calculate the coverage distribution among the theoretical covered regions. If not provided, the coverage is calculated for all the genome.

# Outputs
1. `outputVcf`
- VCF (Variant Call Format) file containing all called SNVs and INDELs. (default)
- GVCF (Genomic Variant Call Format) file containing  all called SNVs and INDELs and all covered regions. (If `gvcfCalling` parameter is used)
   
# Links
- [gatk HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller)
- [FASTA - Reference genome format](https://gatk.broadinstitute.org/hc/en-us/articles/360035531652-FASTA-Reference-genome-format)
- [VCF - Variant Call Format](https://gatk.broadinstitute.org/hc/en-us/articles/360035531692-VCF-Variant-Call-Format)
- [GVCF - Genomic Variant Call Format](https://gatk.broadinstitute.org/hc/en-us/articles/360035531812-GVCF-Genomic-Variant-Call-Format)

# Tool version
- gatk 4.2.5.0


