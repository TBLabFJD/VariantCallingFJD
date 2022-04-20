Calculates BAM coverage quality check for SNV/INDEL and CNV calling using [`mosdepth`](https://github.com/brentp/mosdepth).


# Context

Next-generation sequencing (NGS) coverage describes the average number of reads that align to, or "cover," known reference bases. The sequencing coverage level often determines whether variant discovery can be made with a certain degree of confidence at particular base positions.

# Purpose

This tool performs a sample quality check for SNV/INDEL and CNV calling. It is inspired by the script https://github.com/TBLabFJD/VariantCallingFJD/blob/master/tasks/mosdepth.sh used to calculate coverage on different scenarios generating as output a folder containing the report.

When a BED file is given, the proportion of covered bases with a given read depth is calculated. 

## SNV/Indels
For SNVs/INDELs, if less than 90% of the bases in the target regions have a read depth above 10 reads, the sample is excluded from the analysis. 

## CNVs

For CNVs, as the implemented algorithms are read depth calling based, by default the sample must have a coverage of 40 reads in at least 90% of the targeted bases. The report also includes [quantized](https://github.com/brentp/mosdepth#quantize) metrics.

# Inputs
## Required inputs
This tool has the following **required** inputs:
1. `bamFile`
File containing the mapped reads with the extension `.bam`. A BAM file (*.bam) is the compressed binary version of a SAM file that is used to represent aligned sequences. 
2. `baiFile`
Index file of the `bamFile` with the extension `.bai` or `.bam.bai`. It stores the genomic positions to efficiently access records of the BAM file. 
3. `sample` 
Sample ID.
4. `analysisType`
Type of analysis: "snv" o "cnv".

## Optional inputs
1. `panelFile`
BED file containing the captured regions of the genome based on the library design. It is a TSV file without a header. It has a minimum of 3 columns storing the chromosome, start and stop position for each target. For CNV calling, it needs a 4th column containing the gene name of the target. In this step, the `panelFile` is used to calculate the coverage distribution among the theoretical covered regions. If not provided, the coverage is calculated for all the genome.
2. `intervalPadding`
Number of bases at each side of the theoretical captured regions to extend the analysis to those regions also. In capture strategies for library preparation, it is common to capture the adjoining regions to the targeted ones with high enough coverage to perform variant calling.
3. `minCoverage`
Minimum read depth (number of reads) that at least 90% of the theoretical captured region must have to consider a sample as a good quality one. 10 by default.

# Outputs
1. `report`
Folder containing the following files:
   - `minCovFilterResults_snv.txt` (if analysisType is set to “snvs”): TSV file containing the following columns: sampeID, “10”, proportion of targeted regions covered with a sequencing depth of more than 10 reads, PASSED/FAILED tag if the previous column is higher/lower respectively than 0.9 (90%), run name, `bamFile` path.
   - `minCovFilterResults_cnv.txt` (if analysisType is set to “cnvs”): TSV file containing the following columns: sampeID, `readThrehold`, proportion of targeted regions covered with a sequencing depth of more than `readThrehold` reads, PASSED/FAILED tag if the previous column is higher/lower respectively than 0.9 (90%), run name, `bamFile` path.


# Exit codes
If the coverage requirements are met (at least a 90% of the covered regions must have a coverage greater or equal `minCoverage`), the jobs finalize in `SUCCEEDED` state (exit code 0), otherwise they finalize in `WARNING` state (exit code 1).

# Links

- **mosdepth**: https://github.com/brentp/mosdepth
- **bedtools**: https://github.com/arq5x/bedtools2

# Tools version
- **Bedtools**: 2.30.0
- **mosdepth**: 0.3.3
