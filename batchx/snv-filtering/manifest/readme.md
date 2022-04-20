Performs SNV and INDEL filtering using the [GATK](https://gatk.broadinstitute.org/hc/en-us) suite.

# Details

SNVs and INDELS are filtered separately using GATK recommended parameters for each type of variant. For SNVs the filter are:
```sh
-filter "QD < 2.0" --filter-name "QD2" \
-filter "QUAL < 30.0" --filter-name "QUAL30" \
-filter "SOR > 3.0" --filter-name "SOR3" \
-filter "FS > 60.0" --filter-name "FS60" \
-filter "MQ < 40.0" --filter-name "MQ40" \
-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
```

For INDELs the filters are:
```sh
-filter "QD < 2.0" --filter-name "QD2" \
-filter "QUAL < 30.0" --filter-name "QUAL30" \
-filter "FS > 200.0" --filter-name "FS200" \
-filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
```

Then, both types of variants are merged.

# Inputs
## Required inputs
This tool has the following **required** inputs:
1. `refFolder`
Reference folder containing the reference genome in FASTA format, a dictionary file ending in `.dict` and an index file ending in `.fai`. The dictionary can be created using the command `gatk-launch CreateSequenceDictionary -R ref.fasta` where `ref.fasta` is the reference sequence. The index file can be created using the command `samtools faidx ref.fasta`. These two files are used by gatk HaplotypeCaller (and other tools) to efficiently access records of the genome. The reference must be the same as the one used to create de BAM file.
2. `vcfFile`
VCF (Variant Call Format) file containing all called SNVs and INDELs.
3. `sample`
Sample ID.

# Outputs
1. `filteredVcf`
VCF file containing filtered SNVs and INDELs.
   
# Links
- [gatk HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller)
- [FASTA - Reference genome format](https://gatk.broadinstitute.org/hc/en-us/articles/360035531652-FASTA-Reference-genome-format)
- [VCF - Variant Call Format](https://gatk.broadinstitute.org/hc/en-us/articles/360035531692-VCF-Variant-Call-Format)
- [gatk VariantFiltration](https://gatk.broadinstitute.org/hc/en-us/articles/360037434691-VariantFiltration)
- [filter_vep](https://m.ensembl.org/info/docs/tools/vep/script/vep_filter.html)

# Tool version
- gatk 4.2.5.0



