# Variant-calling-using-GATK4
Variant calling using GATK4 for a single sample

## Objectives
- Take raw DNA sequencing reads and perform variant calling to produce a variant list using GATK4.
- Perform basic exploration of variants.

## Description

This analysis runs through the GATK4 best practices workflow for variant calling. The workflow starts with pairs of sequencing reads and performs a series of steps to determine a set of genetic variants.

**Data**: Illumina HiSeq paired-end (2×100 bp) **chr20** reads in FASTQ format.
**Source**: Sample NA12878 from [IGSR](https://www.internationalgenome.org/data-portal/sample/NA12878)
**Reference**: [Grch38/Hg38](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0/) GATK bundle of [reference](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811)

## Tools and installation

The following tools are used in this analysis. Most of these tools are installed using [Homebrew](https://brew.sh/):

1. [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
``` brew install fastqc```
2. [multiqc](https://multiqc.info/)
```pip install multiqc```
3. [samtools 1.19](https://www.htslib.org/)
```brew install samtools```
4. [bwa 0.7.17](https://github.com/lh3/bwa)
```brew install bwa```
5. [picard 3.1.1](https://broadinstitute.github.io/picard/)
Install [picard.jar](https://github.com/broadinstitute/picard/releases/tag/3.1.1) and most recent Java version using brew
```brew install openjdk@17```
6. [GATK 4.1.3.0](https://gatk.broadinstitute.org/hc/en-us)
This analysis is performed using GATK docker interactively:
```docker pull broadinstitute/gatk:4.1.3.0```



Data preprocessing includes read trimming, alignment, sorting by coordinate, and marking duplicates.
