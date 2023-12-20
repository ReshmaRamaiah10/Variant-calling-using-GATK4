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
```
brew install fastqc
 ```
2. [multiqc](https://multiqc.info/)
```
pip install multiqc
```
3. [samtools 1.19](https://www.htslib.org/)
```
brew install samtools
```
4. [bwa 0.7.17](https://github.com/lh3/bwa)
```
brew install bwa
```
5. [picard 3.1.1](https://broadinstitute.github.io/picard/)

Install [picard.jar](https://github.com/broadinstitute/picard/releases/tag/3.1.1) and most recent Java version using brew
```
brew install openjdk@17
```
6. [GATK 4.1.3.0](https://gatk.broadinstitute.org/hc/en-us)

This analysis is performed using [GATK docker](https://gatk.broadinstitute.org/hc/en-us/articles/360035889991) interactively:
```
docker pull broadinstitute/gatk:4.1.3.0
```

## Data preprocessing

Data preprocessing includes read trimming, alignment, sorting by coordinate, and marking duplicates.

1.  Quality control checks on the raw fastq data

FastQC (Fast Quality Control) is designed to assess the quality of high-throughput sequencing data. Its main purpose is to provide a comprehensive analysis of various quality-related metrics and characteristics of raw sequencing data. Here are the primary purposes and features of FastQC:
- Per-Base Sequence Quality
- Per-Sequence Quality Scores
- Sequence Length Distribution
- GC Content Distribution
- Sequence Duplication Levels
- Overrepresented Sequences
- Adapter Content
- Kmer Content
- Quality per Base Position

Run fastqc on both read 1 and read 2 fastq files:
```
fastqc NA12878.chr20.region_1.fastq.gz -o fatsqc_results/
```

MultiQC is a tool used to aggregate, analyze, and visualize quality control (QC) results from multiple samples. It creates a unified report and interactive HTML reports to view QC metrics of both `NA12878.chr20.region_1.fastq.gz` and `NA12878.chr20.region_2.fastq.gz` files.
```
multiqc fatsqc_results/
```
**Results**: `multiqc_report.html` report depicts good sequence/read quality with good quality scores, GC content and sequence length distribution.

2. Map raw mapped reads to reference genome

Map the raw sequencing data to the Homo sapiens genome assembly GRCh38 (hg38) using the BWA-MEM algorithms. Use Samtools to convert the output to the BAM format.
```
bwa mem -M -t 2 -R "@RG\tID:SRR622461.7\tSM:NA12878\tLB:ERR194147\tPL:ILLUMINA" \
reference/Homo_sapiens_assembly38.fasta \
NA12878.chr20.region_1.fastq.gz \
NA12878.chr20.region_2.fastq.gz | \
samtools view -b -h -o NA12878.bam -
```
At the end of this step you should have a file called `NA12878.bam`.

3. Sort SAM/BAM

The alignment file NA12878.bam is sorted using the Picard tools.
```
java -jar -Xmx7g picard.jar SortSam -I NA12878.bam -O NA12878.sort.bam
```
The above command will create a coordinate sorted BAM file and an index (.bai) file.

4. Mark duplicate reads

The aim of this step is to locate and tag duplicate reads in the BAM file.
```
java -jar -Xmx7g picard.jar MarkDuplicates -I NA12878.sort.bam -O NA12878.sort.dup.bam -M marked_dup_metrics.txt
```

## Base quality recalibration

The last step of pre-processing mapped reads is the base quality score recalibration (BQSR) stage. The GATK tools detects systematic errors made by the sequencing machine while estimating the accuracy of each base. The systematic errors can have various sources ranging from technical machine errors to the variability in the sequencing chemical reactions.

1. We need the GATK tool on its docker from here on. After pulling the docker image as mentioned in the installation step, start a GATK container interactively by mounting your directory with required files:
```
docker run -v /Users/ramaiahr/Documents/wgs:/gatk/my_data -it broadinstitute/gatk:4.1.3.0
```
2. The following two step BQSR process applies machine learning to model the possible errors and adjust the base quality scores accordingly.
```
gatk BaseRecalibrator -I my_data/NA12878.sort.dup.bam \
-R my_data/reference/Homo_sapiens_assembly38.fasta \
--known-sites my_data/reference/dbsnp_146.hg38.vcf.gz \
-O recal_data.table
```
```
gatk ApplyBQSR -I my_data/NA12878.sort.dup.bam \
-R my_data/reference/Homo_sapiens_assembly38.fasta \
--bqsr-recal-file my_data/recal_data.table \
-O my_data/NA12878.sort.dup.bqsr.bam
```
We now have a pre-processed BAM file (NA12878.sort.dup.bqsr.bam) ready for variant calling.

3. The command below uses Picard to generate QC metrics. Run multiQC to aggregate it with fastq data and produce an HTML report.
```
java -jar -Xmx7g picard.jar CollectMultipleMetrics \
-I my_data/NA12878.sort.dup.bqsr.bam \
-R my_data/reference/Homo_sapiens_assembly38.fasta \
-O my_data/NA12878.sort.dup.bqsr.CollectMultipleMetrics
```

## Variant Calling using HaplotypeCaller

HaplotypeCaller is the focal tool within GATK4 to simultaneously call germline SNVs and small Indels using local de novo assembly of haplotype regions.
Briefly, the HaplotypeCaller works by: 1. Identify active regions or regions with evidence of variations. 2. Re-asssemble the active regions. 3. Re-align active region reads to the assembled regions to identify allele.
```
gatk HaplotypeCaller -I my_data/NA12878.sort.dup.bqsr.bam \
-R my_data/reference/Homo_sapiens_assembly38.fasta \
-L chr20 \
-O my_data/NA12878.vcf.gz
```

## Filter and prepare analysis ready variants

1. Consider the following method to filter a single sample VCF file. Here we will go through the Convolutional Neural Net based protocol to annotate and filter the VCF file.

This is a two step process:

(i) CNNScoreVariants will annotate the variant with pre-computed single-sample derived model scores in the INFO field CNN_1D (the neural network performs convolutions over the reference sequence surrounding the variant and combines those features with a multilayer perceptron on the variant annotations).
```
gatk CNNScoreVariants -V my_data/NA12878.vcf.gz \
-R my_data/reference/Homo_sapiens_assembly38.fasta \
-O my_data/filtered_output/NA12878.cnns.vcf
```

2. FilterVariantTranches takes as input the percent sensitivities (0-100) to known sites to apply the filter. Variants with scores higher than for e.g. 99th percentile of variants in the resources pass through the filter and will have PASS in the filter. Others will have a filter values like ‘CNN_1D_INDEL_Tranche_99.40_100.00’ or ‘CNN_1D_SNP_Tranche_99.95_100.00’.
```
gatk FilterVariantTranches -V my_data/filtered_output/NA12878.cnns.vcf \
--resource my_data/reference/hapmap_3.3.hg38.vcf.gz \
--resource my_data/reference/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
--info-key CNN_1D --snp-tranche 99.95 --indel-tranche 99.4 \
-O my_data/filtered_output/NA12878.cnns.cnnfilter.vcf
```

2. Additional filtering
The VariantFiltration tools is designed for hard-filtering variant calls based on custom quality criteria such as sequencing depth, mapping quality etc. The two parameters are the filter-name and filter-expression. The parameter filter-name is the name of the filter to be used in the FILTER column if the expression in filter-expression is true. In the example below, if the sequencing depth at the variant site (VCF field DP) is less than 10, the FILTER field will be populated with the value ‘Low_depth10’. Users can add multiple filter expression/name combinations.
```
gatk VariantFiltration -V my_data/NA12878.vcf.gz \
-R my_data/reference/Homo_sapiens_assembly38.fasta \
-O my_data/output.vqsr.varfilter.vcf \
--filter-name "Low_depth10" \
--filter-expression "DP < 10"
```

## Exporting variant data and visualisation

Given we have a filter annotated VCF files (), we can now create an analysis ready VCF file.

```
gatk VariantsToTable -V my_data/output.vqsr.varfilter.vcf -R my_data/reference/Homo_sapiens_assembly38.fasta -F CHROM -F POS -F FILTER -F TYPE -GF AD -GF DP --show-filtered -O my_data/output.vqsr.varfilter.pass.tsv
```
