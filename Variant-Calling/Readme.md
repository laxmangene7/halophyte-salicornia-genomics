# Variant Calling

This directory contains the scripts used for whole-genome variant discovery from paired-end Illumina resequencing data. The workflow includes reference genome indexing, read trimming, genome alignment, extraction of uniquely mapped reads, BAM sorting and indexing, SNP calling, and variant filtering.

---

## Workflow

### 01. Build HISAT2 Reference Index

Build a HISAT2 index from the reference genome before alignment.

**Script**

- `01_hisat2_build_index.sh`

**Input**

- `Salicornia_ramosissima.ref.fa`

**Output**

- HISAT2 index files (`Salicornia_ramosissima.ref.*.ht2`)

---

### 02. Read Quality Control

Trim adapters and low-quality bases using **fastp**.

**Script**

- `02_fastp_trim.sh`

**Input**

- Paired-end FASTQ files (`*_1.fq.gz`, `*_2.fq.gz`)

**Output**

- Trimmed FASTQ files
- HTML quality reports
- JSON reports

---

### 03. Read Alignment

Align trimmed paired-end reads to the **Salicornia ramosissima** reference genome using HISAT2.

**Script**

- `03_hisat2_alignment.sh`

**Input**

- Trimmed FASTQ files
- HISAT2 reference index

**Output**

- SAM alignment files

---

### 04. Extract Uniquely Mapped Concordant Reads

Retain only uniquely mapped concordant read pairs for downstream variant calling.

**Script**

- `04_extract_unique_concordant_reads.sh`

**Input**

- SAM files

**Output**

- BAM files containing uniquely mapped concordant reads

---

### 05. Sort and Index BAM Files

Sort BAM files and create BAM index files.

**Script**

- `05_sort_and_index_bam.sh`

**Input**

- BAM files

**Output**

- Sorted BAM files (`*.s.bam`)
- BAM index files (`*.s.bam.csi`)

---

### 06. Variant Calling

Call SNPs using **bcftools mpileup** and **bcftools call**.

**Script**

- `06_variant_calling_bcftools.sh`

**Input**

- Sorted BAM files
- Reference genome

**Output**

- Raw VCF file

---

### 07. Basic Variant Filtering

Filter variants based on sequencing quality and allele frequency.

Filtering criteria:

- QUAL ≥ 30
- INFO/DP ≥ 20
- 0.01 < AF < 0.99

**Script**

- `07_basic_variant_filtering.sh`

**Input**

- Raw VCF

**Output**

- `chr3_4.filtered.vcf.gz`

---

### 08. Population-Level Variant Filtering

Apply additional filtering based on population genetic criteria.

Filtering criteria:

- Minor allele frequency (MAF) ≥ 0.02
- Missing genotypes ≤ 20%
- Heterozygosity ≤ 20%

**Scripts**

- `08_population_variant_filtering.sh`
- `vcf_MAF_Missing_Het_filter.pl`

**Input**

- `chr3_4.filtered.vcf`

**Output**

- Population-filtered VCF used for downstream analyses.

---

## Software

- HISAT2
- fastp
- samtools
- bcftools
- Perl
- tabix

---

## Notes

- Adjust the `--array` parameter in each SLURM script to match the number of samples in your dataset.
- Update file paths to match your computing environment.
- The final filtered VCF generated in Step 08 was used for downstream population genomic analyses, including phylogenetic reconstruction, population structure, nucleotide diversity (π), Tajima's D, genetic differentiation (FST), haplotype analysis, and read-depth analyses.
