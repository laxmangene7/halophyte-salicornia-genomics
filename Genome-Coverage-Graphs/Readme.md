# Genome Coverage Graphs

This directory contains the scripts used to generate normalized genome coverage plots from Illumina whole-genome resequencing data. These plots were used to evaluate chromosome-wide read-depth patterns, identify large structural differences, and infer ploidy and chromosome copy-number variation.

---

# Workflow

## 01. Build HISAT2 Reference Index

Build a HISAT2 index from the reference genome.

**Script**

- `01_hisat2_build_index.sh`

**Input**

- `Salicornia_France.S.bigelovii.combined.ref.noChrun.ref.fa`

**Output**

- HISAT2 index files (`Salicornia_France.S.bigelovii.combined.ref.noChrun.ref.*.ht2`)

---

## 02. Read Quality Control

Trim adapters and low-quality bases using **fastp**.

**Script**

- `02_fastp_trim.sh`

**Input**

- Paired-end FASTQ files

**Output**

- Trimmed paired-end FASTQ files
- HTML and JSON quality reports

---

## 03. Genome Alignment

Align paired-end reads to the combined **_Salicornia sp. France–S. bigelovii_** reference genome using HISAT2.

**Script**

- `03_hisat2_alignment.sh`

**Input**

- Trimmed paired-end FASTQ files
- HISAT2 reference index

**Output**

- SAM alignment files

---

## 04. Generate 100-kb Genome Coverage

Generate genome-wide read counts in non-overlapping 100-kb bins using uniquely mapped concordant paired-end alignments.

**Script**

- `04_generate_100kb_read_counts.sh`

**Input**

- SAM alignment files

**Pipeline**

- Remove SAM header lines.
- Retain concordant paired-end reads (`YT:Z:CP`).
- Retain uniquely mapped reads (`NH:i:1`).
- Extract chromosome names and genomic positions.
- Assign reads to non-overlapping 100-kb genomic bins.
- Count the number of reads within each genomic bin.

**Output**

- `<sample>_100kb.txt`

Output format:

- Column 1: Read count
- Column 2: Chromosome
- Column 3: 100-kb bin start position

Each output file contains the number of uniquely mapped concordant read pairs within each non-overlapping 100-kb genomic interval.

---

## 05. Normalize Read Counts

Normalize read counts to account for differences in sequencing depth among samples.

**Normalization procedure**

1. Calculate the average read count of the top 10% highest-coverage bins.
2. Define a threshold equal to 10% of this average.
3. Calculate the average coverage using bins above this threshold.
4. Normalize each genomic bin by this average coverage.

**Script**

- `05_normalize_read_counts.sh`

**Input**

- `*_100kb.txt`

**Output**

- Normalized genome coverage files (`*.norm.txt`)

---

## 06. Plot Genome Coverage

Generate chromosome-wide normalized genome coverage plots.

**Scripts**

- `R_script.plotting.R`
- `Run.R.sh`

**Input**

- `*.norm.txt`

**Output**

- `*_Normalized.reads.pdf`

Each figure displays normalized read depth across all chromosomes. Normalized read depth is plotted across successive 100-kb genomic bins, enabling visualization of chromosome copy-number variation, large structural differences, and ploidy.

---

# Software

- HISAT2
- fastp
- samtools
- GNU awk
- R
- data.table
- ggplot2

---

# Notes

- Adjust the SLURM `--array` parameter according to the number of samples in your dataset.
- Update input directories and reference genome paths to match your computing environment.
- Genome coverage was calculated using only uniquely mapped concordant paired-end reads (`NH:i:1` and `YT:Z:CP`), minimizing mapping artifacts in repetitive genomic regions.
- Read-depth normalization corrects for variation in sequencing depth among samples, enabling direct comparison of chromosome-wide coverage profiles.
