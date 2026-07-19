# Variant Calling

This directory contains the scripts used for preprocessing sequencing reads, mapping reads to the *Salicornia ramosissima* reference genome, and calling SNPs from the global *Salicornia* resequencing panel.

---

## Workflow

### Step 1. Quality Control and Read Trimming

Run:

```bash
trim_raw.fq.sh
```

This script uses **fastp** to:

- trim low-quality bases
- detect and remove paired-end adapters
- discard reads shorter than 150 bp
- remove duplicate reads (`--dedup`)
- generate HTML and JSON quality reports

**Note:** Adjust the SLURM array size to match the total number of samples.

---

### Step 2. Read Alignment

Run:

```bash
hisat2.sh
```

Reads are aligned to the **Salicornia ramosissima** reference genome using **HISAT2**.

The alignment was performed using the following key parameters:

- paired-end alignment
- `--no-spliced-alignment` (genomic DNA alignment)
- `--no-unal` (do not report unmapped reads)
- maximum fragment length `-X 1000`

Output:

- SAM alignment files
- alignment log files

**Note:** Adjust the SLURM array size according to the number of samples.

---

### Step 3. Extract Uniquely Mapped Concordant Reads

Run:

```bash
unique-concor.sh
```

This script retains only uniquely mapped concordant read pairs using the HISAT2 alignment tags:

- `YT:Z:CP`
- `NH:i:1`

The filtered alignments are written as BAM files.

**Note:** Update the SLURM array size according to the number of SAM files.

---

### Step 4. Sort and Index BAM Files

Run:

```bash
bam.to.sorted.bam1.sh
```

This script:

- sorts BAM files using **samtools sort**
- creates CSI index files using **samtools index**

Input BAM files are listed in:

```text
all.bamfile.list2.txt
```

**Note:** Adjust the SLURM array size to match the number of BAM files.

---

### Step 5. SNP Calling

Run:

```bash
snp_Salicornia.sh
```

Variants are called using **bcftools mpileup** followed by **bcftools call**.

The workflow:

- computes genotype likelihoods
- annotates allele depth (AD) and read depth (DP)
- excludes indels
- performs multi-sample SNP calling

Modify the chromosome list (`-r`) as needed to perform chromosome-wise variant calling.

---

## Software

- fastp
- HISAT2
- samtools
- bcftools

---

## Input Files

- Raw paired-end FASTQ files
- HISAT2 reference index for *Salicornia ramosissima*
- `bamFile.list.txt`
- `all.bamfile.list2.txt`

---

## Output

The pipeline produces:

- Quality-controlled FASTQ files
- HISAT2 alignment (SAM) files
- Filtered BAM files containing uniquely mapped concordant reads
- Sorted and indexed BAM files
- Chromosome-specific VCF files for downstream population genomic analyses
