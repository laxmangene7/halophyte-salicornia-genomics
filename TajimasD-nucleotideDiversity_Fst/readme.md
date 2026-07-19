# Tajima's D, Nucleotide Diversity (π), and Fst Analyses

This directory contains the scripts and input files used to calculate and visualize genome-wide nucleotide diversity (π), Tajima's D, and genetic differentiation (Weir and Cockerham's FST) for the Salicornia population genomic analyses.

## Input Files

### Genome-wide Tajima's D and π (_Salicornia_ sp. France reference) for all chromosomes

- `S_sp.France.sites.with.S.France.Ref_50kb.window_All.chr.Tajima.D`
- `S_sp.France.sites.with.S.France.Ref_50kb.window_All.chr.pi`
- `S_procumbens.sites.with.S.France.Ref_50kb.window_All.chr.Tajima.D`
- `S_procumbens.sites.with.S.France.Ref_50kb.window_All.chr.pi`

These files were generated from variants mapped to the **Salicornia_sp. France** reference genome and were used to produce genome-wide Tajima's D and nucleotide diversity (π) plots (Supplementary Figs. 17 and 18).

### Chr6A Tajima's D and π (_S_. _ramosissima_ reference)

- `S_sp.France.sites.with.S_ramosissima.ref.50kb.window_Chr6A_tajimaD.txt`
- `S_sp.France.sites.with.S_ramosissima.ref.50kb.window_Chr6A_pi`
- `S_procumbens.sites.with.S_ramosissima.ref.50kb.window_Chr6A_tajima_D.txt`
- `S_procumbens.sites.with.S_ramosissima.ref.50kb.window_Chr6A_pi`

These files were generated from variants remapped to the **_Salicornia ramosissima_** reference genome and were used for the Chr6A regional analyses (Extended Data Fig. 10a,b).

### Sliding-window Fst

- `Fst_10kb.window.5kb.step.S_procumbens.vs.S.sp.France.weir.pop.fst`

The file contains Weir and Cockerham's FST estimates calculated using 10-kb windows with a 5-kb step size across Chr6A. The sites were called on S. ramosissima

## Scripts

- `S_France.and.S.procumbens_All.Chr.TajimasD_and.Pi.R`
  - Generates genome-wide Tajima's D and nucleotide diversity (π) plots for chromosomes Chr1A–Chr9A and Chr1C–Chr9C (Supplementary Figs. 17 and 18).

- `extended.data.fig.10a-b.TajimasD.Pi.R`
  - Generates Chr6A regional Tajima's D and nucleotide diversity (π) plots using variants mapped to the **S. ramosissima** reference genome (Extended Data Fig. 10a,b).

- `Fst.R`
  - Generates the sliding-window Fst plot for Chr6A (Supplementary Fig. 20).

### Sample-size validation analysis

- `S_sp.France.subset1.random.sample.called.on.S.ramosissima..50kb.window.Tajima.D.txt`
- `tajimasD.S.France.subset1_S.procumbens.R`

Generates Tajima's D plots for Chr6A using a random subset of 12 *Salicornia* sp. France accessions and 12 *S. procumbens* accessions. Variants were called against the *Salicornia ramosissima* reference genome. This analysis was performed to verify that the observed patterns were robust to differences in sample size between populations.

## Notes

- Genome-wide π and Tajima's D were calculated in non-overlapping 50-kb windows.
- Chr6A regional π and Tajima's D were also calculated using non-overlapping 50-kb windows.
- Fst was calculated using 10-kb sliding windows with a 5-kb step size.
- Variants used for the Chr6A regional analyses were called against the **Salicornia ramosissima** reference genome, whereas genome-wide analyses were performed using variants mapped to the **Salicornia sp. France** reference genome.
