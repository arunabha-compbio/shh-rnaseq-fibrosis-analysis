# shh-rnaseq-fibrosis-analysis
# SHH RNA-seq Pathway Analysis in Fibrosis & Repair

RNA-seq pathway-focused transcriptomic analysis of the Sonic Hedgehog (SHH) signaling pathway using mouse featureCounts data (Control vs GR1).

---

## Project Overview

This project investigates dysregulation of the Sonic Hedgehog (SHH) signaling pathway during fibrosis-associated tissue remodeling and repair using bulk RNA-seq count data.

The analysis was performed in RStudio using a featureCounts-generated count matrix derived from mouse RNA-seq samples.

Key objectives:

- Extract SHH pathway genes
- Analyze pathway activation/suppression
- Connect SHH signaling with fibrosis and repair
- Perform GO and KEGG enrichment analysis
- Visualize transcriptomic changes using multiple plot types

---

## Dataset

Input data:

- Mouse RNA-seq featureCounts output
- Samples:
  - Control
  - GR1

RNA-seq workflow:

Biological sample  
↓  
RNA extraction  
↓  
cDNA library preparation  
↓  
NGS sequencing  
↓  
FASTQ files  
↓  
Genome alignment  
↓  
BAM files  
↓  
featureCounts  
↓  
counts.txt

---

## Tools & Packages

### Software

- RStudio
- featureCounts

### R Packages

- ggplot2
- pheatmap
- ggrepel
- clusterProfiler
- enrichplot
- org.Mm.eg.db

---

## SHH Pathway Genes Analyzed

### Ligands
- Shh
- Ihh
- Dhh

### Receptors / Co-receptors
- Ptch1
- Ptch2
- Smo
- Gas1
- Cdon
- Boc

### Regulators
- Sufu
- Kif7
- Hhip

### GLI Transcription Factors
- Gli1
- Gli2
- Gli3

### Downstream Targets
- Mycn
- Snai1
- Ccnd1

---

## Key Findings

- Canonical SHH signaling markers (Gli1, Ptch1, Shh) were suppressed in GR1.
- Gli2, Mycn and Snai1 were upregulated, suggesting non-canonical SHH-associated remodeling.
- Fibrosis-associated genes including Tgfb1, Col1a1, Col1a2, Fn1 and Acta2 were elevated.
- GO and KEGG enrichment analyses indicated inflammatory and extracellular matrix remodeling processes.

Overall interpretation:

> GR1 samples exhibit inflammatory fibrotic remodeling associated with altered SHH pathway regulation.

---

## Generated Visualizations

- SHH pathway heatmap
- Fold-change bar plots
- Functional category plots
- Bubble plots
- Fibrosis/repair pathway plots
- GO enrichment analysis
- KEGG pathway enrichment analysis

---

## Repository Contents

| File | Description |
|------|-------------|
| `Shh_fibrosis.R` | Full R analysis pipeline |
| `counts.txt` | featureCounts RNA-seq count matrix |
| `SHH_analysis_report.pdf` | Internship assignment report |
| `plots/` | Generated figures |
| `results/` | Exported processed results |

---

## Limitations

- Only one control and one GR1 sample were available.
- Results are exploratory and not statistically conclusive due to lack of biological replicates.

---

## Author

Arunabha Pal  
St. Xavier’s College (Autonomous), Kolkata

Internship Project — Epigenomics Lab, University of Calcutta
