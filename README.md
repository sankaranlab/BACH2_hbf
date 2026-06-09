# BACH2_hbf

Scripts and supplementary files for the manuscript:

> **"Human genetics implicates a BACH2-NRF2 axis in fetal hemoglobin activation"**  
> Preprint: https://www.medrxiv.org/content/10.1101/2023.03.24.23287659v3

## Table of contents

- [Overview](#overview)
- [Repository structure](#repository-structure)
- [`StatGen_analysis/`](StatGen_analysis/README.md) — Statistical genetics (meta-GWAS, fine-mapping, heritability, genetic correlation)
- [`SCAVENGE_analysis/`](SCAVENGE_analysis/README.md) — Trait-relevant cell-type enrichment from scATAC-seq
- [`BulkRNAseq_analysis/`](BulkRNAseq_analysis/README.md) — Bulk RNA-seq processing and differential expression
- [`CUTRUN_analysis/`](CUTRUN_analysis/README.md) — CUT&RUN chromatin profiling (BACH2 and NRF2 binding)

---

## Overview

The analyses in this repository support the identification of a BACH2–NRF2 regulatory axis controlling fetal hemoglobin (HbF) expression. The workflow broadly proceeds as follows:

1. **Statistical genetics** — Meta-GWAS, fine-mapping, and heritability/genetic-correlation analyses identify and characterize HbF-associated loci, implicating BACH2.
2. **CUT&RUN** — Maps BACH2 and NRF2 occupancy genome-wide in erythroid cells, linking GWAS variants to regulatory elements.
3. **Bulk RNA-seq** — Characterizes transcriptional changes upon BACH2 knockdown, highlighting upregulation of fetal globin genes (*HBG1/HBG2*).
4. **SCAVENGE** — Propagates fine-mapped GWAS variant scores across a bone marrow scATAC-seq cell graph to identify the hematopoietic cell types most relevant to HbF regulation.

---

## Repository structure

| Directory | Description |
|-----------|-------------|
| [`StatGen_analysis/`](StatGen_analysis/) | Meta-GWAS, COJO fine-mapping, SNP heritability estimation, and genetic correlation analyses |
| [`SCAVENGE_analysis/`](SCAVENGE_analysis/) | Per-cell trait relevance scoring using SCAVENGE on Granja 2019 bone marrow scATAC-seq data |
| [`BulkRNAseq_analysis/`](BulkRNAseq_analysis/) | STAR alignment, featureCounts quantification, and DESeq2 differential expression for BACH2-sh2 vs. scramble control |
| [`CUTRUN_analysis/`](CUTRUN_analysis/) | End-to-end CUT&RUN pipeline: trimming, spike-in normalization, hg38 alignment, peak calling, bigWig tracks, and FIMO motif scanning |
