# Code Repository for ["Human genetics implicates a BACH2-NRF2 axis in fetal hemoglobin activation"](https://www.nature.com/articles/s41586-026-11113-2)
 
Citation: 
> Guo, CJ., Arora, U.P., Cheng, X. et al. Human genetics implicates a BACH2–NRF2 axis in fetal haemoglobin activation. *Nature* (2026). https://doi.org/10.1038/s41586-026-11113-2 

## Overview

This repository hosts scripts used for the identification of a BACH2–NRF2 regulatory axis regulating fetal hemoglobin (HbF) expression. The analyses are grouped into the following categories:

1. **Statistical genetics** — Meta-GWAS, fine-mapping, and heritability/genetic-correlation analyses identify and characterize HbF-associated loci, implicating BACH2.
2. **SCAVENGE** — Propagates fine-mapped GWAS variant scores across a bone marrow scATAC-seq cell graph to identify the hematopoietic cell types most relevant to HbF regulation.
3. **Bulk RNA-seq** — Characterizes transcriptional changes upon BACH2 knockdown, highlighting upregulation of fetal globin genes (*HBG1/HBG2*).
4. **CUT&RUN** — Maps BACH2 and NRF2 occupancy genome-wide in erythroid cells, linking GWAS variants to regulatory elements.


## Repository structure

| Directory | Description |
|-----------|-------------|
| [`StatGen_analysis/`](StatGen_analysis/) | Meta-GWAS, COJO fine-mapping, SNP heritability estimation, and genetic correlation analyses |
| [`SCAVENGE_analysis/`](SCAVENGE_analysis/) | Per-cell trait relevance scoring using SCAVENGE on Granja 2019 bone marrow scATAC-seq data |
| [`BulkRNAseq_analysis/`](BulkRNAseq_analysis/) | STAR alignment, featureCounts quantification, and DESeq2 differential expression for BACH2-sh2 vs. scramble control |
| [`CUTRUN_analysis/`](CUTRUN_analysis/) | End-to-end CUT&RUN pipeline: trimming, spike-in normalization, hg38 alignment, peak calling, bigWig tracks, and FIMO motif scanning |
