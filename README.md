# Code Repository for ["Human genetics implicates a BACH2-NRF2 axis in fetal hemoglobin activation"](https://www.nature.com/articles/s41586-026-11113-2)
 
Citation: 
> Guo, CJ., Arora, U.P., Cheng, X. et al. Human genetics implicates a BACH2–NRF2 axis in fetal haemoglobin activation. *Nature* (2026). https://doi.org/10.1038/s41586-026-11113-2 

## Overview

This repository hosts scripts used for the identification of a BACH2–NRF2 regulatory axis regulating fetal hemoglobin (HbF) expression. The analyses are grouped into the following categories:

1. **[`StatGen_analysis/`](StatGen_analysis/)** <br>
   Includes GWAS meta-analyses, conditional analysis, fine-mapping, heritability & genetic-correlation analyses, and multi-ancestry fine-mapping of the BACH2 locus.

2. **[`SCAVENGE_analysis/`](SCAVENGE_analysis/)**<br>
   Propagates fine-mapped GWAS variant scores across cell graph constructed based on Granja *et al.* (2019)'s bone marrow scATAC-seq data to identify the hematopoietic cell types most relevant to HbF regulation.

3. **[`BulkRNAseq_analysis/`](BulkRNAseq_analysis/)**<br>
   Characterizes transcriptional changes upon BACH2 knockdown, highlighting upregulation of fetal globin genes (*HBG1/HBG2*). 
   Pipeline includes STAR alignment, featureCounts quantification, and DESeq2 differential expression for BACH2-sh2 vs. scramble control.

4. **[`CUTRUN_analysis/`](CUTRUN_analysis/)**<br>
   Maps BACH2 and NRF2 occupancy genome-wide in erythroid cells, linking GWAS variants to regulatory elements. Pipeline includes End-to-end CUT&RUN pipeline: trimming, spike-in normalization, hg38 alignment, peak calling, bigWig tracks, and FIMO motif scanning.
