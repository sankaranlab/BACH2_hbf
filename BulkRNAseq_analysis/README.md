# Bulk RNA-seq Analysis

Scripts for processing bulk RNA-seq data from a BACH2 knockdown experiment in erythroid cells (day-11 differentiation, scramble control vs. BACH2-sh2) and identifying differentially expressed genes.

---

## [`1_rawData_process.sh`](1_rawData_process.sh) — Alignment and read counting

```bash
$ bash 1_rawData_process.sh <sample_name>
```

Aligns paired-end reads to hg38 (STAR), filters, and quantifies per-gene counts (featureCounts).

| | |
|---|---|
| **Input** | `raw_data/<sample>_R1/R2_001.fastq.gz` |
| **Output** | `output/<sample>/<sample>.count.bed` |
| **Tools** | STAR, samtools, featureCounts (Subread) |

---

## [`2_DiffGene.identification.r`](2_DiffGene.identification.r) — Differential expression

```bash
$ Rscript 2_DiffGene.identification.r
```

Runs DESeq2 on scramble vs. BACH2-sh2 samples (2 replicates each) and generates volcano plots. DEGs called at p < 0.05 and |log2FC| > 0.5.

| | |
|---|---|
| **Input** | `<sample>.count.bed` files from step 1 |
| **Output** | `DEGs.p=<threshold>.2fold.scatter.png`, `genelist.png` |
| **Tools** | R, DESeq2, ggplot2 |
