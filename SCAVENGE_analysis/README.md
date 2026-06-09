# SCAVENGE Analysis

Applies [SCAVENGE](https://github.com/sankaranlab/SCAVENGE) (Single Cell Analysis of Variant Enrichment through Network propagation of GEnomic data) to identify cell types most enriched for HbF GWAS signal in a human bone marrow scATAC-seq dataset.

## Overview

SCAVENGE propagates fine-mapped GWAS variant scores across a mutual k-nearest-neighbor (mKNN) graph of single cells (built from their ATAC-seq profiles), producing a per-cell trait relevance score (TRS). Cells in trait-relevant chromatin states receive high TRS values, enabling cell-type-level enrichment analysis.

## Input data

- **`input_data/hbf-META_finemap_20240712_p5e-8.txt`** — Fine-mapped HbF GWAS variants within genome-wide COJO-identified independent loci with posterior inclusion probabilities (PIPs).
- **`normal_bm.peaks.granja2019.rds`** — Normal bone marrow scATAC-seq peak × cell count matrix from [Granja et al. 2019](https://github.com/GreenleafLab/10x-scATAC-2019), aligned to hg19. This file can be downloaded via:
    ```bash
    # scATAC-seq Hematopoiesis cell x peak SummarizedExperiment (healthy)
    wget https://jeffgranja.s3.amazonaws.com/MPAL-10x/Supplementary_Data/Healthy-Data/scATAC-Healthy-Hematopoiesis-191120.rds \
         -O normal_bm.peaks.granja2019.rds
    ```

> **Note:** The fine-mapped variants were lifted over to hg19 to match the Granja 2019 scATAC-seq data, which was originally aligned to hg19.

## How to use

1. **Prepare inputs:**
   - Place your fine-mapped variants BED file in `input_data/` (e.g., `hbf_meta.finemap.hg19.bed`). The 5th column should contain posterior inclusion probabilities (PIPs).
   - Ensure the scATAC-seq peak × cell count matrix (as a `SummarizedExperiment` or compatible RDS file) is accessible.

2. **Edit script paths:**
   - Update `variant_file` to point to your BED file.
   - Update the `readRDS()` path to point to your scATAC-seq data file.

3. **Run the scripts in order:**
   ```bash
   # Liftover finemapped variants from hg38 to hg19
   $ python 1_liftover.py

   # Perform SCAVENGE analysis
   $ Rscript 2_scavenge.granja2019.R

   # Plotting UMAPs
   $ Rscript 3_make_plot.R

   # Plotting enrichment boxplots
   $ Rscript 4_cell_hierarchy.R
   ```

---

## Step 1: Variant liftover

**Script:** [`1_liftover.py`](1_liftover.py)

Prepares the fine-mapped HbF GWAS variants for use with the hg19-based scATAC-seq dataset. To run it:
```bash
$ python 1_liftover.py
```

**Input:** `input_data/hbf-META_finemap_20240712_p5e-8.txt` (fine-mapping results, hg38)  
**Output:** `input_data/hbf_meta.finemap.hg19.bed`, `input_data/hbf_meta.finemap.hg19.unmapped.txt`  
**Tools:** Python, UCSC liftOver

---

## Step 2: SCAVENGE TRS computation

**Script:** [`2_scavenge.granja2019.R`](scavenge.granja2019.R)

Computes per-cell trait relevance scores (TRS) for HbF using the SCAVENGE framework.

**Input:** `input_data/hbf_meta.finemap.hg19.bed`, Granja 2019 scATAC-seq RDS  
**Output:** `scav_res/hbf_meta.finemap.trs.txt`, `seurat_objects/` (intermediate Seurat object, ~2.5G)  
**Tools:** R — Seurat, Signac, SCAVENGE, chromVAR, gchromVAR, SummarizedExperiment, BSgenome.Hsapiens.UCSC.hg19, dplyr

### Output file: `scav_res/hbf_meta.finemap.trs.txt`

Tab-separated, one row per cell:

| Column | Description |
|--------|-------------|
| `barcode` | Cell barcode |
| `trs_raw` | TRS without outlier capping |
| `trs_cap_01` | TRS with top 0.1% capped |
| `trs_cap_1` | TRS with top 1% capped |
| `trs_cap_5` | TRS with top 5% capped |
| `celltype` | Cell type label from Granja 2019 (`BioClassification`) |

---

## Step 3: Visualization

**Script:** [`3_make_plot.R`](3_make_plot.R)

Visualizes per-cell TRS scores across 23 annotated hematopoietic cell types (excluding unclassified clusters).

**Input:** `seurat_objects/normal_bm.atac.granja2019.rds`, `scav_res/hbf_meta.finemap.trs.txt`  
**Output:** `plots/atac.granja_2019.celltype.umap.{pdf,png}`, `plots/atac.granja_2019.trs.umap.{pdf,png}`, `plots/color_key.pdf`, `plots/atac.granja_2019.hbf_meta.trs.box.pdf`  
**Tools:** R — Seurat, Signac, ggplot2, RColorBrewer, viridis, pheatmap

---

## Step 4: Cell hierarchy plot (optional)

**Script:** [`4_cell_hierarchy.R`](4_cell_hierarchy.R)

An extended version of `3_make_plot.R` that additionally computes per-cell-type mean TRS values and sets up a schematic hematopoietic hierarchy layout (`coords_df`) for a hierarchy-ordered enrichment plot. The hierarchy plotting code is currently in development.

**Input:** same as Step 3  
**Output:** same UMAP and boxplot PDFs/PNGs as Step 3; hierarchy plot output pending  
**Tools:** R — Seurat, Signac, ggplot2, RColorBrewer, viridis
