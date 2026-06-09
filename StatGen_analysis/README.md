# Statistical Genetics Analysis 

This directory primarily focuses on scripts on all stat gen analyses downstream of GWAS, i.e., meta-analyses, SNP heritability, genetic correlaton, conditional analysis, and multi-ancestry fine-mapping. 

-----
## Table of Content

- [0. Data download](#data-download)
- [1. Fixed-effect Meta Analysis (FEMA) via METAL](#fema)
- [2. SNP Heritability across cohorts](#gencorr-h2)
- [3. Genetic correlations with blood traits in Blood Cell Consortium (BCX2) study](#bcx2-gencorr)
- [4. COJO analysis (via GCTA)](#cojo)
- [5. Multi-Ancestry Meta Analysis (MAMA)](#mama)



## 0. Data download

<a id="data-download"></a>

We have uploaded all our GWAS summary statistics to the GWAS catalog. They could be accessed via:

| Cohort          | Access Code |
|-----------------|-------------|
| St. Jude        |             |
| GTEx            |             |
| Interval        |             |
| Sardinia        |             |
| TopMed          |             |
| BIOS (eQTLgen)  |             |
| Sweden          |             |
| Thailand        |             |
| Tanzania        |             |

Detailed meta info of each cohort could be found on the corresponding GWAS catalog page or Supplementary Table 1 of the manuscript.

### 0.1. Cohort GWAS data harmonization and QC

Summary statistics that were initially in hg19 was lifted over to hg38. The outputs are merged and parsed into tab-delimited pain-text files, and QC-ed with `0a_gwas_qc.py`, with which only variants with below conditions are kept:

*  Single-nucleotide polymorphisms (SNP)
*  Minor allele frequency (MAF) > 0.001
*  Minor allele count (MAC), if available, $> \max(2Nf_{\min},~\min(40, \sqrt{N}), ~ 20)$
*  INFO score > 0.6

The QC logs could be found in `re_gwas_qc_mac40_260416.log`. After QC, the column names are:

| Col Name   |   Content                                                                    |
|------------|------------------------------------------------------------------------------|
| chr        |   Chromosome                                                                 |
| pos        |   Physical position in bp (hg38)                                             |
| a1         |   Reference allele                                                           |
| a2         |   Effect allele                                                              |
| ea_freq    |   Sample frequency of the effect allele (A2)                                 |
| ea_count   |   Sample count of the effect allele (A2)                                     |
| n          |   Sample Size                                                                |
| p          |   P value  (`P_BOLT_LMM` for BOLT-LMM output)                                |
| beta       |   Effect size estimate                                                       |
| se         |   Standard error                                                             |
| z          |   Z statistics (signed squareroot of `CHISQ_BOLT_LMM` for BOLT-LMM output)   |
| info       |   INFO score                                                                 |


## 1. Fixed-effect Meta Analysis (FEMA) via METAL

<a id="fema"></a>

### 1.1. Main analysis

We conducted fixed-effect meta-analysis (FEMA) using METAL (v2020-05-05) to combine GWAS summary statistics across all cohorts. The analysis was performed with sample-size weighted scheme and included genomic control correction.

A METAL configuration file was prepared with the following key parameters:

```bash
SCHEME SAMPLESIZE
SEPARATOR TAB

GENOMICCONTROL ON
AVERAGEFREQ ON
MINMAXFREQ ON

TRACKPOSITIONS OFF
```

The configuration file specified all cohort summary statistics to be included. The input files we used for our analyses could be found at `FEMA/METAL_hbf_all_Nscheme_gcOff_info6.txt`. Note that `TRACKPOSITIONS` is set "OFF" as [sometimes cohorts could be missed when it's set on](https://github.com/statgen/METAL/issues/51).

The meta-analysis was executed with:

```bash
$ ./metal <path_to_config_file>.txt
```
The corresponding running logs could be located at `FEMA/logs/`. The hg38 coordinates were later added to the output `METAL_hbf_inv_all_info6_Nscheme_gcOn_PosTrackOff_1.tbl` through looking up the rsIDs on dbSNP build 157 by running `step1b_metal_to_hg38.sh` on SLURM:

```bash
$ metal_output='/path/to/METAL_hbf_inv_all_info6_Nscheme_gcOn_PosTrackOff_1.tbl'
$ outname='METAL_hbf_inv_ALL_info6_mac40_gcOn_hg38'

$ sbatch --export INPUT=${metal_output},OUTNAME=${outname} step1b_metal_to_hg38.sh
```

The Manhattan plot in Figure 1b was created from `METAL_hbf_inv_ALL_info6_mac40_gcOn_hg38.tsv.gz` through running `FEMA/Fig1B_metal_manhattan.R`. Note that it might take 30--60min for the script to process the file to make sure extreme (<1e-300) p-values don't automatically equate to zero.

To make the qqplot for this meta GWAS, use the Rscript `SupFig1_QQplot.R`

```bash
# Extended Data Fig 1A - META GWAS
$ Rscript SupFig1A_QQplot.R \
  /path/to/METAL_hbf_inv_ALL_info6_mac40_gcOn_hg38.tsv.gz \
  SupFig1A_1 
```
### 1.2 Ancestry-specific analyses

<a id="fema-ances"></a>

METAL was also run on subsets of GWAS cohorts that share the same genetic ancestry (AFR or EUR) under inverse-variance scheme without genomic control to serve as input for MAMA (hg38) and MultiSuSiE (hg19). Their inputs are `METAL_hbf_{afr,eur}_SEscheme_gcOff_info6.txt`. The script to annotate hg19 coordinates is `FEMA/step1b_metal_to_hg19.sh`. (PS: Both shell scripts need RefSNP VCF files for GRCh37 [GCF_000001405.25] and 38 [GCF_000001405.40] under dbSNP build 157.)

### 1.3 Additional FEMA analyses for QC purposes

<a id="fema-misc"></a>

Several stratified and sensitivity meta-analyses were performed in response to several comments from the reviewers:

* **Gene expression-based HbF measurements only** (GTEx + BIOS cohorts: N ≈ 2,588), input `other_inputs/METAL_hbf_eur-hplc_SEscheme_gcOff_info6.txt`. 
* **Liquid chromatography-based HbF measurements only** (Sardinia, St. Jude, Swedish, Thai, TopMed, Interval, Tanzania: N ≈ 25,448), input `other_inputs/METAL_hbf_eur-rnaseq_SEscheme_gcOff_info6.txt`.
* **Leave-One-Batch-Out (LOBO) analyses** - excluding each cohort iteratively to assess robustness, input `other_inputs/METAL_hbf_LOBO-{}_{SE,N}scheme_gcOff_info6.txt`

These stratified analyses allowed us to evaluate consistency of genetic effects across different measurement modalities and ancestry groups. The running logs can be found in `FEMA/logs/`

## 2. Genetic correlation and Heritability across cohorts

<a id="gencorr-h2"></a>

We used **LDAK (v5.2/v6.1)** to estimate SNP-heritability from summary statistics. Key scripts and result files are in `Heritability_and_GenCorr/`.

### 2.0 Data munging for LDAK

FEMA and ancestry-stratified summary statistics were prepared for LDAK `--sum-hers` using the custom script `prepare_ldak_sumstats.py`, which intersects variants with a tagging file and aligns alleles:

```bash
$ python prepare_ldak_sumstats.py \
    --sumstats <FEMA_output>.tsv \
    --tagging <tagging_file>.tagging \
    --mode fema-ss \          # or fema-iv for inverse-variance-weighted
    --snp-col Predictor \
    --out <output>_forLDAK.tsv
```

### 2.1 Preparing LDAK tagging files

<a id="tagging-thai"></a>

#### 2.1.1 Tagging files for Thai cohort

```bash

# Calculate LD tagging file (weights and correlations) for the study
$ ldak --calc-tagging LDAK-Thin.height \ 
       --bfile human --weights weights.thin \
       --power -.25 --window-cm 1 \
       --regression-predictors height.snps

# Extract gene-level definitions from LDAK results using reference annotation file
$ ldak --cut-genes gbat_ldak/results  \
       --genefile tools/RefSeq_GRCh38.txt \
       --bfile /home/jupyter/workspaces/myeloiddiseasect/data/femagctarun/plink-panel-gctardy-dupsgone-maf1/chrALLldak \
       --gene-buffer 500

```
#### 2.1.2 Tagging files for European ancestry

Pre-computed tagging files were downloaded from the LDAK official website. Three models were used:

- `bld.ldak.hapmap.gbr.tagging` — BLD-LDAK model, HapMap3 GBR reference (1,168,975 SNPs)
- `bld.ldak.lite.alpha.hapmap.gbr.tagging` — BLD-LDAK lite-alpha model, HapMap3 GBR
- `ldak.thin.hapmap.gbr.tagging` — LDAK-Thin model, HapMap3 GBR

AFR-specific tagging files (`bld.ldak.hapmap.afr.tagging`, `bld.ldak.genotyped.afr.tagging`) were also used for ancestry-stratified analyses.
<!-- 
### 2.2 Cross-cohort genetic correlations (LDSC)

Pairwise genetic correlations across all 11 cohorts were computed with `ldsc.py --rg` using GRCh38 baseline LD scores (v2.2):

```bash
$ ./ldsc.py \
    --rg Tanzania.LDSCmunged.sumstats.gz,GTEx.LDSCmunged.sumstats.gz,\
TopMed.LDSCmunged.sumstats.gz,Sardinia.LDSCmunged.sumstats.gz,\
StJude.LDSCmunged.sumstats.gz,Sweden.LDSCmunged.sumstats.gz,\
Thai.LDSCmunged.sumstats.gz,Interval.LDSCmunged.sumstats.gz,\
BIOS-RS.LDSCmunged.sumstats.gz,BIOS-LL.LDSCmunged.sumstats.gz,\
BIOS-LLS_660Q.LDSCmunged.sumstats.gz \
    --ref-ld-chr resources/GRCh38/baselineLD_v2.2/baselineLD.@ \
    --w-ld-chr resources/GRCh38/weights/weights.hm3_noMHC.@ \
    --out Heritability_and_GenCorr/HbF_cohorts.GenCorr
```

The output log is `Heritability_and_GenCorr/logs/cohorts.GenCorr.log`. Per-cohort heritability estimates from this run have large standard errors (e.g., Tanzania h² = −0.204 ± 0.576) due to limited per-cohort sample sizes; the pairwise rg estimates are the primary output of interest. -->

### 2.2 SNP-heritability (LDAK)

#### 2.2.1 All ancestry (FEMA results)

Summary statistics from the all-ancestry FEMA (see [Sec 1.1](#fema)) were first formatted for LDAK input by running `prepare_ldak_sumstats.py`. SNP-heritability was then estimated using `--sum-hers` under three tagging models to assess robustness:

```bash
# BLD-LDAK model with HapMap3 GBR LD patterns
$ ./ldak5.2.linux \
    --sum-hers <output_prefix>_bld.ldak.hapmap.gbr \
    --tagfile bld.ldak.hapmap.gbr.tagging \
    --summary CloudTanz_METALOUT_inv_a1alt-munged_forLDAK.tsv \
    --check-sums NO

# BLD-LDAK lite-alpha model
$ ./ldak5.2.linux \
    --sum-hers <output_prefix>_bld.ldak.lite.alpha.hapmap.gbr \
    --tagfile bld.ldak.lite.alpha.hapmap.gbr.tagging \
    --summary CloudTanz_METALOUT_inv_a1alt-munged_forLDAK.tsv \
    --check-sums NO

# LDAK-Thin model
$ ./ldak5.2.linux \
    --sum-hers <output_prefix>_ldak.thin.hapmap.gbr \
    --tagfile ldak.thin.hapmap.gbr.tagging \
    --summary CloudTanz_METALOUT_inv_a1alt-munged_forLDAK.tsv \
    --check-sums NO
```

The three models converged to SNP-heritability estimates for HbF of **h² ≈ 0.085** (BLD-LDAK), **0.051** (BLD-LDAK lite-alpha), and **0.143** (LDAK-Thin). Results are saved under `Heritability_and_GenCorr/output/`.

#### 2.2.2 Thai ancestry

A Thai-specific tagging file was computed in [Sec 2.1.1](#tagging-thai) using genotyped Thai samples. SNP-heritability estimation for the Thai cohort follows the same `--sum-hers` approach, replacing the tagging file and input summary statistics accordingly:

```bash
$ ./ldak5.2.linux \
    --sum-hers <output_prefix>_ldak.thai \
    --tagfile <Thai_tagging_file>.tagging \
    --summary <Thai_GWAS_sumstats_forLDAK>.tsv \
    --check-sums NO
```

#### 2.2.3 European and African ancestry

The EUR-ancestry FEMA summary statistics (see [Sec 1.2](#fema-ances)) were run under the BLD-LDAK model with HapMap3 GBR LD patterns:

```bash
$ ./ldak5.2.linux \
    --sum-hers sumHers_EUR_metal_BLD-LDAK.gbr \
    --tagfile bld.ldak.hapmap.gbr.tagging \
    --summary <EUR_METALOUT>-munged_forLDAK.tsv \
    --check-sums NO
```

The estimated SNP-heritability for EUR-ancestry HbF GWAS is **h² = 0.252 ± 0.047**. Results are saved at `Heritability_and_GenCorr/output/sumHers_EUR_metal_BLD-LDAK.gbr.*`.

Similarly, the AFR-ancestry FEMA results were analyzed using an AFR-specific tagging file, yielding **h² = 0.232 ± 0.258** (large SE due to smaller sample size). Results are at `Heritability_and_GenCorr/output/sumHers_AFR_metal_BLD-LDAK.gbr.*`.



## 3. Genetic correlations with blood traits in Blood Cell Consortium (BCX2) study

<a id="bcx2-gencorr"></a>

Trans-ethnic summary statistics from the BCX2 study were obtained from http://www.mhi-humangenetics.org/en/resources/ and processed for LDAK genetic-correlation analysis. We computed their genetic correlation with HbF level using LDAK (v5.2).

### 3.0 Preparing BCX2 summary statistics for LDAK

BCX2 files are in hg19 with variant IDs in `CHR:POS_REF_ALT` format. They must be (1) lifted over to hg38, (2) annotated with rsIDs from dbSNP for cross-referencing, and (3) reformatted into LDAK's `--summary2` input format (space-delimited, columns: `Predictor A1 A2 n BETA SE P rsID`), where `Predictor = CHR:POS` in hg38 coordinates matching the LDAK tagging files.

This is handled by the stand-alone R script `Heritability_and_GenCorr/bcx2_liftover_addRSID.R`, extracted from the original notebook `BCX2_hg38_liftover_addRSID.ipynb`. The script processes all 15 BCX2 trait files in a single run.

**Dependencies:** `data.table`, `rtracklayer` (Bioconductor), `optparse`

**Usage:**

```bash
$ Rscript bcx2_liftover_addRSID.R \
    --input-dir  /path/to/bcx2/ta/ \
    --output-dir /path/to/bcx2-hg38/ta/ \
    --chain      /path/to/hg19ToHg38.over.chain.gz \
    --dbsnp      /path/to/dbSNP/common_all_20180418_first5cols.txt
```

The chain file can be obtained from UCSC (hg19ToHg38.over.chain.gz). The dbSNP file used here is the "common all" build 150 VCF for GRCh38, subsetted to the first five columns (`#CHROM POS ID REF ALT`). The script prints per-trait statistics on mapping rate and rsID coverage.

**What the script does, step by step:**

1. **Parse hg19 coordinates** — splits the `CHR:POS_REF_ALT` variant ID to extract chromosome and position.
2. **Liftover hg19 → hg38** — uses `rtracklayer::liftOver()` with the supplied chain file; variants with no unique mapping (unmapped or multi-mapped) are dropped and reported.
3. **Annotate rsIDs** — joins on `CHR:POS` (hg38) against the dbSNP reference; unmatched variants retain `NA` in the `rsID` column but are still written to the output.
4. **Write LDAK-formatted output** — one file per trait, named `BCX2_<TRAIT>_Trans_GWAMA.hg38.forLDAK.txt`, with columns: `Predictor A1 A2 n BETA SE P rsID`.

### 3.1 SNP-heritability of HbF (as used in genetic correlation)

The same BLD-LDAK GBR hapmap tagging model applied to the all-ancestry FEMA summary statistics (see [Sec 2.3.1](#)) was used as Trait 1 in the `--sum-cors` command below. Results for Trait 1 heritability from this run are consistent with those reported in Sec 2.3.1.

### 3.2. Genetic correlation with BCX2 blood traits

Genetic correlations between HbF and all 15 BCX2 blood traits were computed using LDAK's `--sum-cors` function under **two tagging models** (BLD-LDAK GBR hapmap and LDAK-Thin GBR hapmap) for robustness. The traits analyzed were: basophil count (BAS), eosinophil count (EOS), hematocrit (HCT), hemoglobin concentration (HGB), lymphocyte count (LYM), mean corpuscular hemoglobin concentration (MCHC), mean corpuscular hemoglobin (MCH), mean corpuscular volume (MCV), monocyte count (MON), mean platelet volume (MPV), neutrophil count (NEU), platelet count (PLT), red blood cell count (RBC), red cell distribution width (RDW), and white blood cell count (WBC).

Example commands for HGB under each model:

```bash
# BLD-LDAK model
$ ./ldak5.2.linux \
    --sum-cors GenCorr/HGB.hapmap.gbr \
    --summary METAL_hbf_inv_all_a1alt-munged_forLDAK.tsv \
    --summary2 BCX2_HGB_Trans_GWAMA.hg38.forLDAK.txt \
    --tagfile bld.ldak.hapmap.gbr.tagging \
    --allow-ambiguous YES \
    --check-sums NO

# LDAK-Thin model
$ ./ldak5.2.linux \
    --sum-cors GenCorr/HGB.thin.hapmap.gbr \
    --summary METAL_hbf_inv_all_a1alt-munged_forLDAK.tsv \
    --summary2 BCX2_HGB_Trans_GWAMA.hg38.forLDAK.txt \
    --tagfile ldak.thin.hapmap.gbr.tagging \
    --allow-ambiguous YES \
    --check-sums NO
```

This was repeated for all 15 traits. Results are stored under `Heritability_and_GenCorr/GenCorr/` as `<TRAIT>.{hapmap,thin}.hapmap.gbr.cors`. The full summary table is at `Heritability_and_GenCorr/all_bcx2_genetic_correlation.tsv`. LDAK-Thin estimates (where BLD-LDAK often gave `nan` due to negative per-category heritability) are reproduced below for reference:

| Trait | rg | SE | z-score | p-value |
|-------|---:|---:|--------:|--------:|
| BAS   | 0.126 | 0.052 | 2.43 | 0.015 |
| EOS   | 0.085 | 0.055 | 1.55 | 0.121 |
| HCT   | 0.050 | 0.105 | 0.48 | 0.632 |
| HGB   | 0.024 | 0.069 | 0.34 | 0.733 |
| LYM   | 0.107 | 0.046 | 2.34 | 0.020 |
| MCH   | −0.223 | 0.193 | −1.16 | 0.247 |
| MCHC  | −0.099 | 0.147 | −0.68 | 0.499 |
| MCV   | −0.220 | 0.174 | −1.26 | 0.206 |
| MON   | 0.078 | 0.035 | 2.23 | 0.026 |
| MPV   | 0.008 | 0.032 | 0.24 | 0.807 |
| NEU   | 0.099 | 0.037 | 2.70 | 0.007 |
| PLT   | −0.005 | 0.144 | −0.04 | 0.970 |
| RBC   | 0.183 | 0.185 | 0.99 | 0.322 |
| RDW   | 0.112 | 0.116 | 0.97 | 0.335 |
| WBC   | 0.115 | 0.043 | 2.68 | 0.007 |

<!-- Positive rg indicates that higher HbF is genetically correlated with higher trait levels. Significant associations (p < 0.05) were found for white blood cell subtypes — BAS, LYM, MON, NEU, and WBC — consistent with a role of fetal haemoglobin in myeloid and lymphoid cell biology. -->
Run script `SuppFig2b_ldak_genCorr.R` to visualize this table and reproduce E.D. Fig 2a:

```bash
$ Rscript SuppFig2b_ldak_genCorr.R
```


## 4. Conditional Analysis (GCTA-COJO)

<a id="cojo"></a>

Conditional and joint association analysis (GCTA-COJO) of metaGWAS results was performed on the All of Us (AoU) Research Workbench v5 to identify conditionally independent loci. EUR and AFR panels were analyzed separately.

### 4.0. Preparing AoU LD reference panel

Per-ancestry, per-chromosome PLINK LD reference panels were built from AoU v5 WGS data using `COJO/0_process_gcta_panels.py`, which runs on an AoU Dataproc cluster with [Hail](https://hail.is/).

For each ancestry (`EUR`, `AFR`, `EAS`, `SAS`) and each autosome (chr1–22), the script:
1. Removes QC-flagged samples (`flagged_samples.tsv`) and related individuals (`relatedness_flagged_samples.tsv`)
2. Retains only samples matching the target ancestry (`ancestry_preds.tsv`)
3. Removes sites flagged by AoU site-level QC (`mt.filters`)
4. Splits multi-allelic variants
5. Applies variant QC filters: HWE p > 1×10⁻¹⁵, allele frequency > 0.001, call rate > 0.90
6. Exports per-chromosome PLINK `.bed/.bim/.fam` files to `gs://$WORKSPACE_BUCKET/{ancestry}-plink-AF1pct/chr{N}_filt{ANCESTRY}_p1pct_gcta`

```bash
# Run on AoU Dataproc cluster (requires PySpark + Hail)
python COJO/0_process_gcta_panels.py --ancestry EUR
python COJO/0_process_gcta_panels.py --ancestry AFR
```

The EUR panel contained 51,125 individuals and were randomly downsampled to 10,000. Panels were generated for all 22 autosomes.

### 4.1. Parsing METAL output for GCTA-COJO

METAL fixed-effect meta-analysis (FEMA) output was converted to GCTA-COJO format using `COJO/1_prepare_cojo_input.R`. The script:
1. Extracts the original reference allele from the variant ID (`chr:pos:ref:alt` format)
2. Filters variants with allele frequency < 0.01 or > 0.99 (MAF ≥ 0.01)
3. Orients alleles to the original reference: if Allele1 matches the reference, keep freq; otherwise flip to `1 − freq`
4. Computes effect size and standard error from Z-scores and effective sample size (N):

$$\beta = \frac{Z}{\sqrt{2 \cdot f \cdot (1-f) \cdot (N + Z^2)}}, \quad \text{SE} = \frac{1}{\sqrt{2 \cdot f \cdot (1-f) \cdot (N + Z^2)}}$$

Output columns (tab-delimited): `SNP  A1  A2  freq  b  se  p  N`

### 4.2. Running GCTA-COJO

GCTA-COJO stepwise selection (`--cojo-slct`) was run using GCTA v1.94.0 Beta with a p-value threshold of 1×10⁻⁶. EUR and AFR panels were analyzed in parallel, looping over all 22 autosomes:

```bash
p_thresh=1e-6
gcta="./tools/gcta_v1.94.0Beta_linux_kernel_3_x86_64/gcta_v1.94.0Beta_linux_kernel_3_x86_64_static"
panel_dir="/path/to/eur-plink-AF1pct"   # or afr-plink-AF1pct
cojo_input="CloudTanz_METALOUT_inv_a1alt-munged.cojoReady.tsv"

for chr in {1..22}; do
  $gcta \
    --bfile ${panel_dir}/chr${chr}_filtEUR_p1pct_gcta \
    --cojo-slct \
    --cojo-p ${p_thresh} \
    --cojo-wind 10000 \
    --cojo-collinear 0.9 \
    --diff-freq 0.6 \
    --chr ${chr} \
    --maf 0.001 \
    --cojo-file ${cojo_input} \
    --out ./output/eur_plink/chr${chr}.${p_thresh}
done
```

Key parameters:
- `--cojo-wind 10000`: 10 Mb window for conditional analysis
- `--cojo-collinear 0.9`: collinearity threshold for stepwise selection
- `--diff-freq 0.6`: maximum allele frequency difference between GWAS and LD reference
- `--maf 0.001`: minimum MAF in LD reference panel

Per-chromosome `.jma.cojo` files were concatenated into the final outputs:
- **`eur_plink_1e6-eur.jma.cojo`** — 131 conditionally independent loci (EUR panel)
- **`afr_plink_1e6-afr.jma.cojo`** — 253 conditionally independent loci (AFR panel)

Output columns: `Chr  SNP  bp  refA  freq  b  se  p  n  freq_geno  bJ  bJ_se  pJ  LD_r`  
where `bJ`/`bJ_se`/`pJ` are the jointly estimated effect, SE, and p-value.

### 4.3. Post-processing

Two R scripts filter inflated conditional Z-scores and merge EUR/AFR results:

**`COJO/3_filter_inflated_zJ.R`** — filters loci where the conditional Z-score (`bJ/bJ_se`) is inflated relative to the marginal Z-score (`b/se`). Loci with `|ZJ/Z| > (1 + fold-change cutoff)` are removed. Usage:
```bash
Rscript COJO/3_filter_inflated_zJ.R <cojo_file> <foldchange_cutoff> <output_prefix>
```

**`COJO/3_merge_eur-afr_cojos.R`** — merges EUR and AFR COJO results, applies the same fold-change filter to each ancestry, and annotates each locus with the nearest coding gene using a UCSC RefSeq gene reference (`hg38_ncbiRefSeq_AllGenes_2026UCSC.tsv.gz`). Usage:
```bash
Rscript COJO/3_merge_eur-afr_cojos.R \
  eur_plink_1e6-eur.jma.cojo \
  afr_plink_1e6-afr.jma.cojo \
  <output_prefix> \
  <foldchange_cutoff> \
  [ref_genes.tsv.gz]
```



## 5. Multi-Ancestry Meta Analysis (MAMA)

<a id="mama"></a>
<!-- 
(content in `MAMA.ipynb` )
Section 5.0 Data Munging -->

### 5.0 Preparing plink LD panels

We assembled a reference LD panel merging variant calls from whole genome sequences of 197 Thai, 660 1000GenomesProject (1KG) EUR, and 891 1KG AFR individuals. Their IDs are recorded in `MAMA/input/ancespath-file.txt`.

We first merged these samples per-chromosome and recorded all SNPs that have high missing rates. These SNPs were recorded in `MAMA/warned_variants_1stmerge/`.

Then, within samples of each ancestry (EUR, AFR, or THAI), these variants were excluded via `--exclude MAMA/bad_variants_1stmerge/eur_afr_thai_merged_chr*-merge.missnp`. Afterwards, plink2 was deployed to perform per-chromosome variant QC:
```bash
# $ances could be "eur" or "afr"
# Remove warned variants
$ plink \
--exclude bad_variants_1stmerge/eur_afr_thai_merged_chr${chr}-merge.missnp \
    --keep-allele-order \
    --make-bed \
    --out MAMA/PLINK_panel/${ances}/gnomad_${ances}_chr${chr}

# For each ancestry {eur,afr,thai} & each chromosome (chr):
$ plink2 \
--bfile MAMA/PLINK_panel/${ances}/gnomad_${ances}_chr${chr} \
--snps-only \
--autosome \
--maf 0.001 \
--hwe 1e-15 \
--geno 0.02 \
--make-bed --out MAMA/PLINK_panel/${ances}/${ances}_chr${chr}
```

After generating cleaned bedfiles of all chromosomes in all ancestries, their .bim files are pooled to generate a SNP ancestry file `snp-ances-file.txt` that looks like:

```text
SNP	AFR	EUR	THAI
1:14464:A:T	1	0	0
1:14717:G:A	0	0	1
1:15616:G:A	0	1	0
1:17005:A:G	0	0	1
1:17385:G:A	0	0	1
1:17398:C:A	0	0	1
1:17407:G:A	0	0	1
1:17479:G:A	0	1	0
1:17556:C:T	0	0	1
...
```
Suppose all bedfiles are saved under path `MAMA/PLINK_panel/${ances}/${ances}_chr{1..22}.{bim,bam,bed}`, then this SNP ancestry file could be generated by running
```bash
$ python MAMA/0_compile_bims_to_snpfile.py
```
The file `MAMA/input/snp-ances-file.txt.zst` is the compressed output. One could decompress it using `zstd -d <file>` command in a Unix environment.

Meanwhile, for each chromosome, bedfiles of all three ancestries are merged into one using plink for LD score generation.

### 5.1 Step 1
```bash
# by each chromosome $CHR
$ python /path/to/mama/mama_ldscore.py \
        --out MAMA/by-chr-ldscores-kb/mama_1kg-thai_1000kb_chr${chr} \
        --stream-stdout  \
        --snp-ances MAMA/input/snp-ances-file.txt \
        --ances-path MAMA/input/ancespath-file.txt \
        --bfile-merged-path MAMA/PLINK_panel/merged/1kg-thai_chr${chr} \
        --std-geno-ldsc  \
        --ld-wind-kb 1000
```

### 5.2 Step 2

To generate harmonized input for Step 2, run the script `MAMA/2_make_mama_input.py`, which will read the METAL outputs from [Sec 1.2](#fema-ances) for EUR and AFR, as well as the cleaned sumstats for Thai cohort, harmonize variants across the three, and polarize the SNPs according to labels in `snp-ances-file.txt`.



```bash
# File paths
$ INPUT_PREFIX="MAMA/input/thai-eurFema-afrFema_SNPaligned"
$ MAMA_LDSC="MAMA/by-chr-ldscores-kb/mama_1kg-thai_1000kb_chr*.l2.ldscore.gz"
$ OUTNAME="MAMA_SNPaligned"

$ python mama.py \
        --sumstats "${INPUT_PREFIX}.EUR.tsv,EUR,HBF" \
        "${INPUT_PREFIX}.AFR.tsv,AFR,HBF" \
        "${INPUT_PREFIX}.THAI.tsv,THAI,HBF" \
        --ld-scores "${MAMA_LDSC}" \
        --out "MAMA/${OUTNAME}"\
        --use-standardized-units \
        --reg-int-zero \
        --reg-ld-set-corr 1.0 \
        --out-harmonized \
        --freq-bounds 0.0001 0.9999 \
        --verbose
```



### 5.3. Plotting

Manhattan plots for each ancestry group were generated using the R script `MAMA/Fig1cde_MAMA_manhattan.R`. The script writes both PNG and PDF outputs, standardizes several common summary-stat column names, and by default annotates lead peaks with the nearest gene. Peak annotation can be disabled with `--skip-annotation`.

To run it, use:

```bash
$ Rscript Fig1cde_MAMA_manhattan.R <input_file> <output_prefix> <ancestry_label> [figwidth] [figheight] [annotation_cutoff] [--skip-annotation]
```

Optional arguments:

* `figwidth`, `figheight`: figure size in inches. Defaults are `10` and `6`.
* `annotation_cutoff`: threshold for peak annotation. Accepted values are `sig` / `significant` / `5e-8` for genome-wide significant peaks, or `sugg` / `suggestive` / `sug` / `1e-6` for suggestive peaks. Default is `sig`.
* `--skip-annotation`: suppress nearest-gene labels.

To recreate Figure 1c, d, and e, run (pseudocode):
```bash
# Figure 1C - African ancestry
$ Rscript 3_Fig1cde_MAMA_manhattan.R \
        /path/to/XC_SNPaligned_16865_AFR_HBF.res \
        Fig1C-AFR \
        AFR \
        5 \
        4 

# Figure 1D - European ancestry
$ Rscript 3_Fig1cde_MAMA_manhattan.R \
        /path/to/XC_SNPaligned_16865_EUR_HBF.res \
        Fig1D-EUR \
        EUR \
        5 \
        4 

# Figure 1E - Thai ancestry 
$ Rscript 3_Fig1cde_MAMA_manhattan.R \
        /path/to/XC_SNPaligned_16865_THAI_HBF.res \
        Fig1E-THAI \
        THAI \
        5 \
        4 

# Skip automatic peak annotation
$ Rscript 3_Fig1cde_MAMA_manhattan.R \
    /path/to/XC_SNPaligned_16865_EUR_HBF.res \
    Fig1D-EUR-no-labels \
    EUR \
    --skip-annotation
```

To re-create their perspective qqplots, run
```bash
# Extended Data Fig 1B - MAMA AFR
$ Rscript 4_SupFig1_QQplot.R \
    /path/to/XC_SNPaligned_16865_AFR_HBF.res \
    SupFig1B AFR

# Extended Data Fig 1C - MAMA EUR
$ Rscript 4_SupFig1_QQplot.R \
    /path/to/XC_SNPaligned_16865_EUR_HBF.res \
    SupFig1C EUR

# Extended Data Fig 1D - MAMA Thai
$ Rscript 4_SupFig1_QQplot.R \
    /path/to/XC_SNPaligned_16865_THAI_HBF.res \
    SupFig1D THAI

```

To plot the correlation between their effect sizes:

```bash
$ Rscript SupFig1_betaCorr.R \
    XC_SNPaligned_EUR_HBF.res \
    XC_SNPaligned_AFR_HBF.res \
    SupFig1e EUR AFR --sig-only

$ Rscript SupFig1_betaCorr.R \
    XC_SNPaligned_EUR_HBF.res \
    XC_SNPaligned_THAI_HBF.res \
    SupFig1f EUR Thai --sig-only

$ Rscript SupFig1_betaCorr.R \
    XC_SNPaligned_AFR_HBF.res \
    XC_SNPaligned_THAI_HBF.res \
    SupFig1E AFR Thai --sig-only
```


## 6. Multi-Ancestry Finemapping of BACH2 locus

<a id="finemap"></a>

Multi-ancestry statistical fine-mapping of the BACH2 locus was performed using [MultiSuSiE](https://github.com/jordanlab/MultiSuSiE) (RSS formulation), combining EUR, AFR, and Thai GWAS summary statistics with ancestry-matched LD matrices. 

All fine-mapping inputs use **hg19/GRCh37** coordinates. Computing was performed on the **E3** BCH HPC cluster (SLURM). Scripts contain E3- and lab-share–specific paths; adapt these to your own environment.

### 6.0 Downloading Pan-UKBB LD matrices

Pan-UKBB reference LD sub-matrices for the BACH2 region were extracted from the public S3 release using Hail (`MultiSuSiE/scripts/0_extract_hail_LDs.py`). The script supports region lookup either by explicit coordinates or by gene name (resolved against a local UCSC RefSeq reference):

```bash
# By explicit hg19 coordinates
python MultiSuSiE/scripts/0_extract_hail_LDs.py \
    --chrom 6 --start 90000000 --end 92000000 \
    --output-dir /path/to/ld/ref_ld 

# By gene name with a flanking window
python MultiSuSiE/scripts/0_extract_hail_LDs.py \
    --gene BACH2 --flank 1Mb \
    --output-dir /path/to/ld/ref_ld 
```

Reads from S3:
- `s3://pan-ukb-us-east-1/ld_release/UKBB.{POP}.ldadj.variant.ht`
- `s3://pan-ukb-us-east-1/ld_release/UKBB.{POP}.ldadj.bm`

Outputs per population (`EUR`, `AFR`, `EAS`[as proxy for Thai]):
- `UKBB.{POP}.{chrom}_{start}_{end}.ld.npy` — NumPy LD matrix
- `UKBB.{POP}.{chrom}_{start}_{end}.variants.tsv.gz` — variant metadata

Note that the running environment must have Hail and Spark enabled.

### 6.1 Liftover Thai GWAS summary statistics to hg19

Thai cohort GWAS summary statistics (hg38) were lifted over to hg19 for fine-mapping using `liftOver` and validated against dbSNP build 157 GRCh37.

**Step 1 — Coordinate liftover** (`MultiSuSiE/scripts/1_liftover_thai_gwas.sh`, SLURM):
```bash
# Prepare a BED file of hg38 positions and run liftOver
zcat Thai_chrALL_MAF0p1pct_hg38_info6.clean.tsv.gz \
  | awk 'BEGIN{OFS="\t"} NR>1{print "chr"$2, $3-1, $3}' \
  > raw/hg19_sumstats/Thai_gwas_hg38_info6.bed

liftOver -minMatch=0.98 -multiple \
    raw/hg19_sumstats/Thai_gwas_hg38_info6.bed \
    hg38ToHg19.over.chain.gz \
    raw/hg19_sumstats/Thai_gwas_hg19_info6.bed \
    raw/hg19_sumstats/Thai_gwas_info6_unmapped.bed
```

**Step 2 — Merge and validate** (`MultiSuSiE/scripts/1_merge_liftover_coord.py`):
Merges liftover BED coordinates back into the GWAS summary stats, performs sanity checks (non-chromosomal contigs, missing matches), and batch-queries dbSNP build 157 via `bcftools` to validate lifted positions and backfill missing rsIDs. Variants failing coordinate validation are written to a filtered-variants log.

```bash
python MultiSuSiE/scripts/1_merge_liftover_coord.py \
    --gwas-file  raw/hg38_sumstats/Thai_chrALL_MAF0p1pct_hg38_info6.clean.tsv.gz \
    --bed-file   raw/hg19_sumstats/Thai_gwas_hg19_info6.bed \
    --output-file raw/hg19_sumstats/Thai_gwas_hg19_info6_dbsnpValidated.tsv.gz \
    --filtered-file raw/hg19_sumstats/Thai_gwas_hg19_info6_filtered.tsv.gz
```

### 6.2 Constructing fine-mapping inputs

Aligned Z-score files and LD matrices for MultiSuSiE were built using `MultiSuSiE/scripts/2_prepare_{2,3}way_matrices_n_diagnostics.py`, called via the SLURM wrapper `MultiSuSiE/scripts/2b_prepare_matrices_n_diagnostics.sh`.

Two analysis configurations were prepared:

**2-way (EUR + AFR)** — uses Pan-UKBB reference LD for both populations:
```bash
python MultiSuSiE/scripts/2_prepare_2way_matrices_n_diagnostics.py \
    --gene BACH2 --flank 1Mb \
    --eur-sumstats fema/hg19_output/FEMA_EUR_hg19.tsv.gz \
    --afr-sumstats fema/hg19_output/FEMA_AFR_hg19.tsv.gz \
    --eur-ld       ld/ref_ld/UKBB.EUR.6_<start>_<end>.ld.npy \
    --afr-ld       ld/ref_ld/UKBB.AFR.6_<start>_<end>.ld.npy \
    --output-prefix loci/finemap_input/BACH2_2way
```

**3-way (EUR + AFR + Thai)** — use EAS LD from PanUKB as a proxy for Thai:
```bash
python MultiSuSiE/scripts/2_prepare_3way_matrices_n_diagnostics.py \
    --gene BACH2 --flank 1Mb \
    --eur-sumstats  fema/hg19_output/FEMA_EUR_hg19.tsv.gz \
    --afr-sumstats  fema/hg19_output/FEMA_AFR_hg19.tsv.gz \
    --thai-sumstats raw/hg19_sumstats/Thai_gwas_hg19_info6_dbsnpValidated.tsv.gz \
    --output-prefix loci/finemap_input/BACH2_3way
```

Both scripts resolve the BACH2 locus from a UCSC RefSeq gene reference, intersect variants across cohorts and LD sources, align alleles, and produce per-population output pairs:
- `<prefix>.{EUR,AFR,THAI}.z` — tab-delimited Z-score file (`rsid chromosome position allele1 allele2 maf beta se n z`)
- `<prefix>.{EUR,AFR,THAI}.ld` — LD correlation matrix (`.npy` or space-delimited)

Diagnostic plots (LD heatmaps, Z-score distributions, cross-population Manhattan) are also written.

### 6.3 Running MultiSuSiE

Multi-ancestry fine-mapping was run using the [MultiSuSiE](https://github.com/jordanlab/MultiSuSiE) `multisusie_rss` function. Separate scripts handle the 2-way and 3-way configurations:

```bash
# EUR + AFR (2-way)
python MultiSuSiE/scripts/3a_run_2way_multisusie.py \
    --prefix loci/finemap_input/BACH2_2way \
    --out    multisusie/BACH2_2way

# EUR + AFR + Thai (3-way)
python MultiSuSiE/scripts/3a_run_3way_multisusie.py \
    --prefix loci/finemap_input/BACH2_3way \
    --out    multisusie/BACH2_3way
```

Each script loads the `.z` and `.ld` files for all populations, runs `multisusie_rss`, and writes:
- `*_credible_sets.tsv` — credible set membership
- `*_pip.tsv` — per-variant posterior inclusion probabilities (PIPs)
- `*_effects.tsv` (per population) — posterior effect sizes

Optional single-ancestry and alternative fine-mapping runs are provided as SLURM batch scripts:
- `MultiSuSiE/scripts/3b_run_susie.sbatch` — SuSiE (single-ancestry)
- `MultiSuSiE/scripts/3b_run_finemap.sbatch` — FINEMAP
- `MultiSuSiE/scripts/3b_run_carma.sbatch` — CARMA

### 6.4 Plotting

Manuscript figures for the BACH2 locus were generated by four R scripts in `MultiSuSiE/scripts/`, which all source shared utility functions from `step4_utility.R`:

| Script | Figure | Description |
|---|---|---|
| `4_plot_fig2a.R` | Fig. 2a | Regional fine-mapping plot: −log₁₀p + PIP tracks + ATAC-seq signal per cell type (UCSC hg19 coords) |
| `4_plot_fig2a_wMAMA.R` | Fig. 2a variant | Same as above but overlaying MAMA p-values alongside FEMA |
| `4_plot_regional_pval.R` | Regional p-value | −log₁₀p locus plot for EUR/AFR/THAI with optional credible-set and LD coloring |
| `4_plot_regional_MAMA-pval.R` | MAMA regional | Regional p-value plot using MAMA summary statistics |
| `4_plot_variant_tiles.R` | Variant tile plot | Per-variant credible-set tile visualization across ancestries |

`4_plot_fig2a.R` and `4_plot_fig2a_wMAMA.R` integrate ATAC-seq tracks from two published data sources (Corces 2016 or Weng 2024), selectable via the `--atac-source` argument. Example usage:

```bash
Rscript MultiSuSiE/scripts/4_plot_fig2a.R \
    loci/finemap_input/BACH2_3way \   # z-file prefix
    6 90500000 91500000 \              # chrom start end (hg19)
    multisusie/BACH2_3way \           # SuSiE result prefix
    BACH2 \                           # gene name for annotation
    "CD14,CD34,HSC" \                 # cell types to plot (comma-sep, or leave blank for all)
    "rs12203592,rs72934505" \         # credible-set rsIDs to label
    rs12203592 \                      # leading SNP for LD coloring
    Weng2024                          # ATAC source
```

