# CUT&RUN Analysis

Scripts for processing CUT&RUN chromatin profiling data. All scripts take a sample base name as the first positional argument.

## Overview

The pipeline maps paired-end CUT&RUN reads to hg38, uses *E. coli* spike-in reads for cross-sample calibration, calls peaks, generates bigWig tracks, and then scans peak sequences for transcription factor motif enrichment. It was used to profile BACH2 and NRF2 occupancy in erythroid cells.

## How to use

Run the four scripts in order for each sample:

```bash
# Step 0: Trim adapters and normalize read length
$ bash 0_trim_fastq.sh <sample_name>

# Step 1: Map to E. coli spike-in genome (for normalization)
$ bash 1_run.spikein.mapping.sh <sample_name>

# Step 2: Map to hg38, call peaks, generate bigWig tracks
$ bash 2_run.cutrun.sh <sample_name>

# Step 3: Scan peaks for TF motif enrichment
$ bash 3_motif_calling.sh <sample_name>
```

Steps 0–2 should be run for all samples before Step 3. After Step 1, collect the per-sample spike-in read counts from `spikein.read` and compute normalization scale factors (e.g., `min(spike-in counts) / sample spike-in count`). These scale factors are then used to downsample the peak files before running Step 3 — the downsampled files (`<sample>.downsampled_peaks.narrowPeak`) are the expected input for `3_motif_calling.sh`.

---

## [`0_trim_fastq.sh`](0_trim_fastq.sh) — Adapter trimming and read length normalization

```bash
$ bash 0_trim_fastq.sh <sample_name>
```

Two-stage trimming of raw paired-end reads:

1. **Trimmomatic (v0.36)** — Removes TruSeq3 PE adapters (`ILLUMINACLIP:Truseq3.PE.fa:2:15:4:4:true`), quality-trims both ends (LEADING/TRAILING 20, SLIDINGWINDOW 4:15), and discards reads shorter than 25 bp. Paired and unpaired reads are written separately to `trimmed/`.
2. **Length normalization** — Trims reads to a fixed length of 150 bp using `kseq_python.py`. Output goes to `trimmed3/`, which is used by all downstream steps.

| | |
|---|---|
| **Input** | `<HOME>/group/fastq/<sample>_R1/R2_001.fastq.gz` |
| **Output** | `trimmed3/<sample>_1/2.paired.fastq.gz` |
| **Tools** | Trimmomatic 0.36, Python (`kseq_python.py`) |

---

## [`1_run.spikein.mapping.sh`](1_run.spikein.mapping.sh) — Spike-in mapping

```bash
$ bash 1_run.spikein.mapping.sh <sample_name>
```

Maps trimmed reads to the *E. coli* K-12 MG1655 genome to quantify spike-in reads for calibration. The number of properly paired spike-in reads (`samtools view -c -F 12`) is appended to `spikein.read` for all samples. After running this for all samples, use the spike-in counts to compute per-sample normalization scale factors before proceeding to Step 3.

| | |
|---|---|
| **Input** | `trimmed3/<sample>_1/2.paired.fastq.gz` |
| **Output** | `spikein/<sample>_aligned_reads.bam`; spike-in read counts appended to `spikein.read` |
| **Tools** | bowtie2, samtools |

---

## [`2_run.cutrun.sh`](2_run.cutrun.sh) — Alignment, peak calling, and track generation

```bash
$ bash 2_run.cutrun.sh <sample_name>
```

Main processing pipeline, aligned to hg38:

1. **Alignment** — Bowtie2 (`--dovetail`, 10 threads) maps trimmed reads to hg38; raw BAM written to `aligned.aug10/`.
2. **Sorting and indexing** — `samtools sort` + `samtools index` produce a sorted BAM in `sorted/`.
3. **Duplicate removal** — Picard `MarkDuplicates` (`REMOVE_DUPLICATES=true`) writes a deduplicated BAM to `dedup/`.
4. **Peak calling** — MACS2 `callpeak` on the deduplicated BAM (`-g hs -f BAMPE -q 0.01 --SPMR --keep-dup all`). Peak files written to `peak/`.
5. **Monoclonal read extraction** — Pairs are converted to BEDPE, collapsed to unique fragments, and deduplicated via `sorted_bed_merge_redundant_lines_V2.pl`. Output BED written to `monoclonal/`.
6. **Genome coverage** — `bedtools genomecov -bg` computes fragment coverage across hg38 chromosome sizes.
7. **BigWig track** — `bedGraphToBigWig` converts the bedgraph to a bigWig file in `bigwig/` for visualization in IGV or the UCSC Genome Browser.

| | |
|---|---|
| **Input** | `trimmed3/<sample>_1/2.paired.fastq.gz` |
| **Output** | `dedup/<sample>.bam`, `peak/<sample>_peaks.narrowPeak`, `bigwig/<sample>.bw` |
| **Tools** | bowtie2, samtools, Picard, MACS2, bedtools, bedGraphToBigWig |

---

## [`3_motif_calling.sh`](3_motif_calling.sh) — Motif enrichment analysis

```bash
$ bash 3_motif_calling.sh <sample_name>
```

Scans CUT&RUN peaks for TF motif enrichment using FIMO. Expects spike-in–normalized, downsampled narrow peak files as input (see note in [How to use](#how-to-use) above).

1. **Sequence extraction** — `bedtools getfasta` extracts peak sequences from the hg38 FASTA; `fix_sequence.py` cleans the output FASTA.
2. **FIMO motif scanning** — For every `.meme` file in `motif_for_use/`, FIMO is run at p < 0.05 with genomic coordinate parsing (`--parse-genomic-coord`). Results written per motif to `motif_calling/<sample>/fimo2.<motif>/`.
3. **BED conversion** — FIMO's GFF output is converted to BED format via `gff2bed` for downstream intersection analyses.

| | |
|---|---|
| **Input** | `peak/<sample>.downsampled_peaks.narrowPeak` |
| **Output** | `motif_calling/<sample>/fimo2.<motif>/fimo.bed` for each motif |
| **Tools** | FIMO (MEME Suite), bedtools, BEDOPS (`gff2bed`) |
