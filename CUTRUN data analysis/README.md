## CUT&RUN Data Analysis Pipeline

This section describes the CUT&RUN data analysis workflow implemented in this repository.

### Step 0: Adapter Trimming and Quality Control

**Script:** `0_trim_fastq.sh`  

- Perform adapter trimming and basic quality control on raw CUT&RUN FASTQ files.
- Remove low-quality bases and short reads prior to alignment.

**Input:** raw FASTQ files  
**Output:** trimmed FASTQ files  
**Tools:** trimmomatic 

---

### Step 1: Spike-in Mapping and Normalization

**Script:** `1_run.spikein.mapping.sh`  

- Align reads to the spike-in reference genome *E. coli*.

**Input:** trimmed FASTQ files  
**Output:** spike-in–mapped BAM files  
**Tools:** bowtie2, samtools  

---

### Step 2: CUT&RUN Alignment and Peak Calling

**Script:** `2_run.cutrun.sh`  

- Align CUT&RUN reads to the reference genome *hg38*.
- Filter reads and generate normalized coverage tracks.
- Call peaks for CUT&RUN signal.

**Input:** trimmed FASTQ files  
**Output:** filtered BAM files, bigWig tracks, peak files  
**Tools:** bowtie2, samtools, picard, MACS2  

---

### Step 3: Motif Analysis

**Script:** `3_motif_calling.sh`  

- Perform motif calling analysis on CUT&RUN peaks.

**Input:** CUT&RUN peak files  
**Output:** motif enrichment results and motif logos  
**Tools:** MEME / FIMO  

---
