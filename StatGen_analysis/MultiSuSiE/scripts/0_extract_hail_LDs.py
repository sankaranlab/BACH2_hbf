#!/usr/bin/env python3
"""
step2a_extract_hail_LDs.py

Stand-alone script to download a Pan-UKBB regional LD sub-matrix from the
public S3 release and save it locally for offline fine-mapping use.

Reads from S3:
  s3://pan-ukb-us-east-1/ld_release/UKBB.{POP}.ldadj.variant.ht
  s3://pan-ukb-us-east-1/ld_release/UKBB.{POP}.ldadj.bm

Saves locally:
    {output_dir}/UKBB.{POP}.{chrom}_{start}_{end}.ld.npy       - numpy LD matrix
    {output_dir}/UKBB.{POP}.{chrom}_{start}_{end}.variants.tsv.gz - variant metadata

Usage examples:
  # by explicit coordinates
  python fetch_ukbb_ld.py --chrom 6 --start 90000000 --end 92000000 \\
      --output-dir /data/ld_cache --gcp-project my-gcp-project

  # by gene name (requires gene reference file)
  python fetch_ukbb_ld.py --gene BACH2 --flank 1Mb \\
      --output-dir /data/ld_cache --gcp-project my-gcp-project
"""

import argparse
import gzip
import os
import time

import hail as hl
from hail.linalg import BlockMatrix
import numpy as np
import pandas as pd


# paths

ROOT       = "/lab-share/Hem-Sankaran-e2/Public"
b37_GENE_REF = os.path.join(ROOT, "ref_genomes/human/hg19_ncbiRefSeq_RefSeqAll_2026UCSC.tsv.gz")
b38_GENE_REF = os.path.join(ROOT, "ref_genomes/human/hg38_ncbiRefSeq_AllGenes_2026UCSC.tsv.gz")


# helpers
def normalize_chr(value):
    c = str(value).replace("chr", "").replace("CHR", "")
    if c == "M":
        c = "MT"
    return c


def parse_flank_size(value):
    value = value.strip()
    units = {"b": 1, "kb": 1e3, "mb": 1e6, "gb": 1e9, "k": 1e3, "m": 1e6, "g": 1e9}
    try:
        return int(float(value))
    except ValueError:
        pass
    vl = value.lower()
    for unit, mult in units.items():
        if vl.endswith(unit):
            try:
                return int(float(vl[: -len(unit)].strip()) * mult)
            except ValueError:
                pass
    raise argparse.ArgumentTypeError(
        f"Cannot parse flank size '{value}'. Use: 500000, 1e6, 1Mb, 500kb, 2m"
    )


def get_gene_coordinates(gene_name, gene_ref_path=b37_GENE_REF):
    print(f"Looking up coordinates for gene: {gene_name}")
    if not os.path.exists(gene_ref_path):
        raise FileNotFoundError(f"Gene reference file not found: {gene_ref_path}")
    open_func = gzip.open if gene_ref_path.endswith(".gz") else open
    gene_records = []
    with open_func(gene_ref_path, "rt") as f:
        f.readline()  # skip header
        for line in f:
            fields = line.strip().split("\t")
            if len(fields) < 13 or fields[12] != gene_name:
                continue
            chrom = fields[2].replace("chr", "")
            if chrom not in [str(i) for i in range(1, 23)] + ["X", "Y"]:
                continue
            gene_records.append((chrom, int(fields[4]), int(fields[5])))
    if not gene_records:
        raise ValueError(f"Gene '{gene_name}' not found in {gene_ref_path}")
    chrom = gene_records[0][0]
    print(f"Found {len(gene_records)} transcript(s) for {gene_name} on chr{chrom}")
    return chrom, min(r[1] for r in gene_records), max(r[2] for r in gene_records)


# core fetching

def fetch_pop_ld(pop, chrom, start, end, output_dir):
    """
    Fetch the regional LD sub-matrix for one population from S3.

    Skips download if local output files already exist.

    Returns:
        (matrix_path, variants_path) - absolute paths to the saved files.
    """
    tag           = f"{normalize_chr(chrom)}_{start}_{end}"
    matrix_path   = os.path.join(output_dir, f"UKBB.{pop}.{tag}.ld.npy")
    variants_path = os.path.join(output_dir, f"UKBB.{pop}.{tag}.variants.tsv.gz")

    if os.path.exists(matrix_path) and os.path.exists(variants_path):
        nv = len(pd.read_csv(variants_path, sep="\t", compression="gzip"))
        print(f"[{pop}] Local cache already exists ({nv} variants) — skipping download.")
        return matrix_path, variants_path

    chrom_str = normalize_chr(chrom)
    # chrom_alt = f"chr{chrom_str}"
    interval = hl.parse_locus_interval(f"{chrom_str}:{start}-{end}")

    # step 1: read and filter the variant table
    print(f"[{pop}] Reading variant table from S3 ...")
    ht = hl.read_table(f"s3://pan-ukb-us-east-1/ld_release/UKBB.{pop}.ldadj.variant.ht")
    ht = ht.filter(interval.contains(ht.locus))
    # ht = ht.filter(
    #     ((ht.locus.contig == chrom_str) | (ht.locus.contig == chrom_alt))
    #     & (ht.locus.position >= int(start))
    #     & (ht.locus.position <= int(end))
    # )
    # Keep a deterministic order so matrix rows/cols align with the variant table.
    ht = ht.order_by(ht.idx)

    row_fields = set(ht.row_value.dtype.fields.keys())
    if "rsid" in row_fields:
        snp_expr = hl.str(ht.rsid)
    elif "variant" in row_fields:
        snp_expr = hl.str(ht.variant)
    else:
        snp_expr = hl.str(ht.locus)

    ht = ht.annotate(
        chr=ht.locus.contig,
        pos=ht.locus.position,
        snp=snp_expr,
        a1=hl.if_else(hl.len(ht.alleles) > 0, ht.alleles[0], hl.missing(hl.tstr)),
        a2=hl.if_else(hl.len(ht.alleles) > 1, ht.alleles[1], hl.missing(hl.tstr)),
    )
    ht = ht.select("idx", "chr", "pos", "snp", "a1", "a2")

    idx_list = ht.idx.collect()
    if len(idx_list) == 0:
        raise RuntimeError(
            f"[{pop}] No variants found in region {chrom_str}:{start}-{end}. "
            "Check chromosome name and coordinates."
        )
    print(f"[{pop}] {len(idx_list)} variants in region.")

    # step 2: extract the BlockMatrix sub-matrix
    print(f"[{pop}] Fetching BlockMatrix from S3 (this may take several minutes) ...")
    bm     = BlockMatrix.read(f"s3://pan-ukb-us-east-1/ld_release/UKBB.{pop}.ldadj.bm")
    # matrix = bm.filter_rows(idx_list).filter_cols(idx_list).to_numpy()
    bm = bm.filter(idx_list, idx_list)

    # step 3: build variant DataFrame
    var_df        = ht.to_pandas()
    var_df        = var_df.sort_values("idx").reset_index(drop=True)
    var_df["chr"] = var_df["chr"].map(normalize_chr)
    var_df["a1"]  = var_df["a1"].astype(str).str.upper()
    var_df["a2"]  = var_df["a2"].astype(str).str.upper()
    var_df["snp"] = var_df["snp"].astype(str)
    # Re-index idx column to match sorted order
    var_df["idx"] = np.arange(len(var_df), dtype=int)

    # step 4: save
    np.save(matrix_path, bm.to_numpy())
    var_df.to_csv(variants_path, sep="\t", index=False, compression="gzip")
    print(f"[{pop}] Saved matrix   -> {matrix_path}  ({bm.shape[0]}x{bm.shape[1]})")
    print(f"[{pop}] Saved variants -> {variants_path}")
    return matrix_path, variants_path


# main

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Fetch Pan-UKBB EUR/AFR regional LD sub-matrices from S3 and save locally.\n"
            "Provide either (--chrom + --start + --end) or (--gene [--flank])."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    region = parser.add_argument_group("Region specification")
    region.add_argument("--chrom", help="Chromosome (e.g. 6 or chr6)")
    region.add_argument("--start", type=int, help="Region start bp (1-based, inclusive)")
    region.add_argument("--end",   type=int, help="Region end   bp (1-based, inclusive)")
    region.add_argument("--gene",  help="Gene name; coordinates looked up in the reference file")
    region.add_argument(
        "--flank", type=parse_flank_size, default=500_000,
        help="Flanking bp around gene (default: 500000). Supports 1Mb, 500kb, etc.",
    )

    parser.add_argument(
        "--pops", nargs="+", choices=["EUR", "AFR"], default=["EUR", "AFR"],
        help="Populations to fetch (default: EUR AFR)",
    )
    parser.add_argument(
        "--output-dir", required=True,
        help="Directory where .npy and .tsv.gz files will be written",
    )
    parser.add_argument(
        "--ref-genome", choices=["GRCh37", "GRCh38"], default="GRCh37",
        help="Reference genome version (default: GRCh37)",
    )
    parser.add_argument(
        "--gcp-project", default="xcheng-trial",
        help=(
            "GCP project ID for requester-pays GCS buckets. "
            "Required when the job runs outside Google Cloud."
        ),
    )
    args = parser.parse_args()

    # resolve region
    if args.gene:
        gene_ref = b37_GENE_REF if args.ref_genome == "GRCh37" else b38_GENE_REF
        chrom, g_start, g_end = get_gene_coordinates(args.gene, gene_ref)
        start = max(1, g_start - args.flank)
        end   = g_end + args.flank
        print(
            f"Gene {args.gene} -> chr{chrom}:{g_start}-{g_end}; "
            f"region with +/- {args.flank:,} bp flank: {start}-{end}"
        )
    elif args.chrom and args.start and args.end:
        chrom = normalize_chr(args.chrom)
        start, end = args.start, args.end
    else:
        parser.error("Provide either --gene [--flank] or all of --chrom, --start, --end")

    if start < 1:
        parser.error("--start must be >= 1")
    if end < start:
        parser.error("--end must be >= --start")

    os.makedirs(args.output_dir, exist_ok=True)

    # initialise Hail (single call)
    hail_tmp = os.path.join(args.output_dir, "hail_tmp")
    os.makedirs(hail_tmp, exist_ok=True)

    spark_conf = {
        "spark.driver.host":        "127.0.0.1",
        "spark.driver.bindAddress": "127.0.0.1",
        "spark.ui.enabled":         "false",
    }
    if args.gcp_project:
        spark_conf["spark.hadoop.fs.gs.requester.pays.mode"]       = "AUTO"
        spark_conf["spark.hadoop.fs.gs.requester.pays.project.id"] = args.gcp_project

    hl.init(
        default_reference=args.ref_genome,
        master="local[*]",
        app_name="fetch_ukbb_ld",
        quiet=True,
        log=os.path.join(args.output_dir, "fetch_ukbb_ld.hail.log"),
        tmp_dir=hail_tmp,
        spark_conf=spark_conf,
    )

    # fetch each population
    t0 = time.time()
    results = {}
    try:
        for pop in args.pops:
            results[pop] = fetch_pop_ld(pop, chrom, start, end, args.output_dir)
    finally:
        hl.stop()

    print(f"\nDone in {time.time() - t0:.1f}s")
    print("Files ready for use with step2_prepare_matrices_n_diagnostics.py --local-ld-dir:")
    for pop, (mat_p, var_p) in results.items():
        print(f"  {pop}: {mat_p}")
        print(f"       {var_p}")


if __name__ == "__main__":
    main()
