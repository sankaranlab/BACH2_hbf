#!/usr/bin/env python3
"""
Process GCTA panels for different ancestry groups.
Converts genome-wide genotype data to PLINK format for each chromosome and ancestry.

Usage:
    python process_gcta_panels.py --ancestry EUR
    python process_gcta_panels.py --ancestry AFR
    python process_gcta_panels.py --ancestry EAS
    python process_gcta_panels.py --ancestry SAS
"""

import argparse
import os
import sys
import pandas as pd
import hail as hl


def check_pyspark():
    """Verify PySpark installation."""
    try:
        import pyspark
    except ModuleNotFoundError:
        print("!" * 100 + "\n\n"
              "In the Researcher Workbench, Hail can only be used on a Dataproc cluster.\n"
              "Please use the 'Cloud Analysis Environment' side panel to update your runtime compute type.\n\n" +
              "!" * 100)
        sys.exit(1)


def load_datasets():
    """Load ancestry, QC, and relatedness data."""
    print("Loading ancestry predictions...")
    anc_pred = pd.read_csv(
        "gs://fc-aou-datasets-controlled/v5/wgs/vcf/aux/ancestry/ancestry_preds.tsv",
        sep="\t"
    )
    
    print("Loading QC flagged samples...")
    qc_flag = pd.read_csv(
        "gs://fc-aou-datasets-controlled/v5/wgs/vcf/aux/qc/flagged_samples.tsv",
        sep="\t"
    )
    
    print("Loading relatedness information...")
    relate_remove = pd.read_csv(
        "gs://fc-aou-datasets-controlled/v5/wgs/vcf/aux/relatedness/relatedness_flagged_samples.tsv",
        sep="\t"
    )
    relate_remove = relate_remove.rename(
        columns={"Unnamed: 0": "nonsense", "sample_id.s": "s"}
    )
    
    return anc_pred, qc_flag, relate_remove


def get_ancestry_samples(anc_pred, ancestry):
    """Extract sample IDs for a specific ancestry."""
    ancestry_lower = ancestry.lower()
    ancestry_samples = anc_pred[anc_pred.ancestry_pred == ancestry_lower]['research_id']
    return set(map(str, ancestry_samples))


def get_qc_samples(qc_flag):
    """Extract QC-flagged samples to remove."""
    return set(map(str, qc_flag['s']))


def get_related_samples(relate_remove):
    """Extract related samples to remove."""
    return set(map(str, relate_remove['s']))


def filter_matrix_table(mt, anc_pred, ancestry, qc_flag, relate_remove):
    """Apply all filtering steps to matrix table."""
    # Remove QC-flagged samples
    qc_samples = get_qc_samples(qc_flag)
    qc_set = hl.literal(qc_samples)
    mt = mt.filter_cols(~qc_set.contains(mt['s']))
    
    # Remove related samples
    related_samples = get_related_samples(relate_remove)
    related_set = hl.literal(related_samples)
    mt = mt.filter_cols(~related_set.contains(mt['s']))
    
    # Keep only specific ancestry
    ancestry_samples = get_ancestry_samples(anc_pred, ancestry)
    ancestry_set = hl.literal(ancestry_samples)
    mt = mt.filter_cols(ancestry_set.contains(mt['s']))
    
    # Filter sites marked as bad
    mt = mt.filter_rows(hl.is_missing(mt.filters))
    
    # Split multi-allelic variants
    bi = mt.filter_rows(hl.len(mt.alleles) == 2)
    bi = bi.annotate_rows(a_index=1, was_split=False)
    multi = mt.filter_rows(hl.len(mt.alleles) > 2)
    split = hl.split_multi_hts(
        multi,
        keep_star=False,
        left_aligned=False,
        vep_root='vep',
        permit_shuffle=False
    )
    mt = split.union_rows(bi)
    
    # Variant QC filtering
    mt = hl.variant_qc(mt)
    mt = mt.filter_rows(
        (mt.variant_qc.p_value_hwe > 1e-15) &
        (mt.variant_qc.AF[1] > 0.001) &
        (mt.variant_qc.call_rate > 0.90)
    )
    
    return mt


def process_ancestry(ancestry):
    """Process all chromosomes for a specific ancestry."""
    # Validate ancestry
    valid_ancestries = ['eur', 'afr', 'eas', 'sas']
    if ancestry.lower() not in valid_ancestries:
        raise ValueError(
            f"Invalid ancestry: {ancestry}. Must be one of {valid_ancestries}"
        )
    
    # Get environment variables
    bucket = os.getenv("WORKSPACE_BUCKET")
    wgs_path = os.getenv("WGS_HAIL_STORAGE_PATH")
    
    if not bucket or not wgs_path:
        raise RuntimeError(
            "Missing required environment variables:\n"
            "  WORKSPACE_BUCKET: {}\n"
            "  WGS_HAIL_STORAGE_PATH: {}".format(bucket, wgs_path)
        )
    
    print(f"\nProcessing {ancestry.upper()} panel")
    print(f"  Bucket: {bucket}")
    print(f"  WGS path: {wgs_path}")
    
    # Load shared datasets
    anc_pred, qc_flag, relate_remove = load_datasets()
    print(f"\nDataset summary:")
    print(f"  Ancestry counts:\n{anc_pred.ancestry_pred.value_counts()}")
    print(f"  QC-flagged samples: {len(qc_flag)}")
    print(f"  Related samples to remove: {len(relate_remove)}")
    
    # Process each chromosome
    output_prefix = f"{ancestry.lower()}-plink-AF1pct"
    
    for chrom in range(1, 23):
        print(f"\n[Chr{chrom}] Processing chromosome {chrom}...")
        
        # Read and filter matrix table
        mt = hl.read_matrix_table(wgs_path)
        mt = hl.filter_intervals(
            mt,
            [hl.parse_locus_interval(f"chr{chrom}")]
        )
        
        print(f"[Chr{chrom}] Filtering samples and variants...")
        mt = filter_matrix_table(mt, anc_pred, ancestry, qc_flag, relate_remove)
        
        # Count samples and variants after filtering
        n_samples = mt.count_cols()
        n_variants = mt.count_rows()
        print(f"[Chr{chrom}] After filtering: {n_variants} variants, {n_samples} samples")
        
        # Export to PLINK format
        output_path = f"{bucket}/{output_prefix}/chr{chrom}_filt{ancestry.upper()}_p1pct_gcta"
        print(f"[Chr{chrom}] Exporting to {output_path}...")
        hl.export_plink(mt, output_path)
        print(f"[Chr{chrom}] Complete")


def main():
    parser = argparse.ArgumentParser(
        description="Process GCTA panels for different ancestry groups",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python process_gcta_panels.py --ancestry EUR
  python process_gcta_panels.py --ancestry AFR
  python process_gcta_panels.py --ancestry EAS
  python process_gcta_panels.py --ancestry SAS
        """
    )
    
    parser.add_argument(
        '--ancestry',
        required=True,
        choices=['EUR', 'AFR', 'EAS', 'SAS'],
        help='Ancestry population to process (case-insensitive)'
    )
    
    parser.add_argument(
        '--chromosomes',
        type=str,
        default='1-22',
        help='Chromosomes to process (e.g., "1-22" or "1,2,3") [default: 1-22]'
    )
    
    args = parser.parse_args()
    
    # Check PySpark availability
    check_pyspark()
    
    # Initialize Hail
    print("Initializing Hail...")
    hl.init(default_reference='GRCh38')
    
    try:
        # Process the specified ancestry
        process_ancestry(args.ancestry)
        print("\n" + "="*80)
        print(f"Successfully completed processing for {args.ancestry.upper()}")
        print("="*80)
        
    except Exception as e:
        print(f"\nError: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
