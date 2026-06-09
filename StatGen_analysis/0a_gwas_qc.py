#!/usr/bin/env python3

import pandas as pd
import numpy as np
import argparse
import os
import time
import resource
import matplotlib.pyplot as plt

########################################
# Column mapping
########################################

COLUMN_MAP = {
    "snp": ["SNP"], #,"RSID","ID","VARIANT_ID"
    "id": ["ID"],
    "chr": ["CHR"], #,"CHROM"
    "pos": ["BP"], #,"POS"
    "a1": ["A1"],
    "a2": ["A2"],
    "beta": ["BETA","SIGNED_SUMSTAT"],
    "se": ["SE"],
    "p": ["P"], #,"PVALUE","LOG10P"
    "n": ["N"],
    "z": ["Z"],
    "ea_freq": ["AF_ALLELE2"],
    "maf": ["MAF"],
    "ea_count": ["AC_ALLELE2", "AC"],
    "mac": ["MAC"],
    "info": ["IMPUTATIONINFO","INFO"],
    "impute": ["IMPUTATION_gen_build"]
}

########################################
# Detect column names
########################################

def find_col(df, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    return None

def standardize_columns(df):

    rename = {}

    for std, syns in COLUMN_MAP.items():
        col = find_col(df, syns)
        if col:
            rename[col] = std

    df = df.rename(columns=rename)

    required = ["snp","chr","pos","a1","a2","beta","se","p","z","n"]
    missing = [c for c in required if c not in df.columns]

    if missing:
        raise ValueError(f"Missing columns {missing}")
    
    # if ea_freq or ea_count is present, convert them to maf and mac
    if "ea_freq" in df.columns and "n" in df.columns and "maf" not in df.columns:
        df["maf"] = np.where(df["ea_freq"] <= 0.5, df["ea_freq"], 1 - df["ea_freq"])
    if "ea_count" in df.columns and "n" in df.columns and "mac" not in df.columns:
        df["mac"] = np.where(df["ea_count"] <= 0.5 * df["n"], df["ea_count"], df["n"] - df["ea_count"])

    return df


###################################################
# Variant filtering
###################################################

def is_snp(a1,a2):
    bases = {"A","C","G","T"}
    return (a1 in bases) and (a2 in bases)


def sort_by_locus(df):
    df = df.copy()
    chr_num = pd.to_numeric(df["chr"], errors="coerce")
    pos_num = pd.to_numeric(df["pos"], errors="coerce")
    df = df.assign(_chr_num=chr_num, _pos_num=pos_num)
    df = df.sort_values(by=["_chr_num", "_pos_num", "chr", "pos"], na_position="last")
    return df.drop(columns=["_chr_num", "_pos_num"])


def compute_gc_lambda(z_values):
    """Compute genomic control lambda from Z values.
    
    lambda = median(Z^2) / 0.4549, where 0.4549 is qchisq(0.5, 1).
    """
    z_clean = z_values[np.isfinite(z_values)]
    if len(z_clean) == 0:
        return np.nan
    return np.median(z_clean**2) / 0.4549


def plot_z_diagnostics(df, outdir, name, inconsistent_mask):
    """Create MAC-stratified diagnostic plot of Z_wald vs Z_reported for inconsistent sites."""
    if "mac" not in df.columns or "beta" not in df.columns or "se" not in df.columns or "z" not in df.columns:
        print("Skipping Z diagnostics plot: missing MAC, beta, se, or z column")
        return
    
    z_wald = df["beta"] / df["se"]
    z_reported = df["z"]
    valid = np.isfinite(z_wald) & np.isfinite(z_reported) & np.isfinite(df["mac"])
    
    if valid.sum() == 0:
        print("Skipping Z diagnostics plot: no valid Z values")
        return
    
    # Subset to inconsistent sites for MAC stratification
    incons_valid = inconsistent_mask & valid
    if incons_valid.sum() == 0:
        print("Skipping Z diagnostics plot: no valid inconsistent sites")
        return
    
    # Create MAC strata based on inconsistent sites only
    mac_vals_incons = df.loc[incons_valid, "mac"]
    mac_quartiles = mac_vals_incons.quantile([0.25, 0.5, 0.75])
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle(f"Z Diagnostics (Wald vs Reported) - {name}\nMAC-stratified view of inconsistent sites")
    
    strata = [
        (mac_vals_incons < mac_quartiles[0.25], "MAC < Q1"),
        ((mac_vals_incons >= mac_quartiles[0.25]) & (mac_vals_incons < mac_quartiles[0.5]), "Q1 <= MAC < Q2"),
        ((mac_vals_incons >= mac_quartiles[0.5]) & (mac_vals_incons < mac_quartiles[0.75]), "Q2 <= MAC < Q3"),
        (mac_vals_incons >= mac_quartiles[0.75], "MAC >= Q3"),
    ]
    
    for ax, (mask, label) in zip(axes.flat, strata):
        subset_indices = incons_valid.values.copy()
        subset_indices[incons_valid.values] = mask.values
        subset_indices = np.where(subset_indices)[0]
        
        if len(subset_indices) > 0:
            z_w = z_wald.iloc[subset_indices]
            z_r = z_reported.iloc[subset_indices]
            ax.scatter(z_w, z_r, alpha=0.5, s=20, edgecolors="black", linewidth=0.5)
            
            # Add diagonal reference line
            all_z = np.concatenate([z_w, z_r])
            lim_min, lim_max = np.min(all_z), np.max(all_z)
            margin = (lim_max - lim_min) * 0.05
            lims = [lim_min - margin, lim_max + margin]
            ax.plot(lims, lims, "r--", alpha=0.75, zorder=0, linewidth=1, label="Identity")
            ax.set_xlim(lims)
            ax.set_ylim(lims)
            
            ax.set_xlabel("Z_wald (beta/se)", fontsize=10)
            ax.set_ylabel("Z_reported", fontsize=10)
            ax.set_title(f"{label} (n={len(subset_indices)})", fontsize=11)
            ax.grid(True, alpha=0.3)
            ax.legend()
        else:
            ax.text(0.5, 0.5, f"{label}\n(n=0)", ha="center", va="center", fontsize=12)
            ax.set_xticks([])
            ax.set_yticks([])
    
    plotfile = os.path.join(outdir, f"{name}.z_diagnostics.png")
    plt.savefig(plotfile, dpi=100, bbox_inches="tight")
    plt.close()
    print(f"Saved Z diagnostics plot: {plotfile}")


def qc_filter(df, pass_z_check=False, se_cap=None, min_maf=0.0001):

    discarded = []

    def add_discarded(rows, reason):
        if len(rows) == 0:
            return
        rows = rows.copy()
        rows["discard_reason"] = reason
        discarded.append(rows)

    print(f"Input variants: {len(df)}")

    df = df.copy()
    numeric_cols = ["beta", "se", "p", "z", "n", "maf", "mac", "info"]
    present_numeric_cols = []
    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
            present_numeric_cols.append(col)

    if present_numeric_cols:
        finite_mask = np.ones(len(df), dtype=bool)
        for col in present_numeric_cols:
            finite_mask &= np.isfinite(df[col])
        bad = df[~finite_mask]
        add_discarded(bad, "non_finite_numeric")
        df = df[finite_mask]
        print(f"Dropped rows with NA/Inf in numeric columns: {len(bad)}")

    # remove non SNP
    mask = df.apply(lambda r: is_snp(r.a1,r.a2),axis=1)
    bad = df[~mask]
    add_discarded(bad, "non_snp")
    df = df[mask]
    print(f"Dropped non-SNP sites: {len(bad)}")

    # remove multi allelic
    dup_mask = df.duplicated(subset=["chr", "pos"], keep="first")
    bad = df[dup_mask]
    add_discarded(bad, "duplicate_position")
    df = df[~dup_mask]
    print(f"Dropped multiallelic/duplicate-position sites: {len(bad)}")

    # sanity SE filters
    bad = df[(df.se<=0)]  # | (df.se<1e-6)
    add_discarded(bad, "invalid_se")
    df = df.drop(bad.index)
    print(f"Dropped SE-filtered sites: {len(bad)}")

    #sanity p-value filters
    bad = df[(df.p<0) | (df.p>1)]
    add_discarded(bad, "invalid_p")
    df = df.drop(bad.index)
    print(f"Dropped p-value-filtered sites: {len(bad)}")

    # If Z is one-sided (all >=0 or all <=0), treat it as unsigned and align to beta sign.
    z_num = pd.to_numeric(df["z"], errors="coerce")
    beta_num = pd.to_numeric(df["beta"], errors="coerce")
    valid_z = np.isfinite(z_num)
    if valid_z.any():
        z_valid_vals = z_num[valid_z]
        if (z_valid_vals >= 0).all() or (z_valid_vals <= 0).all():
            beta_sign = np.sign(beta_num)
            # For beta==0 or missing beta, keep positive sign to avoid introducing NaNs here.
            beta_sign = beta_sign.where(np.isfinite(beta_sign) & (beta_sign != 0), 1.0)
            df.loc[valid_z, "z"] = z_num[valid_z].abs() * beta_sign[valid_z]
            print("Detected one-sided Z values; aligned Z sign to beta before consistency check")

    # Z consistency diagnostics (no filtering, diagnostic only)
    z_inconsistent_mask = None
    if pass_z_check:
        print("Skipped Z consistency check (--pass-z-check set)")
    else:
        z_wald = df.beta / df.se
        z_reported = df["z"]
        z_inconsistent_mask = np.abs(z_wald - z_reported) > 0.1
        
        if z_inconsistent_mask.any():
            n_inconsistent = z_inconsistent_mask.sum()
            print(f"Found {n_inconsistent} Z-inconsistent sites (|Z_wald - Z_reported| > 0.1)")
            
            # Check for systematic bias
            z_diff = z_wald[z_inconsistent_mask] - z_reported[z_inconsistent_mask]
            print(f"Z bias statistics for inconsistent sites (Z_wald - Z_reported):")
            print(f"  mean: {z_diff.mean():.4g}")
            print(f"  median: {z_diff.median():.4g}")
            print(f"  std: {z_diff.std():.4g}")
            print(f"  min: {z_diff.min():.4g}")
            print(f"  max: {z_diff.max():.4g}")
        else:
            print("No Z-inconsistent sites detected")

    # p=0 fix
    # df.loc[df.p==0,"p"]=1e-300
    # report any remaining p=0 sites
    bad = df[df.p==0]
    add_discarded(bad, "p_zero")
    # df = df.drop(bad.index)
    # print(f"Dropped p=0 sites: {len(bad)}")

    # extreme Z
    bad = df[np.abs(df.z)>80]
    add_discarded(bad, "extreme_z")
    df = df.drop(bad.index)
    print(f"Dropped extreme-Z sites: {len(bad)}")

    # maf
    if "maf" in df.columns:
        bad = df[df.maf < min_maf]
        add_discarded(bad, "low_maf")
        df = df.drop(bad.index)
        print(f"Dropped low-MAF sites (<{min_maf}): {len(bad)}")
    
    # mac
    if "mac" in df.columns and "n" in df.columns:
        # Compute per-row MAC cutoff using each row's N value
        mac_cutoff_per_row = np.maximum(
            np.maximum(min_maf * 2 * df["n"], np.minimum(40.0, np.sqrt(df["n"]))),
            20.0,
        )
        print(
            f"Applied per-row MAC cutoff: max(min_maf*2*N, min(40, sqrt(N)), 20) for each row\n"
            f"  min_maf={min_maf}\n"
            f"  MAC cutoff range: {mac_cutoff_per_row.min():.4g} to {mac_cutoff_per_row.max():.4g}"
        )

        bad = df[df["mac"] < mac_cutoff_per_row]
        add_discarded(bad, "low_mac")
        df = df.drop(bad.index)
        print(f"Dropped low-MAC sites (per-row cutoff): {len(bad)}")


    # info
    if "info" in df.columns:
        bad = df[df["info"] < 0.6]
        add_discarded(bad, "low_info")
        df = df.drop(bad.index)
        print(f"Dropped low-INFO sites (<0.6): {len(bad)}")

    # # ambiguous SNP
    # amb = {("A","T"),("T","A"),("C","G"),("G","C")}
    # mask = df.apply(lambda r:(r.a1,r.a2) in amb,axis=1)
    # bad = df[mask]
    # discarded.append(bad)
    # df = df[~mask]
    # print(f"Dropped ambiguous strand sites: {len(bad)}")

    # report the distribution of se and N for the retained variants
    for col in ("se", "n"):
        if col in df.columns:
            vals = df[col].dropna()
            if len(vals) == 0:
                print(f"{col.upper()} distribution: no non-missing values among retained variants")
                continue
            pct = vals.quantile([0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99])
            print(
                f"{col.upper()} distribution (n={len(vals)}):\n"
                f"\tmin={vals.min():.4g} "
                f"p1={pct[0.01]:.4g} p5={pct[0.05]:.4g} "
                f"p25={pct[0.25]:.4g} median={pct[0.50]:.4g} "
                f"p75={pct[0.75]:.4g} p95={pct[0.95]:.4g} "
                f"p99={pct[0.99]:.4g} max={vals.max():.4g}"
            )
    if se_cap is not None:
        print(f"Applied SE cap filter: SE <= {se_cap}")
        bad = df[df.se > se_cap]
        add_discarded(bad, "se_cap")
        df = df.drop(bad.index)
        print(f"Dropped SE-cap-filtered sites: {len(bad)}")
        # report SE distribution after cap filter
        vals = df["se"].dropna()
        if len(vals) > 0:
            pct = vals.quantile([0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99])
            print(
                f"SE distribution after SE cap filter (n={len(vals)}):\n "
                f"\tmin={vals.min():.4g} "
                f"p1={pct[0.01]:.4g} p5={pct[0.05]:.4g} "
                f"p25={pct[0.25]:.4g} median={pct[0.50]:.4g} "
                f"p75={pct[0.75]:.4g} p95={pct[0.95]:.4g} "
                f"p99={pct[0.99]:.4g} max={vals.max():.4g}"
            )
    # if max se > 5, print a warning about potential SE issues
    # and print the first few variants with se > 5 for inspection
    if "se" in df.columns:
        max_se = df["se"].max()
        if max_se > 5:
            print(f"Warning: Maximum SE > 5 ({max_se:.4g}). Potential SE issues detected.")
            n_large_se = (df["se"] > 5).sum()
            print(f"Number of variants with SE > 5: {n_large_se}")
            print("First few variants with SE > 5:")
            print(df[df["se"] > 5].head(min(10, n_large_se)).to_string(index=False))

    if discarded:
        discarded = pd.concat(discarded, ignore_index=True)
    else:
        discarded = pd.DataFrame(columns=list(df.columns) + ["discard_reason"])

    print(f"Retained variants after QC: {len(df)}")

    return df, discarded, z_inconsistent_mask

###################################################
# Load + QC
###################################################

def process_file(file, outdir, suffix="", pass_z_check=False, se_cap=None, min_maf=0.0001):

    name = os.path.basename(file).split(".")[0]

    df = pd.read_csv(file, sep=None, engine="python")

    df = standardize_columns(df)

    # Report genomic control lambda before QC
    print("\n=== Genomic Control Lambda Before QC ===")
    if "z" in df.columns:
        z_vals_before = pd.to_numeric(df["z"], errors="coerce")
        gc_lambda_before = compute_gc_lambda(z_vals_before)
        print(f"Genomic control lambda (before QC): {gc_lambda_before:.6f}\n")
    else:
        print("Z column not found for pre-QC lambda calculation")
        gc_lambda_before = np.nan

    clean, discarded, z_inconsistent_mask = qc_filter(
        df,
        pass_z_check=pass_z_check,
        se_cap=se_cap,
        min_maf=min_maf,
    )

    # Report genomic control lambda after QC
    print("\n=== Genomic Control Lambda After QC ===")
    if "z" in clean.columns:
        z_vals_after = pd.to_numeric(clean["z"], errors="coerce")
        gc_lambda_after = compute_gc_lambda(z_vals_after)
        print(f"Genomic control lambda (after QC): {gc_lambda_after:.6f}\n")
        if not np.isnan(gc_lambda_before):
            ratio = gc_lambda_after / gc_lambda_before
            print(f"Lambda ratio (after/before): {ratio:.4f}")
    else:
        print("Z column not found for post-QC lambda calculation")
        gc_lambda_after = np.nan
    
    # Plot Z diagnostics if inconsistent sites were found in original data
    if z_inconsistent_mask is not None and z_inconsistent_mask.any():
        print(f"\nGenerating MAC-stratified Z diagnostics plot...")
        plot_z_diagnostics(df, outdir, name, z_inconsistent_mask)

    clean = sort_by_locus(clean)
    discarded = sort_by_locus(discarded)

    cleanfile = os.path.join(outdir, name + suffix + ".clean.tsv.gz")
    discardfile = os.path.join(outdir, name + suffix + ".discarded.tsv.gz")

    clean.to_csv(cleanfile, sep="\t", index=False, compression="gzip")
    discarded.to_csv(discardfile, sep="\t", index=False, compression="gzip")

    return cleanfile


########################################
# Main
########################################

def main():

    start_time = time.time()

    parser = argparse.ArgumentParser()

    parser.add_argument("--input", "-i", required=True, 
                        help="Single summary stats file (one per SLURM array task)")
    parser.add_argument("--outdir", "-o", required=True, 
                        help="Output directory")
    parser.add_argument("--suffix", default="",
                        help="Suffix to add to output files (before .clean.tsv.gz)")
    parser.add_argument("--pass-z-check", action="store_true",
                        help="Bypass Z consistency check")
    parser.add_argument("--min-maf", type=float, default=0.0001,
                        help="MAF cutoff for filtering (default: 0.0001)")
    parser.add_argument("--se-cap", type=float, default=None,
                        help="If not None, also filter out rows with SE > se_cap")

    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    process_file(args.input, args.outdir, suffix=args.suffix,
                 pass_z_check=args.pass_z_check, se_cap=args.se_cap,
                 min_maf=args.min_maf)

    elapsed = time.time() - start_time

    print(f"Total runtime (wall): {elapsed:.2f} seconds")


if __name__ == "__main__":
    main()