#!/usr/bin/env python3
import argparse
import pandas as pd
import numpy as np
import os
import glob
import gzip
import time
import re
from pathlib import Path
import matplotlib.pyplot as plt
from scipy import stats

# set global parameters
ROOT = "/lab-share/Hem-Sankaran-e2/Public"
b37_GENE_REF = os.path.join(ROOT, "ref_genomes/human/hg19_ncbiRefSeq_RefSeqAll_2026UCSC.tsv.gz")
b38_GENE_REF = os.path.join(ROOT, "ref_genomes/human/hg38_ncbiRefSeq_AllGenes_2026UCSC.tsv.gz")
THAI_LD_PATH = os.path.join(ROOT, "projects/xhcheng/HbF/bach2_multi_finemap/ld/thai_QCed_BACH2_1Mb_flanking.ld")


def get_gene_coordinates(gene_name, gene_ref_path=b37_GENE_REF):
    """
    Extract gene coordinates from reference genome file.
    ref header: #bin	name	chrom	strand	txStart	txEnd	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	score	name2	cdsStartStat	cdsEndStat	exonFrames
    Returns: (chrom, start, end) tuple with the widest range for the gene.
    """
    print(f"Looking up coordinates for gene: {gene_name}")
    
    if not os.path.exists(gene_ref_path):
        raise FileNotFoundError(f"Gene reference file not found: {gene_ref_path}")
    
    open_func = gzip.open if gene_ref_path.endswith('.gz') else open
    
    gene_records = []
    with open_func(gene_ref_path, 'rt') as f:
        header = f.readline().strip().split('\t')
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 13:
                continue
            
            # Header format: #bin name chrom strand txStart txEnd cdsStart cdsEnd exonCount exonStarts exonEnds score name2 ...
            record_gene_name = fields[12]
            
            if record_gene_name == gene_name:
                chrom = fields[2].replace('chr', '')  # Remove 'chr' prefix if present
                # only consider assembled chromosomes (not contigs)
                if chrom not in [str(i) for i in range(1, 23)] + ['X', 'Y']:
                    continue
                start = int(fields[4])  # txStart
                end = int(fields[5])    # txEnd
                gene_records.append((chrom, start, end))
    
    if not gene_records:
        raise ValueError(f"Gene '{gene_name}' not found in reference file")
    print(f"Found {len(gene_records)} records for gene {gene_name}")
    print(gene_records)

    # Get the widest range
    chrom = gene_records[0][0]
    min_start = min(rec[1] for rec in gene_records)
    max_end = max(rec[2] for rec in gene_records)
    
    print(f"Found gene {gene_name} on chr{chrom}: {min_start}-{max_end}")
    return chrom, min_start, max_end


def parse_flank_size(value):
    """
    Parse flank size argument that can be in various formats.
    
    Supported formats:
        - Plain integers: "500000"
        - Scientific notation: "1e6", "5e5"
        - With unit suffixes: "1Mb", "500kb", "2m", "1000b"
    
    Args:
        value: String input for flank size
    
    Returns:
        Integer flank size in base pairs
    
    Raises:
        argparse.ArgumentTypeError: If value cannot be parsed
    """
    value = value.strip()
    
    # Unit multipliers
    units = {
        'b': 1,
        'kb': 1e3,
        'mb': 1e6,
        'gb': 1e9,
        'k': 1e3,
        'm': 1e6,
        'g': 1e9,
    }
    
    try:
        # Try parsing as plain number (handles both integers and scientific notation)
        return int(float(value))
    except ValueError:
        pass
    
    # Try parsing with unit suffix
    value_lower = value.lower()
    for unit, multiplier in units.items():
        if value_lower.endswith(unit):
            try:
                numeric_part = value_lower[:-len(unit)].strip()
                numeric_value = float(numeric_part)
                return int(numeric_value * multiplier)
            except ValueError:
                pass
    
    # If we get here, couldn't parse the value
    raise argparse.ArgumentTypeError(
        f"Cannot parse flank size '{value}'. "
        f"Use formats like: 500000, 1e6, 1Mb, 500kb, 2m"
    )



def normalize_chr(value):
    c = str(value).strip().replace("chr", "").replace("CHR", "")
    # Pandas may read integer chromosome columns as floats (e.g. "1.0"); strip the .0.
    if c.endswith(".0") and c[:-2].lstrip("-").isdigit():
        c = c[:-2]
    if c == "M":
        c = "MT"
    return c


def chr_sort_key(chrom):
    c = normalize_chr(chrom)
    if c.isdigit():
        return int(c)
    if c == "X":
        return 23
    if c == "Y":
        return 24
    if c == "MT":
        return 25
    return 999


def read_table(path):
    compression = "gzip" if path.endswith(".gz") else None
    return pd.read_csv(path, sep="\t", compression=compression, low_memory=False)


def parse_input_files(items):
    files = {}
    for token in items:
        if "," not in token:
            raise ValueError(f"Invalid --input-files token '{token}'. Expected 'path,ANC'.")
        path, anc = token.rsplit(",", 1)
        anc = anc.strip().upper()
        path = path.strip()
        if anc in files:
            raise ValueError(f"Duplicate ancestry '{anc}' in --input-files.")
        files[anc] = path
    required = {"EUR", "AFR", "THAI"}
    missing = required - set(files)
    if missing:
        raise ValueError(f"Missing required ancestries: {sorted(missing)}")
    return files


def standardize_gwas(df, pop):
    """
    Standardize GWAS summary statistics for a given population.

    Args:
        df: DataFrame containing GWAS summary statistics
        pop: Population identifier ("EUR", "AFR", "THAI")

    Returns:
        DataFrame with standardized columns
    """
    pop = pop.upper()
    input_n = len(df)
    if pop in {"EUR", "AFR"}:
        need = ["CHR", "GENPOS", "MarkerName", "Allele1", "Allele2", 
                "Freq1", "FreqSE", "Effect", "StdErr", "P-value",
                "HetISq", "HetChiSq", "HetPVal","Direction","TotalN"]
        miss = [c for c in need if c not in df.columns]
        if miss:
            raise ValueError(f"{pop} missing required columns: {miss}")
        out = pd.DataFrame({
            "snp": df["MarkerName"].astype(str),
            "chr": df["CHR"].astype(str).map(normalize_chr),
            "pos": pd.to_numeric(df["GENPOS"], errors="coerce"),
            "a1": df["Allele1"].astype(str).str.upper(),
            "a2": df["Allele2"].astype(str).str.upper(),
            "maf": pd.to_numeric(df["Freq1"], errors="coerce"),
            "beta": pd.to_numeric(df["Effect"], errors="coerce"),
            "se": pd.to_numeric(df["StdErr"], errors="coerce"),
            "p": pd.to_numeric(df["P-value"], errors="coerce"),
            "n": pd.to_numeric(df["TotalN"], errors="coerce"),
            "I2": pd.to_numeric(df["HetISq"], errors="coerce") if "HetISq" in df.columns else np.nan,
            "ChiSq": pd.to_numeric(df["HetChiSq"], errors="coerce") if "HetChiSq" in df.columns else np.nan,
            "HetPVal": pd.to_numeric(df["HetPVal"], errors="coerce") if "HetPVal" in df.columns else np.nan,
            "direction": df["Direction"].astype(str) if "Direction" in df.columns else np.nan,
        })
        # only keep variants with at most one missing cohort (i.e. direction has at most one "?")
        if "direction" in out.columns:
            out = out[out["direction"].apply(lambda s: s.count("?") <= 1)].copy()
    elif pop == "THAI":
        need = ["snp", "hg19_chr", "hg19_pos", "a1", "a2", "maf", "beta", "p", "n", "se"]
        miss = [c for c in need if c not in df.columns]
        if miss:
            raise ValueError(f"THAI missing required columns: {miss}")
        out = pd.DataFrame({
            "snp": df["snp"].astype(str),
            "chr": df["hg19_chr"].astype(str).map(normalize_chr),
            "pos": pd.to_numeric(df["hg19_pos"], errors="coerce"),
            "a1": df["a1"].astype(str).str.upper(),
            "a2": df["a2"].astype(str).str.upper(),
            "maf": pd.to_numeric(df["maf"], errors="coerce"),
            "beta": pd.to_numeric(df["beta"], errors="coerce"),
            "se": pd.to_numeric(df["se"], errors="coerce"),
            "p": pd.to_numeric(df["p"], errors="coerce"),
            "n": pd.to_numeric(df["n"], errors="coerce"),
            "I2": np.nan,
            "HetPVal": np.nan,
        })
    else:
        raise ValueError(f"Unsupported population: {pop}")

    before = len(out)
    out = out.dropna(subset=["snp", "chr", "pos", "a1", "a2", "beta", "se", "p"])
    print(f"[{pop}] standardize dropna required fields: {before} -> {len(out)} (dropped {before - len(out)})")

    before = len(out)
    out = out[out["se"] > 0].copy()
    print(f"[{pop}] standardize se>0 filter: {before} -> {len(out)} (dropped {before - len(out)})")

    before = len(out)
    out = out[out["chr"].map(chr_sort_key) <= 25].copy()
    print(f"[{pop}] standardize autosome/sex chr filter: {before} -> {len(out)} (dropped {before - len(out)})")

    out["maf"] = out["maf"].where((out["maf"] >= 0) & (out["maf"] <= 1), np.nan)
    out = out.sort_values(["chr", "pos"], key=lambda s: s.map(chr_sort_key) if s.name == "chr" else s)
    before = len(out)
    out = out.drop_duplicates(subset=["snp"], keep="first")
    print(f"[{pop}] standardize duplicate SNP filter: {before} -> {len(out)} (dropped {before - len(out)})")
    print(f"[{pop}] standardize summary: input {input_n} -> final {len(out)}")
    return out


def harmonize_three_way(gwas_by_pop):
    ref = gwas_by_pop["EUR"][["snp", "chr", "pos", "a1", "a2", "beta", "se", "p", "maf", "n"]].copy()
    ref = ref.rename(columns={
        "a1": "a1_ref", "a2": "a2_ref", "beta": "beta_EUR", "se": "se_EUR", "p": "p_EUR", "maf": "maf_EUR", "n": "n_EUR"
    })

    out = ref
    print(f"[harmonize] EUR baseline SNPs: {len(out)}")
    for pop in ["AFR", "THAI"]:
        d = gwas_by_pop[pop][["snp", "a1", "a2", "beta", "se", "p", "maf", "n"]].copy()
        m = out.merge(d, on="snp", how="inner", suffixes=("", f"_{pop}"))
        print(f"[harmonize] merge with {pop}: {len(out)} EUR-ref SNPs + {len(d)} {pop} SNPs -> {len(m)} shared SNP IDs")

        same = (m["a1"] == m["a1_ref"]) & (m["a2"] == m["a2_ref"])
        flip = (m["a1"] == m["a2_ref"]) & (m["a2"] == m["a1_ref"])
        same_n = int(same.sum())
        flip_n = int(flip.sum())
        mismatch_n = int(len(m) - same_n - flip_n)
        print(f"[harmonize] {pop} allele status among shared IDs: same={same_n}, flip={flip_n}, mismatch_drop={mismatch_n}")
        m = m[same | flip].copy()

        m[f"beta_{pop}"] = np.where(flip.loc[m.index], -m["beta"], m["beta"])
        m[f"maf_{pop}"] = np.where(flip.loc[m.index], 1.0 - m["maf"], m["maf"])
        m[f"se_{pop}"] = m["se"]
        m[f"p_{pop}"] = m["p"]
        m[f"n_{pop}"] = m["n"]
        out = m.drop(columns=["a1", "a2", "beta", "se", "p", "maf", "n"])
        print(f"[harmonize] retained after harmonizing with {pop}: {len(out)} SNPs")

    out = out.rename(columns={"a1_ref": "a1", "a2_ref": "a2"})
    print(f"[harmonize] final 3-way harmonized SNPs: {len(out)}")
    return out


def compute_heterogeneity(df):
    b = df[["beta_EUR", "beta_AFR", "beta_THAI"]].to_numpy(dtype=float)
    se = df[["se_EUR", "se_AFR", "se_THAI"]].to_numpy(dtype=float)
    w = 1.0 / np.square(se)
    sum_w = np.sum(w, axis=1)
    beta_meta = np.sum(w * b, axis=1) / sum_w
    q = np.sum(w * np.square(b - beta_meta[:, None]), axis=1)
    k = np.sum(np.isfinite(b) & np.isfinite(se), axis=1)
    df_het = np.maximum(k - 1, 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        i2 = np.where(q > 0, np.maximum(0.0, (q - (k - 1)) / q) * 100.0, 0.0)
        c = sum_w - (np.sum(np.square(w), axis=1) / sum_w)
        tau2 = np.where(c > 0, np.maximum(0.0, (q - (k - 1)) / c), 0.0)
    het_p = stats.chi2.sf(q, df=np.maximum(k - 1, 1))

    out = df.copy()
    out["beta_meta"] = beta_meta
    out["Q"] = q
    out["df_het"] = df_het
    out["I2"] = i2
    out["tau2"] = tau2
    out["HetP"] = het_p
    return out


def global_diagnostics(df_harmonized, raw_by_pop, output_prefix):
    """
    Generate comprehensive global diagnostics plots.
    
    Args:
        df_harmonized: Harmonized multi-ancestry DataFrame
        raw_by_pop: Dict of raw DataFrames by population (EUR, AFR, THAI)
        output_prefix: Output prefix for plot file
    """
    # Validate inputs
    if df_harmonized.empty:
        print("Warning: harmonized DataFrame is empty; skipping global diagnostics")
        return
    
    required_pops = {"EUR", "AFR", "THAI"}
    missing_pops = required_pops - set(raw_by_pop.keys())
    if missing_pops:
        print(f"Warning: missing populations {missing_pops}; skipping global diagnostics")
        return
    
    fig, ax = plt.subplots(4, 3, figsize=(18, 14))

    # Panel A, B, C: Sign-adjusted beta vs MAF for each population
    for i, (pop, color) in enumerate([("EUR", "tab:blue"), ("AFR", "tab:green"), ("THAI", "tab:brown")]):
        if f"maf_{pop}" not in df_harmonized.columns or f"beta_{pop}" not in df_harmonized.columns:
            ax[0, i].text(0.5, 0.5, f"Missing {pop} columns", ha="center", va="center", transform=ax[0, i].transAxes)
            continue
        x = df_harmonized[f"maf_{pop}"]
        y = np.sign(df_harmonized["beta_EUR"]) * df_harmonized[f"beta_{pop}"]
        s = np.clip(-np.log10(np.maximum(df_harmonized[f"p_{pop}"], 1e-300)), 0, 60) * 4 + 6
        ax[0, i].scatter(x, y, s=s, alpha=0.45, label=pop, color=color)
        ax[0, i].set_xlabel("MAF")
        ax[0, i].set_ylabel("Sign-adjusted beta")
        ax[0, i].set_title(f"{pop}: Sign-adjusted Beta vs MAF")
        ax[0, i].legend()

    # Panel D, E, F: Pairwise beta correlations
    pairs = [("AFR", "EUR"), ("THAI", "EUR"), ("THAI", "AFR")]
    for (p1, p2), a in zip(pairs, [ax[1, 0], ax[1, 1], ax[1, 2]]):
        required_cols = [f"beta_{p1}", f"beta_{p2}"]
        if not all(c in df_harmonized.columns for c in required_cols):
            a.text(0.5, 0.5, f"Missing {p1}/{p2} columns", ha="center", va="center", transform=a.transAxes)
            continue
        x = df_harmonized[f"beta_{p1}"]
        y = df_harmonized[f"beta_{p2}"]
        # Skip if insufficient variance
        if x.std() == 0 or y.std() == 0:
            a.text(0.5, 0.5, f"No variance in {p1}/{p2}", ha="center", va="center", transform=a.transAxes)
            continue
        r = np.corrcoef(x, y)[0, 1]
        a.scatter(x, y, s=12, alpha=0.45)
        # plot regression line
        m, b = np.polyfit(x, y, 1)
        a.plot(x, m * x + b, color="red", lw=1)
        a.text(0.05, 0.95, f"r={r:.3f}", transform=a.transAxes, ha="left", va="top", color="red")
        a.set_xlabel(f"beta {p1}")
        a.set_ylabel(f"beta {p2}")
        a.set_title(f"{p1} vs {p2} (r={r:.3f})")
    
    # Panel G: I^2 distribution by population
    has_i2_data = all(pop in raw_by_pop and "I2" in raw_by_pop[pop].columns for pop in ["EUR", "AFR"])
    if has_i2_data:
        eur = raw_by_pop["EUR"]
        afr = raw_by_pop["AFR"]
        eur_i2 = eur["I2"].dropna()
        afr_i2 = afr["I2"].dropna()
        
        if len(eur_i2) > 0 and len(afr_i2) > 0:
            ax[2, 0].hist([eur_i2, afr_i2], bins=np.arange(0, 110, 10), density=True, alpha=0.5, 
                          label=["EUR", "AFR"], color=["tab:blue", "tab:green"])
            ax[2, 0].set_xlabel(r"$I^2$")
            ax[2, 0].set_ylabel("Density")
            # Annotate statistics
            for pop, color in [("EUR", "tab:blue"), ("AFR", "tab:green")]:
                d = raw_by_pop[pop]
                i2_vals = d["I2"].dropna()
                if len(i2_vals) > 0:
                    median_i2 = np.median(i2_vals)
                    prop_i2_gt_50 = (i2_vals > 50).mean()
                    y_pos = 0.95 - 0.1 * list(sorted(raw_by_pop.keys())).index(pop)
                    ax[2, 0].text(0.95, y_pos, f"{pop}: median={median_i2:.1f}, >50%={prop_i2_gt_50:.1%}",
                                transform=ax[2, 0].transAxes, ha="right", va="top", color=color, fontsize=9)
            ax[2, 0].set_title(r"$I^2$ distribution")
            ax[2, 0].legend()
        else:
            ax[2, 0].text(0.5, 0.5, "Insufficient I² data", ha="center", va="center", transform=ax[2, 0].transAxes)
    else:
        ax[2, 0].text(0.5, 0.5, "Missing I² columns", ha="center", va="center", transform=ax[2, 0].transAxes)

    # Panel H: I^2 vs association strength
    if has_i2_data:
        for pop, color in [("EUR", "tab:blue"), ("AFR", "tab:green")]:
            d = raw_by_pop[pop]
            if "I2" in d.columns and "p" in d.columns:
                valid = ~(d["I2"].isna() | d["p"].isna())
                if valid.sum() > 0:
                    ax[2, 1].scatter(d.loc[valid, "I2"], -np.log10(np.maximum(d.loc[valid, "p"], 1e-300)), 
                                   s=10, alpha=0.4, label=pop, color=color)
        ax[2, 1].set_xlabel(r"$I^2$")
        ax[2, 1].set_ylabel(r"-\log_{10}(p)")
        ax[2, 1].set_title(r"$I^2$ vs association strength")
        if len(ax[2, 1].lines) > 0 or len(ax[2, 1].collections) > 0:
            ax[2, 1].legend()
    else:
        ax[2, 1].text(0.5, 0.5, "Missing I² data", ha="center", va="center", transform=ax[2, 1].transAxes)

    # Panel I: Heterogeneity p-value QQ plot
    qq_ax = ax[2, 2]
    qq_lim = 0.0
    has_points = False
    for pop, color in [("EUR", "tab:blue"), ("AFR", "tab:green")]:
        if pop not in raw_by_pop or "HetPVal" not in raw_by_pop[pop].columns:
            continue
        pvals = raw_by_pop[pop]["HetPVal"].dropna()
        pvals = pvals[(pvals > 0) & (pvals <= 1)]
        if len(pvals) >= 2:
            obs = -np.log10(np.sort(pvals.values))
            exp = -np.log10((np.arange(1, len(obs) + 1) - 0.5) / len(obs))
            qq_ax.scatter(exp, obs, s=10, alpha=0.6, color=color, label=pop)
            qq_lim = max(qq_lim, float(max(exp.max(), obs.max())))
            has_points = True

    if has_points:
        qq_ax.plot([0, qq_lim], [0, qq_lim], "k--", lw=1)
        qq_ax.set_xlabel("Expected -log10(p)")
        qq_ax.set_ylabel("Observed -log10(p)")
        qq_ax.set_title("EUR/AFR HetPVal QQ")
        qq_ax.legend()
    else:
        qq_ax.text(0.5, 0.5, "Insufficient HetPVal data", ha="center", va="center", transform=qq_ax.transAxes)

    # Panel J: Venn diagram of variant overlap
    try:
        from matplotlib_venn import venn3  # type: ignore[import-not-found]
        set_eur = set(df_harmonized["snp"])
        set_afr = set(raw_by_pop["AFR"]["snp"]) if "snp" in raw_by_pop["AFR"].columns else set()
        set_thai = set(raw_by_pop["THAI"]["snp"]) if "snp" in raw_by_pop["THAI"].columns else set()
        venn_ax = ax[3, 0]
        if len(set_eur) > 0 and len(set_afr) > 0 and len(set_thai) > 0:
            venn3([set_eur, set_afr, set_thai], set_labels=["EUR", "AFR", "THAI"], ax=venn_ax)
            venn_ax.set_title("Variant overlap among populations")
        else:
            venn_ax.text(0.5, 0.5, "Insufficient variant overlap", ha="center", va="center", transform=venn_ax.transAxes)
    except ImportError:
        ax[3, 0].text(0.5, 0.5, "matplotlib_venn not installed", ha="center", va="center", transform=ax[3, 0].transAxes)
    
    # Panel K: EUR cohort presence
    if "EUR" in raw_by_pop and "direction" in raw_by_pop["EUR"].columns:
        eur_dir = raw_by_pop["EUR"]["direction"]
        direction_counts = eur_dir.apply(lambda s: 7 - s.count("?") if isinstance(s, str) else 0)
        ax[3, 1].hist(direction_counts, bins=np.arange(-0.5, 8.5, 1), rwidth=0.8, edgecolor="black")
        ax[3, 1].set_xticks(range(1, 8))
        ax[3, 1].set_xlabel("Number of cohorts with data")
        ax[3, 1].set_ylabel("Count")
        ax[3, 1].set_title("EUR: Cohort presence based on direction")
    else:
        ax[3, 1].text(0.5, 0.5, "Missing direction column", ha="center", va="center", transform=ax[3, 1].transAxes)

    # Panel L: AFR cohort presence
    if "AFR" in raw_by_pop and "direction" in raw_by_pop["AFR"].columns:
        afr_dir = raw_by_pop["AFR"]["direction"]
        direction_counts = afr_dir.apply(lambda s: 3 - s.count("?") if isinstance(s, str) else 0)
        ax[3, 2].hist(direction_counts, bins=np.arange(-0.5, 4.5, 1), rwidth=0.8, edgecolor="black")
        ax[3, 2].set_xticks(range(1, 4))
        ax[3, 2].set_xlabel("Number of cohorts with data")
        ax[3, 2].set_ylabel("Count")
        ax[3, 2].set_title("AFR: Cohort presence based on direction")
    else:
        ax[3, 2].text(0.5, 0.5, "Missing direction column", ha="center", va="center", transform=ax[3, 2].transAxes)

    fig.tight_layout()
    fig.savefig(f"{output_prefix}_global_diagnostics.png", dpi=160)
    plt.close(fig)


def parse_thai_bim(bim_path):
    bim = pd.read_csv(
        bim_path,
        sep=r"\s+",
        header=None,
        names=["chr", "snp", "cm", "pos", "a1", "a2"],
        dtype={"chr": str, "snp": str, "cm": float, "pos": int, "a1": str, "a2": str},
    )
    bim["chr"] = bim["chr"].map(normalize_chr)
    bim["a1"] = bim["a1"].str.upper()
    bim["a2"] = bim["a2"].str.upper()
    return bim


def resolve_pan_ukbb_hail_ld_cache(pop, tag, local_ld_dir, chrom, start, end):
    """
    Resolve Pan-UKBB Hail LD cache files.

    Prefer exact tag match first. If not found, fall back to any cache whose
    variants span fully covers the requested [start, end] interval on `chrom`.
    """
    local_dir = Path(local_ld_dir)
    exact_matrix = local_dir / f"UKBB.{pop}.ldadj.{tag}.npy"
    exact_variants = local_dir / f"UKBB.{pop}.ldadj.{tag}.variants.tsv.gz"
    if exact_matrix.exists() and exact_variants.exists():
        return exact_matrix, exact_variants

    chrom_norm = normalize_chr(chrom)
    requested_start = int(start)
    requested_end = int(end)
    candidates = sorted(local_dir.glob(f"UKBB.{pop}.ldadj.*.npy"))
    covering = []

    for matrix_path in candidates:
        base_name = matrix_path.name[:-4]  # strip .npy
        variants_path = matrix_path.with_name(f"{base_name}.variants.tsv.gz")
        if not variants_path.exists():
            continue

        try:
            span_df = pd.read_csv(
                variants_path,
                sep="\t",
                compression="gzip",
                usecols=["chr", "pos"],
                low_memory=False,
            )
        except Exception:
            continue

        if "chr" not in span_df.columns or "pos" not in span_df.columns:
            continue
        span_df["chr"] = span_df["chr"].astype(str).map(normalize_chr)
        span_df["pos"] = pd.to_numeric(span_df["pos"], errors="coerce")

        pos_chr = span_df.loc[
            (span_df["chr"] == chrom_norm) & span_df["pos"].notna(),
            "pos",
        ]
        if pos_chr.empty:
            continue

        cand_start = int(pos_chr.min())
        cand_end = int(pos_chr.max())
        if cand_start <= requested_start and cand_end >= requested_end:
            span_width = cand_end - cand_start
            covering.append((span_width, cand_start, cand_end, matrix_path, variants_path))

    if covering:
        covering.sort(key=lambda x: (x[0], x[1], x[2], str(x[3])))
        _, cand_start, cand_end, matrix_path, variants_path = covering[0]
        print(
            f"[{pop}] Exact LD tag '{tag}' not found; using covering cache "
            f"{matrix_path.name} spanning chr{chrom_norm}:{cand_start}-{cand_end}"
        )
        return matrix_path, variants_path

    raise FileNotFoundError(
        f"Missing local LD cache for {pop}. Expected exact files:\n"
        f"  {exact_matrix}\n"
        f"  {exact_variants}\n"
        f"and found no fallback cache covering chr{chrom_norm}:{requested_start}-{requested_end}.\n"
        "Run step2a_extract_hail_LDs.py for this or a larger region first."
    )


def load_pan_ukbb_hail_ld_subset(pop, chrom, start, end, tag, local_ld_dir):
    """Load local Pan-UKBB Hail LD cache for a given ancestry label."""
    if local_ld_dir is None:
        raise ValueError(
            "--local-ld-dir is required for EUR/AFR LD loading. "
            "Generate local files first with step2a_extract_hail_LDs.py."
        )

    matrix_path, variants_path = resolve_pan_ukbb_hail_ld_cache(
        pop=pop,
        tag=tag,
        local_ld_dir=local_ld_dir,
        chrom=chrom,
        start=start,
        end=end,
    )

    print(f"[{pop}] Loading Pan-UKBB Hail LD from local cache: {matrix_path}")
    matrix = np.load(matrix_path)
    var_df = pd.read_csv(variants_path, sep="\t", compression="gzip")
    if matrix.ndim == 1:
        matrix = matrix.reshape(1, -1)
    if matrix.shape[0] != matrix.shape[1]:
        raise ValueError(f"Cached {pop} LD matrix is not square: {matrix.shape}")
    if len(var_df) != matrix.shape[0]:
        raise ValueError(
            f"Cached {pop} LD matrix size ({matrix.shape[0]}) does not match variants rows ({len(var_df)})."
        )

    var_df = var_df.copy()
    var_df["chr"] = var_df["chr"].astype(str).map(normalize_chr)
    var_df["a1"] = var_df["a1"].astype(str).str.upper()
    var_df["a2"] = var_df["a2"].astype(str).str.upper()
    var_df["snp"] = var_df["snp"].astype(str)
    print(f"[{pop}] Loaded LD for {len(var_df)} variants in region chr{chrom}:{start}-{end}")
    return matrix, var_df


def load_plink_ld_subset(pop, chrom, start, end, ld_path):
    """Load plink-style LD matrix + BIM and subset to a target region."""
    if pop != "THAI":
        raise ValueError(f"plink LD loader currently supports THAI only, got: {pop}")
    if not os.path.exists(ld_path):
        raise FileNotFoundError(f"{pop} LD matrix not found: {ld_path}")

    bim_path = os.path.splitext(ld_path)[0] + ".bim"
    if not os.path.exists(bim_path):
        raise FileNotFoundError(f"{pop} BIM not found (needed for SNP order): {bim_path}")

    matrix = np.loadtxt(ld_path)
    if matrix.ndim == 1:
        matrix = matrix.reshape(1, -1)
    if matrix.shape[0] != matrix.shape[1]:
        raise ValueError(f"{pop} LD matrix is not square.")

    bim_full = parse_thai_bim(bim_path)
    if len(bim_full) != matrix.shape[0]:
        raise ValueError(
            f"{pop} LD matrix size ({matrix.shape[0]}) does not match BIM rows ({len(bim_full)})."
        )

    mask = (bim_full["chr"] == chrom) & (bim_full["pos"] >= int(start)) & (bim_full["pos"] <= int(end))
    if mask.sum() == 0:
        return np.zeros((0, 0)), pd.DataFrame(columns=["idx", "chr", "pos", "snp", "a1", "a2"])

    idx = np.where(mask.to_numpy())[0]
    matrix_sub = matrix[np.ix_(idx, idx)]
    bim = bim_full.loc[mask, ["chr", "pos", "snp", "a1", "a2"]].copy().reset_index(drop=True)
    bim["idx"] = np.arange(len(bim), dtype=int)
    print(f"[{pop}] Loaded LD for {len(bim)} variants in region chr{chrom}:{start}-{end}")
    return matrix_sub, bim[["idx", "chr", "pos", "snp", "a1", "a2"]]


def resolve_ld_tag(local_ld_dir, pop, gene, flank_raw):
    """
    Resolve LD cache tag for a population.

    Accepts either:
      - {gene}_flank{flank}
      - {geneA}-{geneB}_flank{flank}, where gene is either geneA or geneB

    Resolution order:
      1. Exact single-gene tag with the requested flank.
      2. Exact cluster tag with the requested flank.
      3. Smallest available single-gene tag whose flank is >= requested flank.
      4. Smallest available cluster tag whose flank is >= requested flank.

    This allows a request like 300kb to reuse an existing 1Mb cache.
    """
    suffix = f"_flank{flank_raw}"
    gene_upper = str(gene).upper()
    expected_tag = f"{gene}{suffix}"
    requested_flank = parse_flank_size(flank_raw)

    matrix_path = os.path.join(local_ld_dir, f"UKBB.{pop}.ldadj.{expected_tag}.npy")
    variants_path = os.path.join(local_ld_dir, f"UKBB.{pop}.ldadj.{expected_tag}.variants.tsv.gz")
    if os.path.exists(matrix_path) and os.path.exists(variants_path):
        print(f"[{pop}] using exact LD tag: {expected_tag}")
        return expected_tag

    prefix = f"UKBB.{pop}.ldadj."
    pattern = os.path.join(local_ld_dir, f"{prefix}*.npy")
    exact_cluster_candidates = []
    larger_single_candidates = []
    larger_cluster_candidates = []
    all_matching_tags = []

    for npy_path in sorted(glob.glob(pattern)):
        base = os.path.basename(npy_path)
        if not (base.startswith(prefix) and base.endswith(".npy")):
            continue

        tag = base[len(prefix):-4]
        cand_variants = os.path.join(local_ld_dir, f"UKBB.{pop}.ldadj.{tag}.variants.tsv.gz")
        if not os.path.exists(cand_variants):
            continue
        if "_flank" not in tag:
            continue

        gene_part, tag_flank_raw = tag.rsplit("_flank", 1)
        try:
            tag_flank = parse_flank_size(tag_flank_raw)
        except argparse.ArgumentTypeError:
            continue

        if gene_part.upper() == gene_upper:
            match_kind = "single"
        else:
            pieces = gene_part.split("-")
            if len(pieces) != 2:
                continue
            if gene_upper not in {p.upper() for p in pieces}:
                continue
            match_kind = "cluster"

        all_matching_tags.append(tag)

        if tag_flank == requested_flank and match_kind == "cluster":
            exact_cluster_candidates.append(tag)
            continue

        if tag_flank < requested_flank:
            continue

        candidate_record = (tag_flank, tag)
        if match_kind == "single":
            larger_single_candidates.append(candidate_record)
        else:
            larger_cluster_candidates.append(candidate_record)

    if len(exact_cluster_candidates) == 1:
        print(f"[{pop}] using exact cluster LD tag: {exact_cluster_candidates[0]}")
        return exact_cluster_candidates[0]

    if len(exact_cluster_candidates) > 1:
        raise RuntimeError(
            f"Multiple exact LD cluster tags matched for {pop}, gene={gene}, flank={flank_raw}: {exact_cluster_candidates}. "
            "Please keep only one matching cache set or rename files to avoid ambiguity."
        )

    if larger_single_candidates:
        larger_single_candidates.sort(key=lambda x: (x[0], x[1]))
        chosen_flank, chosen_tag = larger_single_candidates[0]
        print(
            f"[{pop}] exact LD tag '{expected_tag}' not found; using larger single-gene cache tag "
            f"'{chosen_tag}' (flank {chosen_flank} bp >= requested {requested_flank} bp)"
        )
        return chosen_tag

    if larger_cluster_candidates:
        larger_cluster_candidates.sort(key=lambda x: (x[0], x[1]))
        chosen_flank, chosen_tag = larger_cluster_candidates[0]
        print(
            f"[{pop}] exact LD tag '{expected_tag}' not found; using larger cluster cache tag "
            f"'{chosen_tag}' (flank {chosen_flank} bp >= requested {requested_flank} bp)"
        )
        return chosen_tag

    available_text = ", ".join(all_matching_tags) if all_matching_tags else "none"
    raise FileNotFoundError(
        f"Could not resolve LD tag for {pop}. Tried exact tag '{expected_tag}' and searched for single/cluster tags "
        f"covering gene '{gene}' with flank >= {flank_raw}. Available matching tags: {available_text}."
    )


def align_ld_to_reference(matrix, var_df, ref_df, label=""):
    if matrix.shape[0] == 0 or var_df.empty or ref_df.empty:
        if label:
            print(f"[align:{label}] empty input encountered (LD rows={matrix.shape[0]}, LD variants={len(var_df)}, ref variants={len(ref_df)})")
        empty = ref_df.iloc[0:0].copy()
        return np.zeros((0, 0)), empty

    merged = ref_df[["snp", "a1", "a2"]].merge(
        var_df[["idx", "snp", "a1", "a2"]],
        on="snp",
        how="inner",
        suffixes=("_ref", "_ld"),
    )
    if label:
        print(f"[align:{label}] SNP ID overlap ref vs LD: {len(ref_df)} vs {len(var_df)} -> {len(merged)} shared IDs")
    same = (merged["a1_ref"] == merged["a1_ld"]) & (merged["a2_ref"] == merged["a2_ld"])
    flip = (merged["a1_ref"] == merged["a2_ld"]) & (merged["a2_ref"] == merged["a1_ld"])
    if label:
        same_n = int(same.sum())
        flip_n = int(flip.sum())
        mismatch_n = int(len(merged) - same_n - flip_n)
        print(f"[align:{label}] allele status among shared IDs: same={same_n}, flip={flip_n}, mismatch_drop={mismatch_n}")
    merged = merged[same | flip].copy()
    if merged.empty:
        if label:
            print(f"[align:{label}] no SNPs remain after allele harmonization")
        empty = ref_df.iloc[0:0].copy()
        return np.zeros((0, 0)), empty

    merged["flip_sign"] = np.where(flip.loc[merged.index], -1.0, 1.0)
    merged = merged.drop_duplicates(subset=["snp"], keep="first")

    idx_map = dict(zip(merged["snp"], merged["idx"].astype(int)))
    flip_map = dict(zip(merged["snp"], merged["flip_sign"].astype(float)))

    keep = ref_df[ref_df["snp"].isin(idx_map)].copy()
    ld_idx = keep["snp"].map(idx_map).to_numpy(dtype=int)
    signs = keep["snp"].map(flip_map).to_numpy(dtype=float)
    sub = matrix[np.ix_(ld_idx, ld_idx)]
    sub = (signs[:, None] * sub) * signs[None, :]
    if label:
        print(f"[align:{label}] retained {len(keep)} SNPs for LD-aligned output")
    return sub, keep.reset_index(drop=True)


def regional_diagnostics(df_region, ld_by_pop, output_prefix):
    if df_region.empty:
        return

    fig, ax = plt.subplots(2, 2, figsize=(12, 10))

    # p-vals vs position with dual y-axes: EUR on left, AFR/THAI on right.
    a_main = ax[0, 0]
    a_sec = a_main.twinx()
    pos = df_region["pos"]
    y_eur = -np.log10(np.maximum(df_region["p_EUR"], 1e-300))
    y_afr = -np.log10(np.maximum(df_region["p_AFR"], 1e-300))
    y_thai = -np.log10(np.maximum(df_region["p_THAI"], 1e-300))

    h_eur = a_main.scatter(pos, y_eur, s=12, alpha=0.55, color="tab:blue", label="EUR")
    h_afr = a_sec.scatter(pos, y_afr, s=12, alpha=0.55, color="tab:orange", label="AFR")
    h_thai = a_sec.scatter(pos, y_thai, s=12, alpha=0.55, color="tab:green", label="THAI")

    eur_max = float(np.nanmax(y_eur)) if np.isfinite(np.nanmax(y_eur)) else 1.0
    afr_thai_max = float(np.nanmax(np.r_[y_afr.to_numpy(dtype=float), y_thai.to_numpy(dtype=float)]))
    if not np.isfinite(eur_max) or eur_max <= 0:
        eur_max = 1.0
    if not np.isfinite(afr_thai_max) or afr_thai_max <= 0:
        afr_thai_max = 1.0

    eur_ylim_max = max(1.0, float(np.ceil(eur_max)))
    afr_thai_ylim_max = max(1.0, float(np.ceil(afr_thai_max)))
    a_main.set_ylim(0, eur_ylim_max)
    a_sec.set_ylim(0, afr_thai_ylim_max)

    a_main.set_xlabel("Position")
    a_main.set_ylabel("EUR -log10(p)", color="tab:blue")
    a_sec.set_ylabel("AFR/THAI -log10(p)", color="tab:orange")
    a_main.set_title("Regional association (dual y-axis)")
    a_main.tick_params(axis="y", labelcolor="tab:blue")
    a_sec.tick_params(axis="y", labelcolor="tab:orange")
    a_main.legend([h_eur, h_afr, h_thai], ["EUR", "AFR", "THAI"], loc="upper right")

    lead = df_region.loc[df_region["p_meta"].idxmin(), "snp"] if "p_meta" in df_region.columns else df_region.loc[df_region["p_EUR"].idxmin(), "snp"]

    for (p1, p2), a in zip([("EUR", "AFR"), ("EUR", "THAI"), ("AFR", "THAI")],
                           [ax[0, 1], ax[1, 0], ax[1, 1]]):
        r1 = ld_by_pop[p1]
        r2 = ld_by_pop[p2]
        if r1.shape[0] == 0 or r2.shape[0] == 0 or r1.shape != r2.shape:
            a.set_title(f"{p1} vs {p2} LD mismatch (insufficient shared SNPs)")
            continue
        if lead not in df_region["snp"].values:
            a.set_title(f"{p1} vs {p2} LD mismatch (lead SNP missing)")
            continue
        lead_idx = int(df_region.index[df_region["snp"] == lead][0])
        plot_df = pd.DataFrame({
            "pos": df_region["pos"].to_numpy(),
            "r2_p1": np.square(r1[:, lead_idx]),
            "r2_p2": np.square(r2[:, lead_idx]),
            "c": (
                -np.log10(np.maximum(df_region["p_meta"], 1e-300)).to_numpy()
                if "p_meta" in df_region.columns
                else -np.log10(np.maximum(df_region["p_EUR"], 1e-300)).to_numpy()
            ),
        })

        # Remove variants that are in LE with lead SNP in both populations.
        plot_df = plot_df[~((plot_df["r2_p1"] == 0) & (plot_df["r2_p2"] == 0))].copy()
        plot_df = plot_df[np.isfinite(plot_df["r2_p1"]) & np.isfinite(plot_df["r2_p2"]) & np.isfinite(plot_df["c"])].copy()
        if plot_df.empty:
            a.set_title(f"{p1} vs {p2} LD mismatch (all points r2=0 in both)")
            continue

        a.scatter(plot_df["r2_p1"], plot_df["r2_p2"], c=plot_df["c"], s=20, alpha=0.7, cmap="viridis")
        a.set_xlabel(f"r2 with {lead} (lead) in {p1}")
        a.set_ylabel(f"r2 with {lead} (lead) in {p2}")
        a.set_title(f"LD mismatch: {p1} vs {p2}")

    fig.tight_layout()
    fig.savefig(f"{output_prefix}_regional_diagnostics.png", dpi=160)
    plt.close(fig)


def write_finemap_files(df, ld_matrix, prefix, pop, 
                        finemap_master=False,
                        finemap_outdir=None):
    out_z = f"{prefix}.{pop}.z"
    out_ld = f"{prefix}.{pop}.ld"
    z = pd.DataFrame({
        "rsid": df["snp"],
        "chromosome": df["chr"],
        "position": df["pos"].astype(int),
        "allele1": df["a1"],
        "allele2": df["a2"],
        "maf": df[f"maf_{pop}"],
        "beta": df[f"beta_{pop}"],
        "se": df[f"se_{pop}"],
        "p": df[f"p_{pop}"],
        "n": df[f"n_{pop}"],
    })
    z.to_csv(out_z, sep="\t", index=False)
    np.savetxt(out_ld, ld_matrix, fmt="%.8g", delimiter="\t")
    print(f"Wrote {out_z} and {out_ld}")
    # also create the .master file for FINEMAP if commanded to:
    if finemap_master:
        assert finemap_outdir is not None, \
            "finemap_outdir must be specified if finemap_master is True"
        prefix_path = Path(prefix)
        master_path = Path(f"{prefix}.{pop}.master")
        prefix_basename = prefix_path.name
        finemap_dir = Path(finemap_outdir).expanduser().resolve() 

        filename_list = [
            Path(f"{prefix}.{pop}.z").expanduser().resolve().as_posix(),
            Path(f"{prefix}.{pop}.ld").expanduser().resolve().as_posix(),
            (finemap_dir / f"{prefix_basename}.{pop}.snp").as_posix(),
            (finemap_dir / f"{prefix_basename}.config").as_posix(),
            (finemap_dir / f"{prefix_basename}.cred").as_posix(),
            (finemap_dir / f"{prefix_basename}.log").as_posix(),
            str(df[f"n_{pop}"].median()),
        ]
        write_header = (not master_path.exists()) or (master_path.stat().st_size == 0)
        with open(master_path, "a") as f:
            if write_header:
                f.write("z;ld;snp;config;cred;log;n_samples\n")
            f.write(";".join(filename_list) + "\n")


def write_single_ancestry_finemap_files(region_by_pop, ld_raw_by_pop, output_prefix, finemap_outdir=None):
    """
    Write single-ancestry finemapping inputs for each population.

    Each output uses that ancestry's own GWAS region table and LD SNP set,
    so SNP retention is maximized before cross-ancestry harmonization.
    """
    for pop in ["EUR", "AFR", "THAI"]:
        ld_mat_raw, ld_var = ld_raw_by_pop[pop]
        region_df = region_by_pop[pop]
        ld_mat, keep = align_ld_to_reference(
            ld_mat_raw,
            ld_var,
            region_df[["snp", "a1", "a2", "chr", "pos"]],
            label=f"single-{pop}",
        )
        if keep.empty or ld_mat.shape[0] == 0:
            print(f"Skipping single-ancestry {pop}: no LD-aligned SNPs")
            continue

        # Keep z rows in the same order as LD matrix rows.
        single_df = keep[["snp"]].merge(region_df, on="snp", how="inner")
        if len(single_df) != ld_mat.shape[0]:
            raise RuntimeError(
                f"Single-ancestry {pop} row mismatch: z rows={len(single_df)}, LD rows={ld_mat.shape[0]}"
            )

        # Reuse the common writer by projecting ancestry-specific fields.
        single_df[f"maf_{pop}"] = single_df["maf"]
        single_df[f"beta_{pop}"] = single_df["beta"]
        single_df[f"se_{pop}"] = single_df["se"]
        single_df[f"p_{pop}"] = single_df["p"]
        single_df[f"n_{pop}"] = single_df["n"]

        pop_output_prefix = output_prefix
        if pop in {"EUR", "AFR"}:
            out_path = Path(output_prefix)
            # Remove ancestry subset marker (e.g., _subCSA) from EUR/AFR single-ancestry filenames.
            cleaned_name = re.sub(r"_sub[^._/]+", "", out_path.name)
            pop_output_prefix = str(out_path.with_name(cleaned_name))

        write_finemap_files(single_df, ld_mat, 
                            f"{pop_output_prefix}.single", pop,
                            finemap_master=True,
                            finemap_outdir=finemap_outdir)


def load_harmonized_table(path):
    """Load and coerce harmonized table columns to expected dtypes."""
    harm = read_table(path)
    for col in harm.columns:
        if col in ["chr"]:
            harm[col] = harm[col].astype(str).map(normalize_chr)
        elif col in ["pos", "n_EUR", "n_AFR", "n_THAI"]:
            harm[col] = pd.to_numeric(harm[col], errors="coerce")
        elif col in [
            "maf_EUR", "maf_AFR", "maf_THAI",
            "beta_EUR", "beta_AFR", "beta_THAI",
            "se_EUR", "se_AFR", "se_THAI",
            "p_EUR", "p_AFR", "p_THAI",
            "beta_meta", "z_meta", "p_meta",
            "Q", "I2", "tau2", "HetP",
        ]:
            harm[col] = pd.to_numeric(harm[col], errors="coerce")
    return harm


def main():
    parser = argparse.ArgumentParser(description="Prepare LD matrices and diagnostics for fine-mapping")
    parser.add_argument("--gene", required=True, help="Gene name to fine-map (e.g. BACH2)")
    parser.add_argument('--flank', default='500000',
                        help='Flank size in bp (default: 500000). Supports: 500000, 1e6, 1Mb, 500kb, etc.')
    parser.add_argument("--input-files", "-i", nargs="+", required=False, help="Input GWAS summary statistics files (one per population, can be gzipped) and the corresponding ancestry group, separated by comma \"filename,ANC\"")
    parser.add_argument("--ref-genome", choices=["GRCh37", "GRCh38"], default="GRCh37", help="Reference genome version (default: GRCh37)")
    parser.add_argument("--output-prefix", required=True, help="Prefix for output files (e.g. BACH2_finemap)")
    parser.add_argument("--thai-ld-path", default=THAI_LD_PATH, help="Path to Thai LD matrix file (.ld)")
    parser.add_argument(
        "--finemap-outdir", default=None, required=True,
        help=(
            "Directory for FINEMAP artifacts referenced by .master files "
            "(.snp/.config/.cred/.log)."
        ),
    )
    parser.add_argument(
        "--local-ld-dir", required=True,
        help=("Directory containing pre-fetched Pan-UKBB LD files produced by step2a_extract_hail_LDs.py "
            "(e.g. UKBB.EUR.<chrom>_<start>_<end>.ld.npy). "
            "Required: this script only loads EUR/AFR LD from local cache."
        ),
    )
    parser.add_argument(
        "--skip-global-diagnostics", action="store_true",
        help="Skip generation of global diagnostic plots (default: generate plots)"
    )
    parser.add_argument(
        "--use-harmonized-only", action="store_true",
        help=("Skip per-population GWAS loading/standardization and load directly from a precomputed "
            "3-way harmonized file."
        ),
    )
    parser.add_argument(
        "--harmonized-path", default=None,
        help=("Path to precomputed 3-way harmonized table (.tsv/.tsv.gz). "
            "Defaults to <output-prefix>.harmonized.tsv.gz."
        ),
    )
    
    args = parser.parse_args()
    args.flank_raw = args.flank          # original string (e.g. "1Mb", "500000")
    args.flank = parse_flank_size(args.flank_raw)

    if args.use_harmonized_only and args.input_files:
        print("[mode] --use-harmonized-only is set; ignoring --input-files.")
    if (not args.use_harmonized_only) and (not args.input_files):
        parser.error("--input-files is required unless --use-harmonized-only is specified.")

    t0 = time.time()
    outdir = os.path.dirname(args.output_prefix)
    if outdir:
        os.makedirs(outdir, exist_ok=True)
        
    os.makedirs(args.finemap_outdir, exist_ok=True)

    gene_chr, gene_start, gene_end = get_gene_coordinates(
        args.gene,
        gene_ref_path=b37_GENE_REF if args.ref_genome == "GRCh37" else b38_GENE_REF,
    )
    region_start = max(1, int(gene_start - args.flank))
    region_end = int(gene_end + args.flank)
    print(f"Using region chr{gene_chr}:{region_start}-{region_end} (gene +/- flank)")

    target_chr = normalize_chr(gene_chr)
    region_by_pop = {}
    raw = {}
    if not args.use_harmonized_only:
        files = parse_input_files(args.input_files)
        for pop in ["EUR", "AFR", "THAI"]:
            print(f"Reading {pop} GWAS: {files[pop]}")
            raw_df = read_table(files[pop])
            print(f"{pop} variants before standardization: {len(raw_df)}")
            print(raw_df.head(2))
            raw[pop] = standardize_gwas(raw_df, pop)
            print(f"{pop}: {len(raw[pop])} variants after standardization")

        # Per-ancestry region tables (before cross-ancestry harmonization).
        for pop in ["EUR", "AFR", "THAI"]:
            d = raw[pop]
            region_by_pop[pop] = d[
                (d["chr"] == target_chr) &
                (d["pos"] >= region_start) &
                (d["pos"] <= region_end)
            ].copy()
            print(f"{pop} regional variants before harmonization: {len(region_by_pop[pop])}")

    # Load regional LD once and reuse for both single- and cross-ancestry outputs.
    local_ld_dir = args.local_ld_dir
    resolved_tags = {}
    for pop in ["EUR", "AFR"]:
        resolved_tags[pop] = resolve_ld_tag(
            local_ld_dir=local_ld_dir,
            pop=pop,
            gene=args.gene,
            flank_raw=args.flank_raw,
        )
    print("Resolved LD cache tags:")
    for pop in ["EUR", "AFR"]:
        print(f"  {pop}: {resolved_tags[pop]}")

    local_ld_dir_norm = os.path.normpath(os.path.abspath(os.path.expanduser(local_ld_dir)))
    thai_ld_dir_norm = os.path.normpath(
        os.path.abspath(os.path.expanduser(os.path.dirname(args.thai_ld_path)))
    )
    use_thai_hail_cache = (thai_ld_dir_norm == local_ld_dir_norm)
    thai_substitute_pop = None
    if use_thai_hail_cache:
        # get substitute pop label from thai_ld_path filename for messaging
        thai_substitute_pop = os.path.basename(args.thai_ld_path).split(".")[1]
        print(
            f"THAI LD mode: Pan-UKBB Hail cache fallback ({thai_substitute_pop}) because "
            f"dirname(--thai-ld-path) matches --local-ld-dir: {local_ld_dir_norm}"
        )
    else:
        print(
            f"THAI LD mode: plink matrix from --thai-ld-path ({args.thai_ld_path}); "
            f"dirname mismatch with --local-ld-dir"
        )
    # tag = f"{target_chr}_{region_start}_{region_end}"
    # tag = f"{args.gene}_flank{args.flank_raw}"
    ld_raw = {}
    for pop in ["EUR", "AFR", "THAI"]:
        print(f"Loading {pop} LD")
        if pop in {"EUR", "AFR"}:
            tag = resolved_tags[pop]
            mat, var = load_pan_ukbb_hail_ld_subset(
                pop,
                chrom=target_chr,
                start=region_start,
                end=region_end,
                tag=tag,
                local_ld_dir=local_ld_dir,
            )
        # If --thai-ld-path is inside --local-ld-dir, 
        # attempt THAI fallback from the same Hail cache format.
        elif use_thai_hail_cache and (pop == "THAI"):
            substitute_pop = thai_substitute_pop
            if not substitute_pop:
                raise ValueError("THAI Hail cache fallback requested but substitute population could not be inferred from --thai-ld-path")
            thai_base = os.path.basename(args.thai_ld_path)
            thai_prefix = f"UKBB.{substitute_pop}.ldadj."
            if not (thai_base.startswith(thai_prefix) and thai_base.endswith(".npy")):
                raise ValueError(
                    f"Invalid --thai-ld-path filename for Hail cache fallback: {thai_base}. "
                    f"Expected pattern '{thai_prefix}<tag>.npy'."
                )
            thai_tag = thai_base[len(thai_prefix):-4]
            print(f"Using {substitute_pop} LD from Hail cache for THAI")
            mat, var = load_pan_ukbb_hail_ld_subset(
                substitute_pop, 
                chrom=target_chr,
                start=region_start,
                end=region_end,
                tag=thai_tag,
                local_ld_dir=local_ld_dir,
            )
        else:
            mat, var = load_plink_ld_subset(
                pop,
                chrom=target_chr,
                start=region_start,
                end=region_end,
                ld_path=args.thai_ld_path,
            )
        ld_raw[pop] = (mat, var)
        print(f"{pop} LD variants loaded: {len(var)}")

    # Write per-ancestry finemap-ready files prior to harmonization.
    if not args.use_harmonized_only:
        write_single_ancestry_finemap_files(region_by_pop, ld_raw, args.output_prefix, finemap_outdir=args.finemap_outdir)
    else:
        print("Skipping single-ancestry file writing in harmonized-only mode")

    # Check for pre-existing harmonized file when skipping diagnostics
    harmonized_path = args.harmonized_path or f"{args.output_prefix}.harmonized.tsv.gz"
    if args.use_harmonized_only:
        if not os.path.exists(harmonized_path):
            raise FileNotFoundError(
                f"Harmonized-only mode requested but file not found: {harmonized_path}"
            )
        print(f"[mode] Harmonized-only: loading {harmonized_path}")
        harm = load_harmonized_table(harmonized_path)
        print(f"Loaded {len(harm)} harmonized variants from file")
    elif args.skip_global_diagnostics and os.path.exists(harmonized_path):
        print(f"Loading pre-existing harmonized file: {harmonized_path}")
        harm = load_harmonized_table(harmonized_path)
        print(f"Loaded {len(harm)} harmonized variants from file")
    else:
        harm = harmonize_three_way(raw)
        if args.skip_global_diagnostics:
            print("Skipping heterogeneity/meta-stat computation (--skip-global-diagnostics); regional diagnostics will use EUR p-values")
        else:
            harm = compute_heterogeneity(harm)
            with np.errstate(divide="ignore"):
                harm["z_meta"] = harm["beta_meta"] / np.sqrt(
                    1.0 / (1.0 / np.square(harm["se_EUR"]) + 1.0 / np.square(harm["se_AFR"]) + 1.0 / np.square(harm["se_THAI"]))
                )
            harm["p_meta"] = 2.0 * stats.norm.sf(np.abs(harm["z_meta"]))
        harm = harm.sort_values(["chr", "pos"], key=lambda s: s.map(chr_sort_key) if s.name == "chr" else s).reset_index(drop=True)
        harm.to_csv(harmonized_path, sep="\t", index=False, compression="gzip")
        print(f"Harmonized variants across EUR/AFR/THAI: {len(harm)}")
    

    if args.use_harmonized_only and not args.skip_global_diagnostics:
        print("Global diagnostics require raw per-population GWAS inputs; skipping in harmonized-only mode")
    elif not args.skip_global_diagnostics:
        global_diagnostics(harm, raw, args.output_prefix)
    else:
        print("Skipping global diagnostic plots (--skip-global-diagnostics)")

    region = harm[
        (harm["chr"] == normalize_chr(gene_chr)) &
        (harm["pos"] >= region_start) &
        (harm["pos"] <= region_end)
    ].copy()
    region.to_csv(f"{args.output_prefix}.region.tsv.gz", sep="\t", index=False, compression="gzip")
    # report how much variants is in the region after harmonization
    print(f"Regional variants in gene +/- flank: {len(region)}")
    print(region.head(2))
    if region_by_pop:
        print("Regional retention after 3-way harmonization:")
        for pop in ["EUR", "AFR", "THAI"]:
            pre_n = len(region_by_pop[pop])
            shared_n = int(region_by_pop[pop]["snp"].isin(region["snp"]).sum())
            print(f"  {pop}: regional pre-harmonization={pre_n}, retained in 3-way region={shared_n}, dropped={pre_n - shared_n}")

    if len(region) == 0:
        raise RuntimeError("No harmonized variants found in target region.")

    aligned = {}
    keep_sets = []
    for pop in ["EUR", "AFR", "THAI"]:
        mat, var = ld_raw[pop]
        sub, keep = align_ld_to_reference(
            mat,
            var,
            region[["snp", "a1", "a2", "chr", "pos"]],
            label=f"joint-{pop}",
        )
        aligned[pop] = (sub, keep)
        keep_sets.append(set(keep["snp"]))
        print(f"{pop} LD aligned SNPs: {len(keep)}")

    common_snps = set(region["snp"])
    print(f"[common] starting from harmonized regional SNPs: {len(common_snps)}")
    for pop, keep_set in zip(["EUR", "AFR", "THAI"], keep_sets):
        overlap_n = len(common_snps & keep_set)
        print(f"[common] overlap with {pop} LD-aligned SNPs: {overlap_n} / {len(common_snps)}")
    for s in keep_sets:
        common_snps &= s
    print(f"[common] final SNPs shared by harmonized region and all 3 LD panels: {len(common_snps)}")
    if len(common_snps) == 0:
        raise RuntimeError("No common SNPs across region GWAS and all three LD panels after allele harmonization.")

    region_common = region[region["snp"].isin(common_snps)].copy()
    region_common = region_common.sort_values(["chr", "pos"], key=lambda s: s.map(chr_sort_key) if s.name == "chr" else s).reset_index(drop=True)

    ld_final = {}
    for pop in ["EUR", "AFR", "THAI"]:
        mat, keep = aligned[pop]
        idx_map = {snp: i for i, snp in enumerate(keep["snp"].tolist())}
        ix = np.array([idx_map[s] for s in region_common["snp"]], dtype=int)
        ld_final[pop] = mat[np.ix_(ix, ix)]

    regional_diagnostics(region_common, ld_final, args.output_prefix)

    region_common.to_csv(f"{args.output_prefix}.region_common.tsv.gz", sep="\t", index=False, compression="gzip")
    for pop in ["EUR", "AFR", "THAI"]:
        write_finemap_files(region_common, ld_final[pop], 
                            args.output_prefix, pop,
                            finemap_master=False)

    dt = time.time() - t0
    print(f"Done in {dt:.1f} seconds")


if __name__ == "__main__":
    main()
