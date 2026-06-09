#!/usr/bin/env python3
import argparse
import pandas as pd
import numpy as np
import os
import gzip
import glob
import time
import matplotlib.pyplot as plt
from scipy import stats

# set global parameters
ROOT = "/lab-share/Hem-Sankaran-e2/Public"
b37_GENE_REF = os.path.join(ROOT, "ref_genomes/human/hg19_ncbiRefSeq_RefSeqAll_2026UCSC.tsv.gz")
b38_GENE_REF = os.path.join(ROOT, "ref_genomes/human/hg38_ncbiRefSeq_AllGenes_2026UCSC.tsv.gz")


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
    required = {"EUR", "AFR"}
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


def harmonize_two_way(gwas_by_pop):
    ref = gwas_by_pop["EUR"][["snp", "chr", "pos", "a1", "a2", "beta", "se", "p", "maf", "n"]].copy()
    ref = ref.rename(columns={
        "a1": "a1_ref", "a2": "a2_ref", "beta": "beta_EUR", "se": "se_EUR", "p": "p_EUR", "maf": "maf_EUR", "n": "n_EUR"
    })

    out = ref
    print(f"[harmonize] EUR baseline SNPs: {len(out)}")
    pop = "AFR"
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
    print(f"[harmonize] final 2-way harmonized SNPs: {len(out)}")
    return out


def compute_heterogeneity(df):
    b = df[["beta_EUR", "beta_AFR"]].to_numpy(dtype=float)
    se = df[["se_EUR", "se_AFR"]].to_numpy(dtype=float)
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
    
    required_pops = {"EUR", "AFR"}
    missing_pops = required_pops - set(raw_by_pop.keys())
    if missing_pops:
        print(f"Warning: missing populations {missing_pops}; skipping global diagnostics")
        return
    
    fig, ax = plt.subplots(3, 3, figsize=(14, 14))

    # Panel A, B: Sign-adjusted beta vs MAF for each population
    for i, (pop, color) in enumerate([("EUR", "tab:blue"), 
                                      ("AFR", "tab:green")]):
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

    # Panel C: EUR vs AFR beta correlation
    pairs = [("AFR", "EUR")]
    for (p1, p2), a in zip(pairs, [ax[0, 2]]):
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
    # ax[1, 1].axis("off")
    
    # Panel D: I^2 distribution by population
    has_i2_data = all(pop in raw_by_pop and "I2" in raw_by_pop[pop].columns for pop in ["EUR", "AFR"])
    if has_i2_data:
        eur = raw_by_pop["EUR"]
        afr = raw_by_pop["AFR"]
        eur_i2 = eur["I2"].dropna()
        afr_i2 = afr["I2"].dropna()
        
        if len(eur_i2) > 0 and len(afr_i2) > 0:
            ax[1, 0].hist([eur_i2, afr_i2], bins=np.arange(0, 110, 10), density=True, alpha=0.5, 
                          label=["EUR", "AFR"], color=["tab:blue", "tab:green"])
            ax[1, 0].set_xlabel(r"$I^2$")
            ax[1, 0].set_ylabel("Density")
            # Annotate statistics
            for pop, color in [("EUR", "tab:blue"), ("AFR", "tab:green")]:
                d = raw_by_pop[pop]
                i2_vals = d["I2"].dropna()
                if len(i2_vals) > 0:
                    median_i2 = np.median(i2_vals)
                    prop_i2_gt_50 = (i2_vals > 50).mean()
                    y_pos = 0.95 - 0.1 * list(sorted(raw_by_pop.keys())).index(pop)
                    ax[1, 0].text(0.95, y_pos, f"{pop}: median={median_i2:.1f}, >50%={prop_i2_gt_50:.1%}",
                                transform=ax[1, 0].transAxes, ha="right", va="top", color=color, fontsize=9)
            ax[1, 0].set_title(r"$I^2$ distribution")
            ax[1, 0].legend()
        else:
            ax[1, 0].text(0.5, 0.5, r"Insufficient $I^2$ data", ha="center", va="center", transform=ax[1, 0].transAxes)
    else:
        ax[1, 0].text(0.5, 0.5, r"Missing $I^2$ columns", ha="center", va="center", transform=ax[1, 0].transAxes)

    # Panel E: I^2 vs association strength
    if has_i2_data:
        for pop, color in [("EUR", "tab:blue"), ("AFR", "tab:green")]:
            d = raw_by_pop[pop]
            if "I2" in d.columns and "p" in d.columns:
                valid = ~(d["I2"].isna() | d["p"].isna())
                if valid.sum() > 0:
                    ax[2, 1].scatter(d.loc[valid, "I2"], -np.log10(np.maximum(d.loc[valid, "p"], 1e-300)), 
                                   s=10, alpha=0.4, label=pop, color=color)
        ax[1, 1].set_xlabel(r"$I^2$")
        ax[1, 1].set_ylabel(r"-\log_{10}(p)")
        ax[1, 1].set_title(r"$I^2$ vs association strength")
        if len(ax[1, 1].lines) > 0 or len(ax[1, 1].collections) > 0:
            ax[1, 1].legend()
    else:
        ax[1, 1].text(0.5, 0.5, r"Missing $I^2$ data", ha="center", va="center", transform=ax[1, 1].transAxes)

    # Panel F: Heterogeneity p-value QQ plot
    qq_ax = ax[1, 2]
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

    # Panel G: Venn diagram of variant overlap
    try:
        from matplotlib_venn import venn2 # type: ignore[import-not-found]
        set_eur = set(df_harmonized["snp"])
        set_afr = set(raw_by_pop["AFR"]["snp"]) if "snp" in raw_by_pop["AFR"].columns else set()
        venn_ax = ax[2, 0]
        if len(set_eur) > 0 and len(set_afr) > 0:
            venn2([set_eur, set_afr], set_labels=["EUR", "AFR"], ax=venn_ax)
            venn_ax.set_title("Variant overlap between populations")
        else:
            venn_ax.text(0.5, 0.5, "Insufficient variant overlap", ha="center", va="center", transform=venn_ax.transAxes)
    except ImportError:
        ax[2, 0].text(0.5, 0.5, "matplotlib_venn not installed", ha="center", va="center", transform=ax[2, 0].transAxes)
    
    # Panel I: EUR cohort presence
    if "EUR" in raw_by_pop and "direction" in raw_by_pop["EUR"].columns:
        eur_dir = raw_by_pop["EUR"]["direction"]
        direction_counts = eur_dir.apply(lambda s: 7 - s.count("?") if isinstance(s, str) else 0)
        ax[2, 1].hist(direction_counts, bins=np.arange(-0.5, 8.5, 1), rwidth=0.8, edgecolor="black")
        ax[2, 1].set_xticks(range(1, 8))
        ax[2, 1].set_xlabel("Number of cohorts with data")
        ax[2, 1].set_ylabel("Count")
        ax[2, 1].set_title("EUR: Cohort presence based on direction")
    else:
        ax[2, 1].text(0.5, 0.5, "Missing direction column", ha="center", va="center", transform=ax[2, 1].transAxes)

    # Panel J: AFR cohort presence
    if "AFR" in raw_by_pop and "direction" in raw_by_pop["AFR"].columns:
        afr_dir = raw_by_pop["AFR"]["direction"]
        direction_counts = afr_dir.apply(lambda s: 3 - s.count("?") if isinstance(s, str) else 0)
        ax[2, 2].hist(direction_counts, bins=np.arange(-0.5, 4.5, 1), rwidth=0.8, edgecolor="black")
        ax[2, 2].set_xticks(range(1, 4))
        ax[2, 2].set_xlabel("Number of cohorts with data")
        ax[2, 2].set_ylabel("Count")
        ax[2, 2].set_title("AFR: Cohort presence based on direction")
    else:
        ax[2, 2].text(0.5, 0.5, "Missing direction column", ha="center", va="center", transform=ax[2, 2].transAxes)

    fig.tight_layout()
    fig.savefig(f"{output_prefix}_global_diagnostics.png", dpi=160)
    plt.close(fig)


def load_ld_subset(pop, chrom, start, end, tag, ref_genome="GRCh37", local_ld_dir=None): # 
    print(f"Loading {pop} LD")
    if pop in {"EUR", "AFR"}:
        if local_ld_dir is None:
            raise ValueError(
                "--local-ld-dir is required for EUR/AFR LD loading. "
                "Generate local files first with step2a_extract_hail_LDs.py."
            )
        # tag = f"{normalize_chr(chrom)}_{start}_{end}"
        matrix_path = os.path.join(local_ld_dir, f"UKBB.{pop}.ldadj.{tag}.npy")
        variants_path = os.path.join(local_ld_dir, f"UKBB.{pop}.ldadj.{tag}.variants.tsv.gz")
        if not os.path.exists(matrix_path) or not os.path.exists(variants_path):
            raise FileNotFoundError(
                f"Missing local LD cache for {pop}. Expected files:\n"
                f"  {matrix_path}\n"
                f"  {variants_path}\n"
                "Run step2a_extract_hail_LDs.py for the same region first."
            )

        print(f"[{pop}] Loading LD from local cache: {matrix_path}")
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

    raise ValueError(f"Unsupported population for LD: {pop}")


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

    # First, prefer exact single-gene tag for backward compatibility.
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

    fig, ax = plt.subplots(1, 2, figsize=(12, 5))

    for pop, color in [("EUR", "tab:blue"), ("AFR", "tab:orange")]:
        ax[0].scatter(
            df_region["pos"],
            -np.log10(np.maximum(df_region[f"p_{pop}"], 1e-300)),
            s=12,
            alpha=0.55,
            color=color,
            label=pop,
        )
    ax[0].set_xlabel("Position")
    ax[0].set_ylabel("-log10(p)")
    ax[0].set_title("Regional association")
    ax[0].legend()

    lead = df_region.loc[df_region["p_meta"].idxmin(), "snp"] if "p_meta" in df_region.columns else df_region.loc[df_region["p_EUR"].idxmin(), "snp"]

    for (p1, p2), a in zip([("EUR", "AFR")], [ax[1]]):
        r1 = ld_by_pop[p1]
        r2 = ld_by_pop[p2]
        if r1.shape[0] == 0 or r2.shape[0] == 0 or r1.shape != r2.shape:
            a.set_title(f"{p1} vs {p2} LD mismatch (insufficient shared SNPs)")
            continue
        if lead not in df_region["snp"].values:
            a.set_title(f"{p1} vs {p2} LD mismatch (lead SNP missing)")
            continue
        lead_idx = int(df_region.index[df_region["snp"] == lead][0])
        x = np.square(r1[:, lead_idx])
        y = np.square(r2[:, lead_idx])
        c = -np.log10(np.maximum(df_region["p_meta"], 1e-300)) if "p_meta" in df_region.columns else -np.log10(np.maximum(df_region["p_EUR"], 1e-300))
        a.scatter(x, y, c=c, s=20, alpha=0.7, cmap="viridis")
        a.set_xlabel(f"r2 {p1} with lead")
        a.set_ylabel(f"r2 {p2} with lead")
        a.set_title(f"LD mismatch: {p1} vs {p2}")

    fig.tight_layout()
    fig.savefig(f"{output_prefix}_regional_diagnostics.png", dpi=160)
    plt.close(fig)


def write_finemap_files(df, ld_matrix, prefix, pop):
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


def main():
    parser = argparse.ArgumentParser(description="Prepare EUR/AFR LD matrices and diagnostics for fine-mapping")
    parser.add_argument("--gene", required=True, help="Gene name to fine-map (e.g. BACH2)")
    parser.add_argument('--flank', default='500000',
                        help='Flank size in bp (default: 500000). Supports: 500000, 1e6, 1Mb, 500kb, etc.')
    parser.add_argument("--input-files", "-i", nargs="+", required=True, help="Input GWAS summary statistics files (one per population, can be gzipped) and the corresponding ancestry group, separated by comma \"filename,ANC\"")
    parser.add_argument("--ref-genome", choices=["GRCh37", "GRCh38"], default="GRCh37", help="Reference genome version (default: GRCh37)")
    parser.add_argument("--output-prefix", required=True, help="Prefix for output files (e.g. BACH2_finemap)")
    parser.add_argument(
        "--local-ld-dir", required=True,
        help=(
            "Directory containing pre-fetched Pan-UKBB LD files produced by step2a_extract_hail_LDs.py "
            "(e.g. UKBB.EUR.<chrom>_<start>_<end>.ld.npy). "
            "Required: this script only loads EUR/AFR LD from local cache."
        ),
    )
    parser.add_argument(
        "--skip-global-diagnostics", action="store_true",
        help="Skip generation of global diagnostic plots (default: generate plots)"
    )
    args = parser.parse_args()
    args.flank_raw = args.flank          # original string (e.g. "1Mb", "500000")
    args.flank = parse_flank_size(args.flank_raw)

    t0 = time.time()
    outdir = os.path.dirname(args.output_prefix)
    if outdir:
        os.makedirs(outdir, exist_ok=True)

    gene_chr, gene_start, gene_end = get_gene_coordinates(
        args.gene,
        gene_ref_path=b37_GENE_REF if args.ref_genome == "GRCh37" else b38_GENE_REF,
    )
    region_start = max(1, int(gene_start - args.flank))
    region_end = int(gene_end + args.flank)
    print(f"Using region chr{gene_chr}:{region_start}-{region_end} (gene +/- flank)")

    files = parse_input_files(args.input_files)
    raw = {}
    for pop in ["EUR", "AFR"]:
        print(f"Reading {pop} GWAS: {files[pop]}")
        raw_df = read_table(files[pop])
        print(f"{pop} variants before standardization: {len(raw_df)}")
        print(raw_df.head(2))
        raw[pop] = standardize_gwas(raw_df, pop)
        print(f"{pop}: {len(raw[pop])} variants after standardization")

    # Per-ancestry region tables (before cross-ancestry harmonization).
    target_chr = normalize_chr(gene_chr)
    region_by_pop = {}
    for pop in ["EUR", "AFR"]:
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
    if len(set(resolved_tags.values())) == 1:
        print("  Note: EUR and AFR use the same tag.")
    else:
        print("  Note: EUR and AFR use different tags.")

    ld_raw = {}
    for pop in ["EUR", "AFR"]:
        tag = resolved_tags[pop]
        mat, var = load_ld_subset(
            pop,
            chrom=target_chr,
            start=region_start,
            end=region_end,
            tag=tag,
            ref_genome=args.ref_genome,
            local_ld_dir=local_ld_dir,
        )
        ld_raw[pop] = (mat, var)
        print(f"{pop} LD variants loaded: {len(var)}")

    # Check for pre-existing harmonized file when skipping diagnostics
    harmonized_path = f"{args.output_prefix}.harmonized.tsv.gz"
    if args.skip_global_diagnostics and os.path.exists(harmonized_path):
        print(f"Loading pre-existing harmonized file: {harmonized_path}")
        harm = read_table(harmonized_path)
        # Ensure proper data types for the loaded harmonized data
        for col in harm.columns:
            if col in ["chr"]:
                harm[col] = harm[col].astype(str).map(normalize_chr)
            elif col in ["pos", "n_EUR", "n_AFR"]:
                harm[col] = pd.to_numeric(harm[col], errors="coerce")
            elif col in ["maf_EUR", "maf_AFR", "beta_EUR", "beta_AFR", "se_EUR", "se_AFR", 
                        "p_EUR", "p_AFR", "beta_meta", "z_meta", "p_meta", "Q", "I2", "tau2", "HetP"]:
                harm[col] = pd.to_numeric(harm[col], errors="coerce")
        print(f"Loaded {len(harm)} harmonized variants from file")
    else:
        harm = harmonize_two_way(raw)
        harm = compute_heterogeneity(harm)
        with np.errstate(divide="ignore"):
            harm["z_meta"] = harm["beta_meta"] / np.sqrt(
                1.0 / (1.0 / np.square(harm["se_EUR"]) + 1.0 / np.square(harm["se_AFR"]))
            )
        harm["p_meta"] = 2.0 * stats.norm.sf(np.abs(harm["z_meta"]))
        harm = harm.sort_values(["chr", "pos"], key=lambda s: s.map(chr_sort_key) if s.name == "chr" else s).reset_index(drop=True)
        harm.to_csv(harmonized_path, sep="\t", index=False, compression="gzip")
        print(f"Harmonized variants across EUR/AFR: {len(harm)}")
    

    if not args.skip_global_diagnostics:
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
    print("Regional retention after 2-way harmonization:")
    for pop in ["EUR", "AFR"]:
        pre_n = len(region_by_pop[pop])
        shared_n = int(region_by_pop[pop]["snp"].isin(region["snp"]).sum())
        print(f"  {pop}: regional pre-harmonization={pre_n}, retained in 2-way region={shared_n}, dropped={pre_n - shared_n}")

    if len(region) == 0:
        raise RuntimeError("No harmonized variants found in target region.")

    aligned = {}
    keep_sets = []
    for pop in ["EUR", "AFR"]:
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
    for pop, keep_set in zip(["EUR", "AFR"], keep_sets):
        overlap_n = len(common_snps & keep_set)
        print(f"[common] overlap with {pop} LD-aligned SNPs: {overlap_n} / {len(common_snps)}")
    for s in keep_sets:
        common_snps &= s
    print(f"[common] final SNPs shared by harmonized region and both LD panels: {len(common_snps)}")
    if len(common_snps) == 0:
        raise RuntimeError("No common SNPs across region GWAS and both LD panels after allele harmonization.")

    region_common = region[region["snp"].isin(common_snps)].copy()
    region_common = region_common.sort_values(["chr", "pos"], key=lambda s: s.map(chr_sort_key) if s.name == "chr" else s).reset_index(drop=True)

    ld_final = {}
    for pop in ["EUR", "AFR"]:
        mat, keep = aligned[pop]
        idx_map = {snp: i for i, snp in enumerate(keep["snp"].tolist())}
        ix = np.array([idx_map[s] for s in region_common["snp"]], dtype=int)
        ld_final[pop] = mat[np.ix_(ix, ix)]

    regional_diagnostics(region_common, ld_final, args.output_prefix)

    region_common.to_csv(f"{args.output_prefix}.region_common.tsv.gz", sep="\t", index=False, compression="gzip")
    for pop in ["EUR", "AFR"]:
        write_finemap_files(region_common, ld_final[pop], args.output_prefix, pop)

    dt = time.time() - t0
    print(f"Done in {dt:.1f} seconds")


if __name__ == "__main__":
    main()
