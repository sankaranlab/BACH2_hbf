#!/usr/bin/env python3
"""
This script integrates all summary statistics from FEMA, MAMA, and 11 cohorts,
  and make supplementary table 3 for the manuscript.

FEMA header: ("MarkerName"--> rsID)
- "Chromosome      Position        MarkerName      Allele1 Allele2 Freq1   FreqSE  MinFreq MaxFreq Weight      Zscore  P-value Direction"
- "CHR      GENPOS        MarkerName      Allele1 Allele2 Freq1   FreqSE  MinFreq MaxFreq Weight      Zscore  P-value Direction"

MAMA header: ("SNP" is not rsID) (FREQ and BETA is for A1)
- "SNP     CHR     BP      A1      A2      FREQ    BETA    SE      Z       P       N_EFF   N_ORIG"

Cohorts header:
- BIOS sumstats:
    "SNP     CHR     BP      A1      A2      ID      AC_ALLELE2      AF_ALLELE2      IMPUTATIONINFO  N  BETA     SE      TSTAT   P       VART    VARTSTAR        Z       IMPUTATION_gen_build"
- GTEx sumstats:
    "SNP     CHR     BP      A1      A2      VARIANT_ID      P       BETA    SE      AF_ALLELE2      alt_alleles Z       N"
- Interval sumstats:
    "SNP     CHR     BP      A1      A2      AF_ALLELE2      INFO    BETA    SE      P       N       Z  IMPUTATION_gen_build"
- Sardinia sumstats:
    "SNP     CHR     BP      A1      A2      ID      BETA    SE      AF_ALLELE2      N       Z       P  IMPUTATION_gen_build"
- StJude sumstats:
    "SNP     CHR     BP      A1      A2      N       AC      CALLRATE        GENOCNT AF_ALLELE2      SIGNED_SUMSTAT      P       BETA    SE      Z"
- Sweden sumstats:
    "SNP     CHR     BP      A1      A2      AC_ALLELE2      AF_ALLELE2      IMPUTATIONINFO  N       BETASE      TSTAT   P       VART    VARTSTAR        ID      Z"
- Tanzania sumstats:
    "SNP     CHR     BP      VARIANT_ID      A1      A2      AF_ALLELE2      N       BETA    SE      LOG10P      hg38_ref        hg38_alt        P       Z"
- Thai sumstats:
    "SNP     CHR     BP      A1      A2      ID      AF_ALLELE2      IMPUTATIONINFO  BETA    P       N  SE       Z       IMPUTATION_gen_build"
- TOPMed sumstats:
    "SNP     CHR     BP      A1      A2      AC_ALLELE2      AF_ALLELE2      IMPUTATIONINFO  N       BETASE      TSTAT   P       VART    VARTSTAR        Z"
Post-QC sumstats header:
cohort_sumstats/BIOS_LL_chrALL_MAF0p1pct_hg38_info6.clean.tsv.gz
snp	chr	pos	a1	a2	id	ea_count	ea_freq	info	n	beta	se	TSTAT	pVART	VARTSTAR	z	impute	maf	mac

cohort_sumstats/BIOS_LLS_660Q_chrALL_MAF0p1pct_hg38_info6.clean.tsv.gz
snp	chr	pos	a1	a2	id	ea_count	ea_freq	info	n	beta	se	TSTAT	pVART	VARTSTAR	z	impute	maf	mac

cohort_sumstats/BIOS_RS_chrALL_MAF0p1pct_hg38_info6.clean.tsv.gz
snp	chr	pos	a1	a2	id	ea_count	ea_freq	info	n	beta	se	TSTAT	pVART	VARTSTAR	z	impute	maf	mac

cohort_sumstats/GTEx_chrALL_MAF0p1pct_hg38_info6.clean.tsv.gz
snp	chr	pos	a1	a2	VARIANT_ID	p	beta	se	ea_freq	alt_alleles	z	nmaf

cohort_sumstats/Interval_chrALL_MAF0p1pct_hg38_info6_pass-z-check.clean.tsv.gz
snp	chr	pos	a1	a2	ea_freq	info	beta	se	p	n	z	impute	maf

cohort_sumstats/Sardinia_chrALL_MAF0p1pct_hg38_info6.clean.tsv.gz
snp	chr	pos	a1	a2	id	beta	se	ea_freq	n	z	p	impute	maf

cohort_sumstats/StJude_chrALL_MAF0p1pct_hg38_info6_se-cap-10.clean.tsv.gz
snp	chr	pos	a1	a2	n	ea_count	CALLRATE	GENOCNT	ea_freq	SIGNED_SUMSTAT	pbeta	se	z	maf	mac

cohort_sumstats/Sweden_chrALL_MAF0p1pct_hg38_info6.clean.tsv.gz
snp	chr	pos	a1	a2	ea_count	ea_freq	info	n	beta	se	TSTAT	p	VART	VARTSTAR	id	z	maf	mac

cohort_sumstats/Tanzania_chrALL_MAF0p1pct_hg38_info6.clean.tsv.gz
snp	chr	pos	VARIANT_ID	a1	a2	ea_freq	n	beta	se	LOG10P	hg38_ref	hg38_alt	p	z	maf

cohort_sumstats/Thai_chrALL_MAF0p1pct_hg38_info6.clean.tsv.gz
snp	chr	pos	a1	a2	id	ea_freq	info	beta	p	n	se	z	impute	maf

cohort_sumstats/TOPMed_chrALL_MAF0p1pct_hg38_info6.clean.tsv.gz
snp	chr	pos	a1	a2	ea_count	ea_freq	info	n	beta	se	TSTAT	p	VART	VARTSTAR	z	maf	mac

    

Output header:
"MarkerName	Chr	Position	A1	A2	AF2_META	BETA_META	SE_META	P_META	EUR_MAMA_AF2	EUR_MAMA_BETA	EUR_MAMA_SE	EUR_MAMA_P	AFR_MAMA_AF2	AFR_MAMA_BETA	AFR_MAMA_SE	AFR_MAMA_P	THAI_MAMA_AF2	THAI_MAMA_BETA	THAI_MAMA_SE	THAI_MAMA_P	INTERVAL_AF2	INTERVAL_BETA	INTERVAL_SE	INTERVAL_P	SARDINIA_AF2	SARDINIA_BETA	SARDINIA_SE	SARDINIA_P	SWEDEN_AF2	SWEDEN_BETA	SWEDEN_SE	SWEDEN_P	GTEx_AF2	GTEx_BETA	GTEx_SE	GTEx_P	BIOS_LL_AF2	BIOS_LL_BETA	BIOS_LL_SE	BIOS_LL_P	BIOS_LLS_AF2	BIOS_LLS_BETA	BIOS_LLS_SE	BIOS_LLS_P	BIOS_RS_AF2	BIOS_RS_BETA	BIOS_RS_SE	BIOS_RS_P	TOPMED_AF2	TOPMED_BETA	TOPMED_SE	TOPMED_P	StJude_AF2	StJude_BETA	StJude_SE	StJude_P	TANZANIA_AF2	TANZANIA_BETA	TANZANIA_SE	TANZANIA_P	THAI_AF2	THAI_BETA	THAI_SE	THAI_P"
"""
import sys
import time
import os
import pandas as pd
import math
import re


# set working directory
LABSHARE = "/lab-share/Hem-Sankaran-e2/Public"
WDIR = f"{LABSHARE}/projects/xhcheng/HbF/Hbf_BACH2_wip/metaGWAS_analysis"


def normalize_chromosome(chrom):
    if pd.isna(chrom):
        return pd.NA

    chrom_str = str(chrom).strip()
    if chrom_str == "" or chrom_str.lower() in {"na", "nan", "none"}:
        return pd.NA

    return re.sub(r"^chr", "", chrom_str, flags=re.IGNORECASE).strip().upper()


def build_fema_position_issue(chrom, bp):
    issues = []
    if pd.isna(chrom):
        issues.append("missing_chr")
    elif str(chrom).upper() == "Y":
        issues.append("chrY")

    if pd.isna(bp):
        issues.append("missing_bp")
 
    return ";".join(issues) if issues else pd.NA


def compute_neglog10p_from_string(p_val):
    """
    Parse p-value from string (including scientific notation) and return -log10(P).
    Handles ultra-small values like 1e-400 without floating-point underflow.
    
    Args:
        p_val: p-value as string or float
    
    Returns:
        float (-log10(P)) or None if invalid
    """
    if pd.isna(p_val):
        return None
    
    p_clean = str(p_val).strip().lower()
    
    if p_clean in ['na', 'nan', 'inf', '-inf', '']:
        return None
    
    # Try scientific notation: mantissa e exponent
    match = re.match(r'^([+-]?[0-9]*\.?[0-9]+)(?:e([+-]?[0-9]+))?$', p_clean)
    if match:
        mantissa_str = match.group(1)
        exponent_str = match.group(2)
        
        try:
            mantissa = float(mantissa_str)
            exponent = float(exponent_str) if exponent_str else 0
            
            if mantissa <= 0:
                return None
            
            # Compute -log10(p) = -(log10(mantissa) + exponent)
            return -(math.log10(mantissa) + exponent)
        except (ValueError, OverflowError):
            return None
    
    # Try regular decimal
    try:
        p_num = float(p_clean)
        if p_num <= 0 or math.isinf(p_num):
            return None
        return -math.log10(p_num)
    except (ValueError, OverflowError):
        return None


def swap_output_allele_orientation(df):
    """
    Swap A1/A2 in the final output table and flip effect direction for every
    harmonized dataset so all reported statistics remain aligned to the new A2.
    """
    swapped = df.copy()

    if {"A1", "A2"}.issubset(swapped.columns):
        original_a1 = swapped["A1"].copy()
        swapped["A1"] = swapped["A2"]
        swapped["A2"] = original_a1

    af_cols = [col for col in swapped.columns if col.endswith("_AF2")]
    beta_cols = [col for col in swapped.columns if col.endswith("_BETA")]
    z_cols = [col for col in swapped.columns if col.endswith("_Z") or col == "Z_META"]

    for col in af_cols:
        swapped[col] = pd.to_numeric(swapped[col], errors="coerce").rsub(1)

    for col in beta_cols + z_cols:
        swapped[col] = -pd.to_numeric(swapped[col], errors="coerce")

    return swapped


FILE_PATHS = {
    # "FEMA": "cohort_sumstats/FEMA_CloudTanz_METALOUT_inv_1.tbl",
    "FEMA": "cohort_sumstats/METAL_hbf_inv_ALL_info6_mac40_SE_gcOn_hg38.tsv.gz",
    # "MAMA_skeleton": "cohort_sumstats/UA_HBF_MAMA_{pop}.res", # pop can be EUR, AFR, THAI
    # "MAMA_skeleton": "cohort_sumstats/XC_newFema_16865_{pop}_HBF.res", # pop can be EUR, AFR, THAI
    "MAMA_skeleton": "MAMA/XC_SNPaligned_16865_{pop}_HBF.res", # pop can be EUR, AFR, THAI
    "BIOS": "cohort_sumstats/BIOS_{panel}_chrALL_MAF0p1pct_hg38_info6.clean.tsv.gz", # panel can be LL, LLS_660Q, RS
    "GTEx": "cohort_sumstats/GTEx_chrALL_MAF0p1pct_hg38_info6.clean.tsv.gz", 
    "Sardinia": "cohort_sumstats/Sardinia_chrALL_MAF0p1pct_hg38_info6.clean.tsv.gz",
    "Sweden": "cohort_sumstats/Sweden_chrALL_MAF0p1pct_hg38_info6.clean.tsv.gz",
    "Interval": "cohort_sumstats/Interval_chrALL_MAF0p1pct_hg38_info6_pass-z-check.clean.tsv.gz",
    "Tanzania": "cohort_sumstats/Tanzania_chrALL_MAF0p1pct_hg38_info6.clean.tsv.gz",
    "StJude": "cohort_sumstats/StJude_chrALL_MAF0p1pct_hg38_info6_se-cap-10.clean.tsv.gz",
    "TOPMed": "cohort_sumstats/TOPMed_chrALL_MAF0p1pct_hg38_info6.clean.tsv.gz",
    "Thai": "cohort_sumstats/Thai_chrALL_MAF0p1pct_hg38_info6.clean.tsv.gz",
}


def read_and_standardize(filepath, source_name, cohort=None, literal_p_as_string=False):
    """
    Read summary stats file and extract SNP, AF_ALLELE2, BETA, SE, P with standardized names.
    Returns DataFrame with columns: SNP, CHR, BP, A1, A2, AF, BETA, SE, P
    """
    full_path = os.path.join(WDIR, filepath)
    if literal_p_as_string:
        # Keep p-value text as string to avoid float underflow for very small values (e.g. <1e-300).
        pvalue_converters = {key: lambda x: str(x).strip() for key in ["P", "P-value"]}
        df = pd.read_csv(full_path, sep="\t", low_memory=False, converters=pvalue_converters)
    else:
        df = pd.read_csv(full_path, sep="\t", low_memory=False)

    print(f"  Loaded {source_name} from {filepath}: {df.shape[0]} rows x {df.shape[1]} cols")
    
    # Standardize key columns based on source
    if source_name == "FEMA":
        df = df.rename(columns={
            "MarkerName": "SNP",
            "Chromosome": "CHR",
            "Position": "BP",
            "Allele2": "A1",
            "Allele1": "A2",
            "Freq1": "AF",
            "Zscore": "Z",
            "P-value": "P",
            "Weight": "N"
        })
        if "GENPOS" in df.columns:
            df = df.rename(columns={"GENPOS": "BP"})
    elif source_name == "MAMA":
        df = df.rename(columns={
            "FREQ": "AF",
            "A1": "A2",
            "A2": "A1"
            })
    else:
        # All cohorts use standard naming or need AF_ALLELE2 -> AF mapping
        if "AF_ALLELE2" in df.columns:
            df["AF"] = df["AF_ALLELE2"]
        elif "ea_freq" in df.columns: # this is probably from post-QC files
            df = df.rename(columns={
                "snp": "SNP",
                "chr": "CHR",
                "pos": "BP",
                "a1": "A1",
                "a2": "A2",
                "id": "ID",
                "info": "INFO",
                "ea_count": "AC",
                "mac": "MAC",
                "ea_freq": "AF",
                "maf": "MAF",
                "beta": "BETA",
                "se": "SE",
                "z": "Z",
                "p": "P",
                "n": "N"
            })
    
    # Extract required columns
    required_cols = ["SNP", "CHR", "BP", "A1", "A2", "AF", "BETA", "SE", "Z", "P", "N"]
    existing_cols = [col for col in required_cols if col in df.columns]
    df = df[existing_cols]

    if "SNP" in df.columns:
        df["SNP"] = df["SNP"].astype("string").str.strip()

    if "CHR" in df.columns:
        df["CHR"] = df["CHR"].apply(normalize_chromosome)

    if "BP" in df.columns:
        df["BP"] = pd.to_numeric(df["BP"], errors="coerce")

    # make sure A1 and A2 are uppercase for consistent merging
    if "A1" in df.columns:
        df["A1"] = df["A1"].astype("string").str.upper()
    if "A2" in df.columns:
        df["A2"] = df["A2"].astype("string").str.upper()
    
    # Preserve original p-value text for auditing/debugging.
    if "P" in df.columns:
        df["P_RAW"] = df["P"].astype("string").str.strip()

    # Convert P to -log10(P) with precision-safe parsing from string
    if "P" in df.columns:
        df["P"] = df["P"].apply(compute_neglog10p_from_string)
    
    # Add cohort identifier for later merging
    if cohort:
        df["COHORT"] = cohort
    
    return df


def extract_cohort_data(df, cohort_name, prefix):
    """
    Extract AF, BETA, SE, P from cohort data and rename with cohort prefix.
    Note: P is already -log10(P) from compute_neglog10p_from_string.
    """
    result = df[["SNP", "AF", "BETA", "SE", "Z", "P"]].copy()
    result = result.rename(columns={
        "AF": f"{prefix}_AF2",
        "BETA": f"{prefix}_BETA",
        "SE": f"{prefix}_SE",
        "Z": f"{prefix}_Z",
        "P": f"{prefix}_log10p"
    })
    return result


def harmonize_and_merge_by_pos(base_df, incoming_df, prefix, incompatible_outfile=None):
    """
    Merge incoming cohort into base by chromosome/position and harmonize allele direction.
    If incoming A1/A2 are swapped relative to base A1/A2, flip AF and BETA.
    """
    required_cols = ["SNP", "CHR", "BP", "A1", "A2", "AF", "BETA", "SE", "Z", "P"]
    missing = [c for c in required_cols if c not in incoming_df.columns]
    if missing:
        raise ValueError(f"Missing required columns in incoming df for {prefix}: {missing}")

    temp = incoming_df[required_cols].copy()
    temp = temp.rename(columns={
        "SNP": "_in_snp",
        "CHR": "_in_chr",
        "BP": "_in_bp",
        "A1": "_in_a1",
        "A2": "_in_a2",
        "AF": f"{prefix}_AF2",
        "BETA": f"{prefix}_BETA",
        "SE": f"{prefix}_SE",
        "Z": f"{prefix}_Z",
        "P": f"{prefix}_log10p",
    })

    # Normalize merge keys so formatting differences do not eliminate matches
    # (e.g., "chr1" vs "1", numeric vs string BP).
    base_work = base_df.copy()
    base_work["_key_chr"] = base_work["Chr"].apply(normalize_chromosome)
    base_work["_key_bp"] = pd.to_numeric(base_work["Position"], errors="coerce")

    temp["_key_chr"] = temp["_in_chr"].apply(normalize_chromosome)
    temp["_key_bp"] = pd.to_numeric(temp["_in_bp"], errors="coerce")

    merged = base_work.merge(
        temp,
        on=["_key_chr", "_key_bp"],
        how="left",
    )

    # For FEMA rows with missing or invalid original coordinates, recover CHR/BP by rsID.
    needs_position_recovery = (
        merged.get("FEMA_position_ambiguous", False).fillna(False)
        & merged["Chr"].isna()
        & merged["Position"].isna()
        & merged["MarkerName"].astype("string").str.match(r"^rs", case=False, na=False)
    )

    rs_lookup = temp[temp["_in_snp"].astype("string").str.match(r"^rs", case=False, na=False)].copy()
    rs_lookup = rs_lookup.loc[~rs_lookup["_in_snp"].duplicated(keep=False)]
    rs_lookup = rs_lookup.rename(columns={
        "_in_snp": "_rs_snp",
        "_in_chr": "_rs_chr",
        "_in_bp": "_rs_bp",
        "_in_a1": "_rs_a1",
        "_in_a2": "_rs_a2",
        f"{prefix}_AF2": f"{prefix}_AF2_rs",
        f"{prefix}_BETA": f"{prefix}_BETA_rs",
        f"{prefix}_SE": f"{prefix}_SE_rs",
        f"{prefix}_Z": f"{prefix}_Z_rs",
        f"{prefix}_log10p": f"{prefix}_log10p_rs",
    })
    rs_lookup = rs_lookup.drop(columns=["_key_chr", "_key_bp"])

    merged = merged.merge(rs_lookup, left_on="MarkerName", right_on="_rs_snp", how="left")

    recovered_by_rsid = needs_position_recovery & merged["_rs_chr"].notna() & merged["_rs_bp"].notna()
    if recovered_by_rsid.any():
        replacement_cols = [
            ("_in_chr", "_rs_chr"),
            ("_in_bp", "_rs_bp"),
            ("_in_a1", "_rs_a1"),
            ("_in_a2", "_rs_a2"),
            (f"{prefix}_AF2", f"{prefix}_AF2_rs"),
            (f"{prefix}_BETA", f"{prefix}_BETA_rs"),
            (f"{prefix}_SE", f"{prefix}_SE_rs"),
            (f"{prefix}_Z", f"{prefix}_Z_rs"),
            (f"{prefix}_log10p", f"{prefix}_log10p_rs"),
        ]
        for target_col, source_col in replacement_cols:
            merged.loc[recovered_by_rsid, target_col] = merged.loc[recovered_by_rsid, source_col]

        merged.loc[recovered_by_rsid, "Chr"] = merged.loc[recovered_by_rsid, "_rs_chr"]
        merged.loc[recovered_by_rsid, "Position"] = merged.loc[recovered_by_rsid, "_rs_bp"]
        if "Position_source" in merged.columns:
            merged.loc[recovered_by_rsid, "Position_source"] = prefix

    same = (merged["A1"] == merged["_in_a1"]) & (merged["A2"] == merged["_in_a2"])
    swapped = (merged["A1"] == merged["_in_a2"]) & (merged["A2"] == merged["_in_a1"])
    matched = merged["_in_chr"].notna() & merged["_in_bp"].notna()

    # Flag allele mismatches among rows that matched by Chr/BP.
    incompatible = matched & ~(same | swapped)
    if incompatible.any() and incompatible_outfile:
        incompatible_rows = merged.loc[
            incompatible,
            ["MarkerName", "Chr", "Position", "A1", "A2", "_in_a1", "_in_a2"],
        ].copy()
        incompatible_rows.insert(0, "COHORT", prefix)
        incompatible_rows = incompatible_rows.rename(columns={"_in_a1": "incoming_A1", "_in_a2": "incoming_A2"})
        incompatible_rows.to_csv(
            incompatible_outfile,
            sep="\t",
            index=False,
            mode="a",
            header=not os.path.exists(incompatible_outfile),
        )

    # Keep only identical or swapped allele orientation; otherwise set to missing.
    merged.loc[
        ~(same | swapped),
        [f"{prefix}_AF2", f"{prefix}_BETA", f"{prefix}_SE", f"{prefix}_Z", f"{prefix}_log10p"],
    ] = pd.NA

    # Flip effect direction/frequency when allele orientation is reversed.
    merged.loc[swapped, f"{prefix}_AF2"] = 1 - merged.loc[swapped, f"{prefix}_AF2"]
    merged.loc[swapped, f"{prefix}_BETA"] = -merged.loc[swapped, f"{prefix}_BETA"]
    merged.loc[swapped, f"{prefix}_Z"] = -merged.loc[swapped, f"{prefix}_Z"]

    print(
        f"    {prefix}: matched={matched.sum()}, same={same[matched].sum()}, "
        f"swapped={swapped[matched].sum()}, incompatible={incompatible.sum()}, "
        f"recovered_by_rsid={recovered_by_rsid.sum()}, unmatched={(~matched).sum()}"
    )

    merged = merged.drop(columns=[
        "_in_snp", "_in_chr", "_in_bp", "_in_a1", "_in_a2", "_key_chr", "_key_bp",
        "_rs_snp", "_rs_chr", "_rs_bp", "_rs_a1", "_rs_a2",
        f"{prefix}_AF2_rs", f"{prefix}_BETA_rs", f"{prefix}_SE_rs", f"{prefix}_Z_rs", f"{prefix}_log10p_rs",
    ])
    return merged


def main():
    t0 = time.time()
    def ts():
        return f"+{time.time() - t0:.1f}s"
    
    input_suffix = sys.argv[1] if len(sys.argv) > 1 else "input"

    incompatible_outfile = f"supplementary_table_3_{input_suffix}_incompatible.tsv"
    ambiguous_position_outfile = f"supplementary_table_3_{input_suffix}_ambiguous_fema_positions.tsv"
    if os.path.exists(incompatible_outfile):
        os.remove(incompatible_outfile)
    if os.path.exists(ambiguous_position_outfile):
        os.remove(ambiguous_position_outfile)
    print(f"[{ts()}] Incompatible allele rows will be written to {incompatible_outfile}")

    print("[1/6] Loading FEMA meta-analysis...")
    # Read FEMA as base
    fema = read_and_standardize(FILE_PATHS["FEMA"], "FEMA", 
                                literal_p_as_string=True)
    fema = fema.rename(columns={"SNP": "MarkerName", "CHR": "Chr", "BP": "Position", "AF": "AF2_META", "P": "P_META_log10p", "N": "N_META"})

    # Print all FEMA rows where P < 1e-300 (equivalent to -log10(P) > 300).
    if "P_META_log10p" in fema.columns:
        very_small_p_mask = fema["P_META_log10p"] > 300
        n_very_small = int(very_small_p_mask.sum())
        print(f"  [{ts()}] FEMA rows with P < 1e-300: {n_very_small}")
        if n_very_small > 0:
            cols_to_show = [
                c for c in ["MarkerName", "Chr", "Position", "A1", "A2", "P_RAW", "P_META_log10p"]
                if c in fema.columns
            ]
            print("\t".join(cols_to_show))
            for _, row in fema.loc[very_small_p_mask, cols_to_show].iterrows():
                print("\t".join("" if pd.isna(v) else str(v) for v in row.values))

    # Extract FEMA meta-analysis results
    fema_data = fema[["MarkerName", "Chr", "Position", "A1", "A2"]].copy()
    fema_data["AF2_META"] = fema["AF2_META"]
    # fema_data["BETA_META"] = fema["BETA"]
    # fema_data["SE_META"] = fema["SE"]
    fema_data["P_META_log10p"] = fema["P_META_log10p"]
    fema_data["N_META"] = fema["N_META"]

    ambiguous_fema_mask = (
        fema_data["Chr"].isna()
        | (fema_data["Chr"] == "Y")
        | fema_data["Position"].isna()
    )
    fema_data["FEMA_original_Chr"] = fema_data["Chr"]
    fema_data["FEMA_original_Position"] = fema_data["Position"]
    fema_data["FEMA_position_ambiguous"] = ambiguous_fema_mask
    fema_data["FEMA_position_issue"] = [
        build_fema_position_issue(chrom, bp)
        for chrom, bp in zip(fema_data["FEMA_original_Chr"], fema_data["FEMA_original_Position"])
    ]
    fema_data["Position_source"] = "FEMA"
    fema_data.loc[ambiguous_fema_mask, ["Chr", "Position"]] = pd.NA

    print(f"  [{ts()}] FEMA loaded: {fema_data.shape[0]} rows x {fema_data.shape[1]} cols")
    print(f"  [{ts()}] FEMA rows with ambiguous coordinates: {ambiguous_fema_mask.sum()}")
    del fema

    # Start merging from FEMA base and harmonize incoming effects by Chr/Position + A1/A2.
    result = fema_data.copy()
    del fema_data

    print(f"[{ts()}] [2/6] Loading and merging BIOS panels...")
    for panel in ["LL", "LLS_660Q", "RS"]:
        filepath = FILE_PATHS["BIOS"].format(panel=panel)
        bios_df = read_and_standardize(filepath, "BIOS")
        print(f"  [{ts()}] Loaded BIOS {panel}: {bios_df.shape[0]} rows x {bios_df.shape[1]} cols")
        panel_name = "LLS" if panel == "LLS_660Q" else panel
        result = harmonize_and_merge_by_pos(result, bios_df, f"BIOS_{panel_name}", incompatible_outfile)
        print(f"  [{ts()}] Merged BIOS {panel}")
        del bios_df

    print(f"[{ts()}] [3/6] Loading and merging remaining cohorts...")
    cohort_mappings = {
        "GTEx": (FILE_PATHS["GTEx"], "GTEx"),
        "Sardinia": (FILE_PATHS["Sardinia"], "SARDINIA"),
        "Sweden": (FILE_PATHS["Sweden"], "SWEDEN"),
        "Interval": (FILE_PATHS["Interval"], "INTERVAL"),
        "Tanzania": (FILE_PATHS["Tanzania"], "TANZANIA"),
        "StJude": (FILE_PATHS["StJude"], "StJude"),
        "TOPMed": (FILE_PATHS["TOPMed"], "TOPMED"),
        "Thai": (FILE_PATHS["Thai"], "THAI"),
    }
    for cohort_key, (filepath, prefix) in cohort_mappings.items():
        cohort_df = read_and_standardize(filepath, cohort_key)
        print(f"  [{ts()}] Loaded {cohort_key}: {cohort_df.shape[0]} rows x {cohort_df.shape[1]} cols")
        result = harmonize_and_merge_by_pos(result, cohort_df, prefix, incompatible_outfile)
        print(f"  [{ts()}] Merged {cohort_key}")
        del cohort_df

    # load MAMA last because they don't have rsIDs
    print(f"[{ts()}] [4/6] Loading and merging MAMA cohorts...")
    for pop in ["EUR", "AFR", "THAI"]:
        filepath = FILE_PATHS["MAMA_skeleton"].format(pop=pop)
        mama_df = read_and_standardize(filepath, "MAMA")
        print(f"  [{ts()}] Loaded MAMA {pop}: {mama_df.shape[0]} rows x {mama_df.shape[1]} cols")
        result = harmonize_and_merge_by_pos(result, mama_df, f"{pop}_MAMA", incompatible_outfile)
        print(f"  [{ts()}] Merged MAMA {pop}")
        del mama_df

    # After all merges, output the FEMA rows with ambiguous original coordinates for transparency, including any recovered CHR/BP from rsID lookup.
    ambiguous_records = result[result["FEMA_position_ambiguous"].fillna(False)].copy()
    ambiguous_record_cols = [
        "MarkerName", "FEMA_original_Chr", "FEMA_original_Position", "FEMA_position_issue",
        "Chr", "Position", "Position_source", "A1", "A2",
    ]
    ambiguous_record_cols += [
        col for col in ambiguous_records.columns
        if col not in ambiguous_record_cols and col not in {"FEMA_position_ambiguous"}
    ]
    ambiguous_records = ambiguous_records[ambiguous_record_cols]
    ambiguous_records.to_csv(ambiguous_position_outfile, sep="\t", index=False)
    recovered_ambiguous_n = (
        ambiguous_records["Chr"].notna() & ambiguous_records["Position"].notna()
    ).sum()
    print(
        f"[{ts()}] Wrote ambiguous FEMA coordinate audit: {len(ambiguous_records)} rows, "
        f"{recovered_ambiguous_n} with recovered CHR/BP"
    )

    print(f"[{ts()}] [5/6] Post-merge filtering and summary...")
    
    # Reorder columns to match output format (with -log10P values)
    output_cols = [
        "MarkerName", "Chr", "Position", "A1", "A2", 
        "AF2_META", "Z_META", "P_META_log10p", "N_META",
        "EUR_MAMA_AF2", "EUR_MAMA_BETA", "EUR_MAMA_SE", "EUR_MAMA_log10p",
        "AFR_MAMA_AF2", "AFR_MAMA_BETA", "AFR_MAMA_SE", "AFR_MAMA_log10p",
        "THAI_MAMA_AF2", "THAI_MAMA_BETA", "THAI_MAMA_SE", "THAI_MAMA_log10p",
        "INTERVAL_AF2", "INTERVAL_BETA", "INTERVAL_SE", "INTERVAL_log10p",
        "SARDINIA_AF2", "SARDINIA_BETA", "SARDINIA_SE", "SARDINIA_log10p",
        "SWEDEN_AF2", "SWEDEN_BETA", "SWEDEN_SE", "SWEDEN_log10p",
        "GTEx_AF2", "GTEx_BETA", "GTEx_SE", "GTEx_log10p",
        "BIOS_LL_AF2", "BIOS_LL_BETA", "BIOS_LL_SE", "BIOS_LL_log10p",
        "BIOS_LLS_AF2", "BIOS_LLS_BETA", "BIOS_LLS_SE", "BIOS_LLS_log10p",
        "BIOS_RS_AF2", "BIOS_RS_BETA", "BIOS_RS_SE", "BIOS_RS_log10p",
        "TOPMED_AF2", "TOPMED_BETA", "TOPMED_SE", "TOPMED_log10p",
        "StJude_AF2", "StJude_BETA", "StJude_SE", "StJude_log10p",
        "TANZANIA_AF2", "TANZANIA_BETA", "TANZANIA_SE", "TANZANIA_log10p",
        "THAI_AF2", "THAI_BETA", "THAI_SE", "THAI_log10p",
    ]
    
    result = result[[col for col in output_cols if col in result.columns]]

    # filter to only keep SNPs that are significant in FEMA meta-analysis (P_META_log10p > -log10(5e-8))
    sig_threshold_log10p = -math.log10(5e-8)  # ~7.301
    suggested_threshold_log10p = -math.log10(1e-6)  # 6
    before_filter_n = len(result)
    # result = result[result["P_META_log10p"] > sig_threshold_log10p]
    result = result[result["P_META_log10p"] > suggested_threshold_log10p]
    # print(f"  [{ts()}] Filtered significant SNPs (log10p > {sig_threshold_log10p:.2f}): {before_filter_n} -> {len(result)}")
    print(f"  [{ts()}] Filtered significant SNPs (log10p > {suggested_threshold_log10p:.2f}): {before_filter_n} -> {len(result)}")

    # Only keep SNPs (A1 and A2 are single nucleotides). Write the rest to incompatible file for transparency.
    before_filter_n = len(result)
    non_snp = ~result["A1"].isin(["A", "C", "G", "T"]) | ~result["A2"].isin(["A", "C", "G", "T"])
    if non_snp.any():
        non_snp_rows = result.loc[non_snp, ["MarkerName", "Chr", "Position", "A1", "A2"]].copy()
        non_snp_rows.insert(0, "COHORT", "Non-SNP Alleles")
        non_snp_rows.to_csv(
            incompatible_outfile,
            sep="\t",
            index=False,
            mode="a",
            header=not os.path.exists(incompatible_outfile),
        )
    result = result[~non_snp]
    print(f"  [{ts()}] Filtered SNPs with non-single-nucleotide alleles: {before_filter_n} -> {len(result)}")

    # Check unique values in Chr; if non-integers exist, print them + the number of harmonized SNPs on those chromosomes for transparency.
    unique_chromosomes = result["Chr"].dropna().unique()
    non_integer_chromosomes = [chrom for chrom in unique_chromosomes if not str(chrom).isdigit()]
    if non_integer_chromosomes:
        print(f"  [{ts()}] Warning: Non-integer chromosome values found: {non_integer_chromosomes}")
        for chrom in non_integer_chromosomes:
            count = (result["Chr"] == chrom).sum()
            print(f"    Chromosome '{chrom}': {count} SNPs")
    # filter out any non-integer chromosomes 
    before_filter_n = len(result)
    result = result[result["Chr"].apply(lambda x: str(x).isdigit())]
    print(f"  [{ts()}] Filtered non-integer chromosomes: {before_filter_n} -> {len(result)}")
    # convert Chr to integer type after filtering
    result["Chr"] = result["Chr"].astype(int)
    # Sort by chromosome and position.
    result = result.sort_values(by=["Chr", "Position"]).reset_index(drop=True)

    # Swap to the opposite reported allele across the final output table.
    result = swap_output_allele_orientation(result)


    # Save a NA-filled version of the result before filtering out rows with missing MAMA data, for transparency and potential future use.
    result.to_csv(f"supplementary_table_3_{input_suffix}_logP_wNA.tsv", sep="\t", index=False)
    print(f"[{ts()}] Done. Supplementary Table 3 (with NAs) written: {len(result)} SNPs with -log10P values")

    # remove rows with missing values in MAMA:
    before_filter_n = len(result)
    result = result.dropna(subset=["EUR_MAMA_log10p", "AFR_MAMA_log10p", "THAI_MAMA_log10p"], how="all")
    print(f"  [{ts()}] Filtered SNPs with missing MAMA data: {before_filter_n} -> {len(result)}")

    # report how many SNPs are in the final result (per chromosome and overall)
    snps_per_chromosome = result.groupby("Chr").size()
    for chrom, count in snps_per_chromosome.items():
        print(f"Chromosome {chrom}: {count} SNPs")
    print(f"Overall: {len(result)} SNPs")
    
    print(f"[{ts()}] [6/6] Writing output file...")
    # Write output with -log10P values
    result.to_csv(f"supplementary_table_3_{input_suffix}_logP.tsv", sep="\t", index=False)
    print(f"[{ts()}] Done. Supplementary Table 3 written: {len(result)} SNPs (all P values converted to -log10P)")


if __name__ == "__main__":
    main()