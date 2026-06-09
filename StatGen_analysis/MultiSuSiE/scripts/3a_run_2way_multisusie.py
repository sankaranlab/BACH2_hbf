#!/usr/bin/env python3
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from MultiSuSiE import multisusie_rss  # type: ignore[import-not-found]


DEFAULT_POPULATIONS = ["EUR", "AFR"]


def load_finemap_files(prefix, populations=DEFAULT_POPULATIONS):
    """
    Load finemap-formatted .z and .ld files for multiple populations.
    
    Args:
        prefix: Output prefix from step2b (e.g., 'BACH2_finemap')
        populations: List of population codes
    
    Returns:
        z_list: List of z-score DataFrames
        R_list: List of correlation matrices (LD)
        snp_info: DataFrame with merged SNP information
    """
    z_list = []
    R_list = []
    snp_dfs = []
    pop_sizes = []
    maf_list = []
    
    for pop in populations:
        print(f"Loading {pop}...")
        
        # Load z-score file
        z_file = f"{prefix}.{pop}.z"
        z_df = pd.read_csv(z_file, sep="\t")
        for col in ["beta", "se", "maf", "n", "position"]:
            if col in z_df.columns:
                z_df[col] = pd.to_numeric(z_df[col], errors="coerce")
        z_list.append(z_df)
        snp_dfs.append(z_df[["rsid", "chromosome", "position", "allele1", "allele2"]])
        pop_sizes.append(z_df["n"].median())  # Assuming 'n' column exists with sample size info
        maf_list.append(z_df["maf"])  # Assuming 'maf' column exists with minor allele frequency info

        if "beta" in z_df.columns and "se" in z_df.columns:
            bad_beta = int((~np.isfinite(z_df["beta"])).sum())
            bad_se = int((~np.isfinite(z_df["se"])).sum())
            zero_se = int((z_df["se"] == 0).sum())
            if bad_beta or bad_se or zero_se:
                print(f"  [{pop}] non-finite beta={bad_beta}, non-finite se={bad_se}, zero se={zero_se}")
        
        # Load LD matrix
        ld_file = f"{prefix}.{pop}.ld"
        R = np.loadtxt(ld_file)
        if R.ndim == 1:
            R = R.reshape(1, -1)
        
        # Ensure symmetry: if matrix has zero lower-diagonal, symmetrize it
        is_symmetric = np.allclose(R, R.T, rtol=1e-10, atol=1e-10)
        if not is_symmetric:
            # Check if only upper triangular is populated (lower triangle is ~0)
            lower_tri = np.tril(R, k=-1)
            upper_tri = np.triu(R, k=1)
            if np.allclose(lower_tri, 0, atol=1e-10):
                print(f"  [{pop}] LD matrix has zero lower-diagonal; symmetrizing...")
                R = np.triu(R)  # Keep upper triangle
                R = R + R.T - np.diag(np.diag(R))  # Make symmetric
            else:
                # If asymmetry is small numerical error, average them
                print(f"  [{pop}] LD matrix not symmetric (max diff={np.max(np.abs(R - R.T)):.2e}); averaging...")
                R = (R + R.T) / 2.0
        R_list.append(R)
        
        print(f"  {pop}: {len(z_df)} variants, LD shape {R.shape}")
    
    # Merge SNP info across populations (keeping common SNPs)
    snp_info = snp_dfs[0].copy()
    snp_info = snp_info.rename(columns={"rsid": "snp"})
    # print pop size
    print(f"Population sizes: {pop_sizes} for populations {populations}")
    # compute correlation matrix of cross-population effect sizes (rho) from betas of all three populations
    beta_mat = np.vstack([z_df["beta"].to_numpy(dtype=float) for z_df in z_list if "beta" in z_df.columns])
    valid_cols = np.all(np.isfinite(beta_mat), axis=0)
    beta_mat_valid = beta_mat[:, valid_cols]
    if beta_mat_valid.shape[1] < 2 or np.any(np.nanstd(beta_mat_valid, axis=1) == 0):
        print("Warning: insufficient/constant beta values for stable rho; using identity matrix")
        corr_matrix = np.eye(len(populations), dtype=float)
    else:
        corr_matrix = np.corrcoef(beta_mat_valid)
    print(f"Cross-population correlation matrix (rho):\n{corr_matrix}")
    
    return z_list, R_list, snp_info, corr_matrix, pop_sizes, maf_list


def run_multisusie(z_list, R_list, rho, pop_sizes, maf_list,L=10, max_iter=2000):
    """
    Run MultiSuSiE RSS with multiple populations.
    
    Args:
        z_list: List of z-score DataFrames
        R_list: List of correlation matrices
        rho: Cross-population correlation matrix
        pop_sizes: List of population sizes
        maf_list: List of minor allele frequency Series for each population
        L: Number of components (credible sets)
        max_iter: Maximum iterations
    
    Returns:
        result: MultiSuSiE result object
    """
    # Extract z-scores from DataFrames
    z_scores = []
    for i, z_df in enumerate(z_list):
        if "beta" in z_df.columns and "se" in z_df.columns:
            beta = z_df["beta"].to_numpy(dtype=float)
            se = z_df["se"].to_numpy(dtype=float)
            with np.errstate(divide="ignore", invalid="ignore"):
                z = beta / se
            bad = int((~np.isfinite(z)).sum())
            if bad:
                print(f"  [pop_idx={i}] replacing {bad} non-finite z-scores with 0")
                z = np.nan_to_num(z, nan=0.0, posinf=0.0, neginf=0.0)
            z_scores.append(z)
        else:
            z = pd.to_numeric(pd.Series(z_df.values[:, 0]), errors="coerce").to_numpy(dtype=float)
            bad = int((~np.isfinite(z)).sum())
            if bad:
                print(f"  [pop_idx={i}] fallback z column has {bad} non-finite values; setting to 0")
                z = np.nan_to_num(z, nan=0.0, posinf=0.0, neginf=0.0)
            z_scores.append(z)
    
    print("\nRunning MultiSuSiE…")
    result = multisusie_rss(
        R_list=R_list,
        z_list=z_scores,
        rho=rho,
        population_sizes=pop_sizes,
        maf_list=maf_list,
        L=L,
        max_iter=max_iter,
        estimate_residual_variance=True
    )
    
    return result


def parse_multisusie_result(result, snp_info, populations=DEFAULT_POPULATIONS):
    """
    Parse MultiSuSiE output once and store reusable structures.

    Args:
        result: MultiSuSiE result object
        snp_info: DataFrame with SNP information
        populations: List of population codes

    Returns:
        dict with parsed arrays/dataframes used downstream
    """
    parsed = {
        "pip": None,
        "coef": None,
        "coef_sd": None,
        "sets": None,
        "diag": {},
    }

    for key in ["pip", "coef", "coef_sd", "sets"]:
        if hasattr(result, key):
            parsed[key] = getattr(result, key)

    for key in ["converged", "niter", "lbf", "KL", "ER2", "sigma2"]:
        if hasattr(result, key):
            parsed["diag"][key] = getattr(result, key)

    if parsed["pip"] is not None:
        pip = np.asarray(parsed["pip"], dtype=float)
        print(f"\nParsed PIP for {pip.size} variants")
        print(f"  {(pip > 0.5).sum()} variants with PIP > 0.5")
        print(f"  {(pip > 0.9).sum()} variants with PIP > 0.9")
        if (pip > 0.5).any():
            print(f"    {', '.join(snp_info['snp'][pip > 0.5])}")

    coef = parsed["coef"]
    if coef is not None and coef.shape[0] != len(populations):
        print(f"Warning: coef has {coef.shape[0]} populations but expected {len(populations)}")

    if "converged" in parsed["diag"]:
        print(f"  Converged: {parsed['diag']['converged']}")
    if "niter" in parsed["diag"]:
        print(f"  Iterations: {parsed['diag']['niter']}")

    return parsed


def extract_credible_sets(parsed, snp_info, output_prefix):
    """
    Extract and save credible sets with variant details.
    
    Args:
        parsed: Parsed MultiSuSiE structures
        snp_info: DataFrame with SNP information
        output_prefix: Output prefix for output file
    """
    print("\nExtracting credible sets…")

    sets_list = parsed.get("sets")
    if sets_list is None:
        print("Warning: result.sets not available")
        return
    
    # sets[0] = indices, sets[1] = purity, sets[2] = coverage, sets[3] = passed_filter
    if len(sets_list) < 4:
        print("Warning: incomplete sets structure")
        return
    
    indices, purity, coverage, passed_filter = sets_list[:4]
    
    cs_records = []
    for i, (idx_list, pur, cov, passed) in enumerate(zip(indices, purity, coverage, passed_filter)):
        if len(idx_list) == 0:
            continue
        idx_clean = [int(x) for x in sorted(idx_list)]
        snp_list = snp_info.iloc[idx_clean]["snp"].tolist() if len(idx_clean) > 0 else []
        cs_records.append({
            "component": i,
            "num_variants": len(idx_clean),
            "purity": pur,
            "coverage": cov,
            "passed_filter": passed,
            "variant_indices": ",".join([str(x) for x in idx_clean]),
            "variants": ",".join(snp_list)
        })
    
    cs_df = pd.DataFrame(cs_records)
    cs_df.to_csv(f"{output_prefix}_credible_sets.tsv", sep="\t", index=False)
    print(f"Saved {len(cs_df)} credible sets to {output_prefix}_credible_sets.tsv")
    if len(cs_df) > 0:
        print(f"  Purity range: [{cs_df['purity'].min():.3f}, {cs_df['purity'].max():.3f}]")
        print(f"  Coverage range: [{cs_df['coverage'].min():.3f}, {cs_df['coverage'].max():.3f}]")
        print(f"  {(cs_df['passed_filter']).sum()} / {len(cs_df)} passed filtering")


def write_pip_effects_table(parsed, snp_info, output_prefix, populations=DEFAULT_POPULATIONS):
    """
    Save one merged table containing PIP and posterior effect sizes.
    
    Args:
        parsed: Parsed MultiSuSiE structures
        snp_info: DataFrame with SNP information
        output_prefix: Output prefix for output files
        populations: List of population codes (K populations)
    """
    print("\nWriting merged PIP + effect size table…")

    pip = parsed.get("pip")
    coef = parsed.get("coef")
    coef_sd = parsed.get("coef_sd")

    if pip is None:
        print("Warning: result.pip not available")
        return

    pip = np.asarray(pip, dtype=float)

    combined_effect_df = pd.DataFrame({
        "snp": snp_info["snp"],
        "position": snp_info["position"],
        "pip": pip,
    })

    if coef is None or coef_sd is None:
        print("Warning: result.coef or result.coef_sd not available; writing PIP-only columns")
        combined_effect_df = combined_effect_df.sort_values("pip", ascending=False)
        combined_effect_df.to_csv(f"{output_prefix}_pip_effects.tsv", sep="\t", index=False)
        print(f"Saved merged table to {output_prefix}_pip_effects.tsv")
        return

    for pop_idx, pop in enumerate(populations):
        beta = coef[pop_idx, :] if pop_idx < coef.shape[0] else np.full(len(snp_info), np.nan, dtype=float)
        se = coef_sd[pop_idx, :] if pop_idx < coef_sd.shape[0] else np.full(len(snp_info), np.nan, dtype=float)

        effect_df = pd.DataFrame({
            "snp": snp_info["snp"],
            "position": snp_info["position"],
            "beta": beta,
            "se": se,
        })
        effect_df["ci_lower"] = effect_df["beta"] - 1.96 * effect_df["se"]
        effect_df["ci_upper"] = effect_df["beta"] + 1.96 * effect_df["se"]

        effect_df = effect_df.rename(columns={
            "beta": f"beta_{pop}",
            "se": f"se_{pop}",
            "ci_lower": f"ci_lower_{pop}",
            "ci_upper": f"ci_upper_{pop}",
        })

        combined_effect_df = combined_effect_df.merge(effect_df, on=["snp", "position"], how="left")

    combined_effect_df = combined_effect_df.sort_values("pip", ascending=False)
    combined_effect_df.to_csv(f"{output_prefix}_pip_effects.tsv", sep="\t", index=False)
    print(f"Saved merged table to {output_prefix}_pip_effects.tsv")


def plot_pip_manhattan(z_list, snp_info, parsed, output_prefix, populations=DEFAULT_POPULATIONS):
    """
    Create Manhattan plots for each population with PIP as color.
    
    Args:
        z_list: List of z-score DataFrames
        snp_info: DataFrame with SNP information
        parsed: Parsed MultiSuSiE structures
        output_prefix: Output prefix for plot
        populations: List of population codes
    """
    print("\nGenerating PIP Manhattan plots…")

    pip = parsed.get("pip")
    if pip is None:
        print("Warning: could not extract PIP")
        return
    pip = np.asarray(pip, dtype=float)
    
    fig, axes = plt.subplots(len(populations), 1, figsize=(14, 3.5 * len(populations)))
    if len(populations) == 1:
        axes = [axes]
    
    pos = pd.to_numeric(snp_info["position"], errors="coerce").to_numpy(dtype=float)
    
    for idx, (pop, z_df, ax) in enumerate(zip(populations, z_list, axes)):
        if "beta" not in z_df.columns or "se" not in z_df.columns:
            ax.text(0.5, 0.5, f"Missing effect size columns for {pop}", 
                   ha="center", va="center", transform=ax.transAxes)
            continue
        
        beta = z_df["beta"].values
        se = z_df["se"].values
        valid = (se > 0) & np.isfinite(beta) & np.isfinite(se)
        
        if valid.sum() == 0:
            ax.text(0.5, 0.5, f"No valid effect sizes for {pop}", 
                   ha="center", va="center", transform=ax.transAxes)
            continue
        
        z = np.zeros(len(beta), dtype=float)
        z[valid] = beta[valid] / se[valid]
        pvals = 2.0 * stats.norm.sf(np.abs(z))
        logp = -np.log10(np.maximum(pvals, 1e-300))
        plot_mask = np.isfinite(pos) & np.isfinite(logp) & np.isfinite(pip)
        if plot_mask.sum() == 0:
            ax.text(0.5, 0.5, f"No finite points for {pop}", ha="center", va="center", transform=ax.transAxes)
            continue

        scatter = ax.scatter(pos[plot_mask], logp[plot_mask], c=pip[plot_mask], s=30, alpha=0.7, cmap="RdYlBu_r", vmin=0, vmax=1)
        ax.set_ylabel(f"-log10(p)\n{pop}", fontsize=10)
        ax.set_xlim(np.nanmin(pos) - 1000, np.nanmax(pos) + 1000)
        ax.grid(alpha=0.3, axis="y")
        plt.colorbar(scatter, ax=ax, label="PIP")
    
    axes[-1].set_xlabel("Position (bp)", fontsize=11)
    fig.suptitle("Manhattan plots colored by posterior inclusion probability", fontsize=12, y=0.995)
    fig.tight_layout()
    fig.savefig(f"{output_prefix}_pip_manhattan.png", dpi=160, bbox_inches="tight")
    print(f"Saved to {output_prefix}_pip_manhattan.png")
    plt.close(fig)


def plot_diagnostics(parsed, output_prefix, populations=DEFAULT_POPULATIONS):
    """
    Create one combined diagnostics figure.

    Panels:
      1) pairwise effect size comparisons
      2) credible set purity vs coverage
      3) credible set size by component
      4) component-wise LBF and KL (if available)
    """
    print("\nGenerating combined diagnostics plot…")

    coef = parsed.get("coef")
    pip = parsed.get("pip")
    sets_list = parsed.get("sets")
    diag = parsed.get("diag", {})

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    ax_scatter = axes[0, 0]
    ax_cs_quality = axes[0, 1]
    ax_cs_size = axes[1, 0]
    ax_conv = axes[1, 1]

    if coef is None or pip is None or coef.shape[0] < 2:
        ax_scatter.text(0.5, 0.5, "Effect comparison unavailable", ha="center", va="center", transform=ax_scatter.transAxes)
        ax_scatter.set_title("Cross-population effect comparison")
    else:
        pip = np.asarray(pip, dtype=float)
        beta_i = coef[0, :]
        beta_j = coef[1, :]
        mask = np.isfinite(beta_i) & np.isfinite(beta_j) & np.isfinite(pip)
        if mask.sum() == 0:
            ax_scatter.text(0.5, 0.5, "No finite points", ha="center", va="center", transform=ax_scatter.transAxes)
        else:
            scatter = ax_scatter.scatter(beta_i[mask], beta_j[mask], c=pip[mask], s=35, alpha=0.6, cmap="RdYlBu_r", vmin=0, vmax=1)
            lim = max(np.nanmax(np.abs(beta_i[mask])), np.nanmax(np.abs(beta_j[mask])))
            ax_scatter.plot([-lim, lim], [-lim, lim], "k--", lw=1, alpha=0.5)
            if mask.sum() > 1 and np.nanstd(beta_i[mask]) > 0 and np.nanstd(beta_j[mask]) > 0:
                r = np.corrcoef(beta_i[mask], beta_j[mask])[0, 1]
                ax_scatter.text(0.05, 0.95, f"r={r:.3f}", transform=ax_scatter.transAxes, ha="left", va="top", fontsize=10, bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))
            plt.colorbar(scatter, ax=ax_scatter, label="PIP")
        ax_scatter.set_xlabel(f"Effect size {populations[0]}", fontsize=10)
        ax_scatter.set_ylabel(f"Effect size {populations[1]}", fontsize=10)
        ax_scatter.set_title(f"Effect comparison: {populations[0]} vs {populations[1]}")
        ax_scatter.grid(alpha=0.3)

    if sets_list is None or len(sets_list) < 4:
        ax_cs_quality.text(0.5, 0.5, "Credible sets unavailable", ha="center", va="center", transform=ax_cs_quality.transAxes)
        ax_cs_size.text(0.5, 0.5, "Credible sets unavailable", ha="center", va="center", transform=ax_cs_size.transAxes)
    else:
        indices, purity, coverage, passed_filter = sets_list[:4]
        colors = ["green" if p else "red" for p in passed_filter]

        ax_cs_quality.scatter(purity, coverage, c=colors, s=90, alpha=0.7, edgecolors="black")
        ax_cs_quality.set_xlabel("Purity", fontsize=10)
        ax_cs_quality.set_ylabel("Coverage", fontsize=10)
        ax_cs_quality.set_title("Credible set quality")
        ax_cs_quality.set_xlim(0, 1)
        ax_cs_quality.set_ylim(0, 1)
        ax_cs_quality.grid(alpha=0.3)

        set_sizes = [len(idx_list) for idx_list in indices]
        ax_cs_size.bar(range(len(set_sizes)), set_sizes, color=colors, alpha=0.7)
        ax_cs_size.set_xlabel("Component", fontsize=10)
        ax_cs_size.set_ylabel("Number of variants", fontsize=10)
        ax_cs_size.set_title("Credible set sizes")
        ax_cs_size.grid(alpha=0.3, axis="y")

    lbf = diag.get("lbf")
    kl = diag.get("KL")
    if lbf is None and kl is None:
        ax_conv.text(0.5, 0.5, "Convergence traces unavailable", ha="center", va="center", transform=ax_conv.transAxes)
    else:
        if lbf is not None:
            lbf_arr = np.asarray(lbf).flatten()
            ax_conv.plot(range(len(lbf_arr)), lbf_arr, marker="o", label="LBF")
        if kl is not None:
            kl_arr = np.asarray(kl).flatten()
            ax_conv.plot(range(len(kl_arr)), kl_arr, marker="s", label="KL")
        ax_conv.legend()
    title_parts = ["Convergence by component"]
    if "converged" in diag:
        title_parts.append(f"converged={diag['converged']}")
    if "niter" in diag:
        title_parts.append(f"iters={diag['niter']}")
    ax_conv.set_title(" | ".join(title_parts))
    ax_conv.set_xlabel("Component", fontsize=10)
    ax_conv.set_ylabel("Value", fontsize=10)
    ax_conv.grid(alpha=0.3)

    fig.suptitle("MultiSuSiE diagnostics", fontsize=13, y=0.995)
    fig.tight_layout()
    fig.savefig(f"{output_prefix}_diagnostics.png", dpi=160, bbox_inches="tight")
    print(f"Saved to {output_prefix}_diagnostics.png")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="Run two-way MultiSuSiE on EUR and AFR fine-mapping data"
    )
    parser.add_argument(
        "--input-prefix", required=True,
        help="Input prefix from step2b (e.g., BACH2_finemap)"
    )
    parser.add_argument(
        "--output-prefix", required=True,
        help="Output prefix for results (e.g., BACH2_finemap)"
    )
    parser.add_argument(
        "--L", type=int, default=10,
        help="Number of components/credible sets (default: 10)"
    )
    parser.add_argument(
        "--max-iter", type=int, default=2000,
        help="Maximum iterations for MultiSuSiE (default: 2000)"
    )
    args = parser.parse_args()
    
    # Load finemap files
    populations = DEFAULT_POPULATIONS
    z_list, R_list, snp_info, corr_matrix, pop_sizes, maf_list = load_finemap_files(
        args.input_prefix,
        populations=populations,
    )
    
    # Run MultiSuSiE
    result = run_multisusie(z_list, R_list, rho=corr_matrix, 
                            pop_sizes=pop_sizes, maf_list=maf_list,
                            L=args.L, max_iter=args.max_iter)

    parsed = parse_multisusie_result(result, snp_info, populations)
    
    # Extract and save results
    extract_credible_sets(parsed, snp_info, args.output_prefix)
    write_pip_effects_table(parsed, snp_info, args.output_prefix, populations)
    
    # Generate visualizations
    plot_pip_manhattan(z_list, snp_info, parsed, args.output_prefix, populations)
    plot_diagnostics(parsed, args.output_prefix, populations)
    
    print("\nDone.")


if __name__ == "__main__":
    main()