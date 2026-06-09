#!/usr/bin/env python3
"""
Prepare summary statistics for LDAK's --sum-hers
Intersect with the tagging files, and output 
    For sample-size-weighted Fema (the all-cohort one doesn't have beta and se):
        Predictor A1 A2 n Z P Direction A1Freq
    For inverse-variance-weighted Fema :
        Predictor A1 A2 n BETA SE P A1Freq
(Predictor is the rsID; A1 is the test allele)

Tagging files are space delimited with header:
    Predictor A1 A2 Neighbours Tagging Weight MAF Categories Exp_Heritability Annotation_1 ... Annotation_65 Base
Make sure the output match the Predictor, A1, and A2. 
    Swap A1/A2 and flip BETA if needed to match the tagging file.
"""

from __future__ import annotations

import argparse
import csv
import glob
import gzip
import sys
import time
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


ALLELE_COMPLEMENT = {
    "A": "T",
    "T": "A",
    "C": "G",
    "G": "C",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Prepare summary statistics for LDAK --sum-hers by intersecting with "
            "tagging files and aligning alleles."
        )
    )
    parser.add_argument(
        "--sumstats",
        required=True,
        help="Path to summary statistics file (supports .gz).",
    )
    parser.add_argument(
        "--tagging",
        required=True,
        nargs="+",
        help=(
            "One or more tagging files (supports globs and .gz). "
            "Expected columns: Predictor A1 A2"
        ),
    )
    parser.add_argument("--out", required=True, help="Output path.")

    parser.add_argument("--snp-col", default="Predictor", help="SNP/rsID column in sumstats.")
    parser.add_argument("--a1-col", default="A1", help="Effect allele column in sumstats.")
    parser.add_argument("--a2-col", default="A2", help="Other allele column in sumstats.")
    parser.add_argument("--beta-col", default="BETA", help="Effect size column in sumstats (fema-iv).")
    parser.add_argument("--se-col", default="SE", help="Standard error column in sumstats (fema-iv).")
    parser.add_argument("--n-col", default="n", help="Sample size column in sumstats.")
    parser.add_argument("--a1freq-col", default="A1Freq", help="Effect (A1) allele frequency column in sumstats.")
    parser.add_argument("--p-col", default="P", help="P-value column in sumstats.")
    parser.add_argument("--z-col", default="Z", help="Z-score column in sumstats (required for fema-ss; used to infer Direction if --direction-col is empty).")
    parser.add_argument("--direction-col", default="", help="Direction column in sumstats (fema-ss, legacy). If empty, Direction is inferred from BETA or Z-score.")
    parser.add_argument(
        "--mode",
        required=True,
        choices=["fema-ss", "fema-iv", "cohort"],
        help=(
            "Output format: 'fema-ss' for sample-size-weighted "
            "(Predictor A1 A2 n A1Freq P Direction); "
            "'fema-iv' for inverse-variance-weighted "
            "(Predictor A1 A2 n A1Freq BETA SE); "
            "'cohort' for cohort-specific summary statistics "
            "(Predictor A1 A2 n A1Freq BETA SE)."
        ),
    )

    parser.add_argument(
        "--sumstats-delim",
        choices=["auto", "whitespace", "tab", "comma"],
        default="auto",
        help="Delimiter for sumstats input.",
    )
    parser.add_argument(
        "--allow-strand-flip",
        action="store_true",
        help="Allow strand-complement matching (A<->T, C<->G).",
    )
    parser.add_argument(
        "--allow-duplicates",
        action="store_true",
        help="If set, keep duplicate Predictor rows from sumstats in output.",
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=1000000,
        help=(
            "Print progress every N sumstats data rows (default: 1000000). "
            "Set to 0 to disable periodic progress logs."
        ),
    )
    return parser.parse_args()


def open_text(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return open(path, "rt", encoding="utf-8", newline="")


def expand_paths(patterns: Iterable[str]) -> List[str]:
    paths: List[str] = []
    for pattern in patterns:
        matches = sorted(glob.glob(pattern))
        if matches:
            paths.extend(matches)
        else:
            paths.append(pattern)
    unique_paths: List[str] = []
    seen = set()
    for p in paths:
        rp = str(Path(p))
        if rp not in seen:
            seen.add(rp)
            unique_paths.append(rp)
    return unique_paths


def detect_delim(header_line: str, mode: str) -> str:
    if mode == "whitespace":
        return "whitespace"
    if mode == "tab":
        return "\t"
    if mode == "comma":
        return ","
    if "\t" in header_line:
        return "\t"
    if "," in header_line:
        return ","
    return "whitespace"


def split_row(line: str, delim: str) -> List[str]:
    s = line.rstrip("\n\r")
    if delim == "whitespace":
        return s.strip().split()
    return next(csv.reader([s], delimiter=delim))


def uppercase_allele(a: str) -> str:
    return a.strip().upper()


def complement_allele(a: str) -> Optional[str]:
    a = uppercase_allele(a)
    if not a:
        return None
    out = []
    for ch in a:
        comp = ALLELE_COMPLEMENT.get(ch)
        if comp is None:
            return None
        out.append(comp)
    return "".join(out)


def clamp_pvalue(p_str: str) -> str:
    """Clamp p-values smaller than 1e-300 to 1e-299 to avoid underflow issues."""
    try:
        p = float(p_str.strip())
        if 0 < p < 1e-300:
            return "1e-299"
    except (ValueError, TypeError):
        pass
    return p_str


def flip_direction(direction: str) -> str:
    """Swap '+' and '-' signs in a Direction string; '?' and other characters are unchanged."""
    return direction.translate(str.maketrans("+-", "-+"))


def infer_direction_from_effect(value_str: str) -> str:
    """Infer Direction from effect size (beta or Z-score) sign. Returns '+', '-', or '?'."""
    try:
        value = float(value_str.strip())
        if value > 0:
            return "+"
        elif value < 0:
            return "-"
        else:
            return "?"
    except (ValueError, TypeError):
        return "?"


def resolve_col_indices(header: List[str], required_cols: List[str]) -> Dict[str, int]:
    col_to_idx = {c.lower(): i for i, c in enumerate(header)}
    out: Dict[str, int] = {}
    missing = []
    for col in required_cols:
        idx = col_to_idx.get(col.lower())
        if idx is None:
            missing.append(col)
        else:
            out[col] = idx
    if missing:
        raise ValueError(
            "Missing columns: " + ", ".join(missing) + "; available: " + ", ".join(header)
        )
    return out


def load_tagging_map(tagging_files: Iterable[str]) -> Tuple[Dict[str, Tuple[str, str]], Dict[str, int]]:
    tag_map: Dict[str, Tuple[str, str]] = {}
    stats = {
        "tag_rows": 0,
        "tag_loaded": 0,
        "tag_conflicts": 0,
        "tag_skipped_malformed": 0,
    }

    for path in tagging_files:
        with open_text(path) as f:
            header_line = None
            for raw in f:
                if raw.strip():
                    header_line = raw
                    break
            if header_line is None:
                continue

            header = split_row(header_line, "whitespace")
            try:
                idx = resolve_col_indices(header, ["Predictor", "A1", "A2"])
            except ValueError as e:
                raise ValueError(f"Tagging file {path}: {e}") from e

            for raw in f:
                if not raw.strip():
                    continue
                parts = split_row(raw, "whitespace")
                stats["tag_rows"] += 1

                max_idx = max(idx.values())
                if len(parts) <= max_idx:
                    stats["tag_skipped_malformed"] += 1
                    continue

                snp = parts[idx["Predictor"]].strip()
                a1 = uppercase_allele(parts[idx["A1"]])
                a2 = uppercase_allele(parts[idx["A2"]])

                if not snp or not a1 or not a2:
                    stats["tag_skipped_malformed"] += 1
                    continue

                prev = tag_map.get(snp)
                if prev is None:
                    tag_map[snp] = (a1, a2)
                    stats["tag_loaded"] += 1
                elif prev != (a1, a2):
                    stats["tag_conflicts"] += 1

    return tag_map, stats


def align_to_tagging(
    sum_a1: str,
    sum_a2: str,
    tag_a1: str,
    tag_a2: str,
    allow_strand_flip: bool,
) -> Tuple[bool, bool]:
    """
    Return (matched, should_flip_beta).
    """
    a1 = uppercase_allele(sum_a1)
    a2 = uppercase_allele(sum_a2)

    if a1 == tag_a1 and a2 == tag_a2:
        return True, False
    if a1 == tag_a2 and a2 == tag_a1:
        return True, True

    if not allow_strand_flip:
        return False, False

    c1 = complement_allele(a1)
    c2 = complement_allele(a2)
    if c1 is None or c2 is None:
        return False, False

    if c1 == tag_a1 and c2 == tag_a2:
        return True, False
    if c1 == tag_a2 and c2 == tag_a1:
        return True, True

    return False, False


def main() -> int:
    args = parse_args()
    t0 = time.time()

    def elapsed_s() -> str:
        return f"{time.time() - t0:.1f}s"

    print(f"[start] loading tagging files at +{elapsed_s()}")

    tagging_paths = expand_paths(args.tagging)
    missing_tags = [p for p in tagging_paths if not Path(p).exists()]
    if missing_tags:
        print("ERROR: Tagging file(s) not found:", file=sys.stderr)
        for p in missing_tags:
            print(f"  {p}", file=sys.stderr)
        return 1

    if not Path(args.sumstats).exists():
        print(f"ERROR: Sumstats file not found: {args.sumstats}", file=sys.stderr)
        return 1

    tag_map, tag_stats = load_tagging_map(tagging_paths)
    if not tag_map:
        print("ERROR: No valid Predictor/A1/A2 rows loaded from tagging files.")
        return 1
    print(
        f"[progress] loaded {tag_stats['tag_loaded']} unique tagging predictors at +{elapsed_s()}",
        file=sys.stderr,
    )

    stats = {
        "sum_rows": 0,
        "in_tagging": 0,
        "aligned": 0,
        "swapped": 0,
        "skipped_unaligned": 0,
        "skipped_bad_beta": 0,
        "skipped_duplicate": 0,
    }

    with open_text(args.sumstats) as fin:
        print(f"[start] reading sumstats at +{elapsed_s()}")
        header_line = None
        for raw in fin:
            if raw.strip():
                header_line = raw
                break
        if header_line is None:
            print("ERROR: Sumstats file is empty.", file=sys.stderr)
            return 1

        delim = detect_delim(header_line.rstrip("\n\r"), args.sumstats_delim)
        header = split_row(header_line, delim)
        common_cols = [args.snp_col, args.a1_col, args.a2_col, args.n_col, args.a1freq_col]
        if args.mode == "fema-ss":
            mode_cols = [args.z_col, args.p_col]
            if args.direction_col:
                mode_cols.append(args.direction_col)
            out_header = "Predictor A1 A2 n Z P Direction A1Freq\n"
        else:  # fema-iv or cohort
            mode_cols = [args.beta_col, args.se_col, args.p_col]
            out_header = "Predictor A1 A2 n BETA SE P A1Freq\n"
        required = common_cols + mode_cols
        try:
            idx = resolve_col_indices(header, required)
        except ValueError as e:
            print(f"ERROR: Sumstats file column check failed: {e}", file=sys.stderr)
            return 1

        seen_predictors = set()
        with open(args.out, "wt", encoding="utf-8", newline="") as fout:
            fout.write(out_header)

            for raw in fin:
                if not raw.strip():
                    continue

                parts = split_row(raw, delim)
                stats["sum_rows"] += 1

                if args.progress_every > 0 and stats["sum_rows"] % args.progress_every == 0:
                    print(
                        (
                            f"[progress] processed {stats['sum_rows']} rows; "
                            f"written {stats['aligned']} rows; "
                            f"elapsed +{elapsed_s()}"
                        ),
                        file=sys.stderr,
                    )

                max_idx = max(idx.values())
                if len(parts) <= max_idx:
                    continue

                snp = parts[idx[args.snp_col]].strip()
                if not snp:
                    continue

                target = tag_map.get(snp)
                if target is None:
                    continue

                stats["in_tagging"] += 1
                sum_a1 = parts[idx[args.a1_col]]
                sum_a2 = parts[idx[args.a2_col]]
                n_str = parts[idx[args.n_col]].strip()
                a1freq_str = parts[idx[args.a1freq_col]].strip()

                matched, should_flip = align_to_tagging(
                    sum_a1=sum_a1,
                    sum_a2=sum_a2,
                    tag_a1=target[0],
                    tag_a2=target[1],
                    allow_strand_flip=args.allow_strand_flip,
                )
                if not matched:
                    stats["skipped_unaligned"] += 1
                    continue

                if not args.allow_duplicates:
                    if snp in seen_predictors:
                        stats["skipped_duplicate"] += 1
                        continue
                    seen_predictors.add(snp)

                if should_flip:
                    try:
                        a1freq_str = format(1.0 - float(a1freq_str), ".15g")
                    except ValueError:
                        stats["skipped_bad_beta"] += 1
                        continue

                if args.mode == "fema-iv" or args.mode == "cohort":
                    beta_str = parts[idx[args.beta_col]].strip()
                    se_str = parts[idx[args.se_col]].strip()
                    p_str = clamp_pvalue(parts[idx[args.p_col]].strip())
                    if should_flip:
                        try:
                            beta_str = format(-float(beta_str), ".15g")
                        except ValueError:
                            stats["skipped_bad_beta"] += 1
                            continue
                    out_line = (
                        f"{snp} {target[0]} {target[1]} {n_str} {beta_str} {se_str} {p_str} {a1freq_str}\n"
                    )
                else:  # fema-ss
                    z_str = parts[idx[args.z_col]].strip()
                    p_str = clamp_pvalue(parts[idx[args.p_col]].strip())
                    if args.direction_col:
                        dir_str = parts[idx[args.direction_col]].strip()
                        if should_flip:
                            dir_str = flip_direction(dir_str)
                    else:
                        dir_str = infer_direction_from_effect(z_str)
                        if should_flip:
                            dir_str = flip_direction(dir_str)
                    if should_flip:
                        try:
                            z_str = format(-float(z_str), ".15g")
                        except ValueError:
                            stats["skipped_bad_beta"] += 1
                            continue
                    out_line = (
                        f"{snp} {target[0]} {target[1]} {n_str} {z_str} {p_str} {dir_str} {a1freq_str}\n"
                    )

                if should_flip:
                    stats["swapped"] += 1
                stats["aligned"] += 1
                fout.write(out_line)

    print(f"[done] processing completed at +{elapsed_s()}")

    print("Tagging summary:")
    print(f"  files: {len(tagging_paths)}")
    print(f"  rows read: {tag_stats['tag_rows']}")
    print(f"  unique predictors loaded: {tag_stats['tag_loaded']}")
    print(f"  conflicting predictor allele entries: {tag_stats['tag_conflicts']}")
    print(f"  malformed rows skipped: {tag_stats['tag_skipped_malformed']}")

    print("Sumstats summary:")
    print(f"  rows read: {stats['sum_rows']}")
    overlap_percent = (stats['in_tagging'] / tag_stats['tag_rows'] * 100) if tag_stats['tag_rows'] > 0 else 0
    print(f"  overlapping with tagging predictors: {stats['in_tagging']} ({overlap_percent:.2f}% of all tagging SNPs)")
    print(f"  aligned and written: {stats['aligned']}")
    print(f"  allele-swapped (fields flipped): {stats['swapped']}")
    print(f"  unaligned alleles skipped: {stats['skipped_unaligned']}")
    print(f"  bad beta rows skipped: {stats['skipped_bad_beta']}")
    print(f"  duplicate predictors skipped: {stats['skipped_duplicate']}")
    print(f"Output written to: {args.out}")
    print(f"Total runtime: {elapsed_s()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
