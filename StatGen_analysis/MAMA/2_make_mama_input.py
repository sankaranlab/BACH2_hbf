#!/usr/bin/env python3
"""Take in Thai summary stats and position-annotated METAL output files for EUR and AFR,
harmonize them, and output MAMA input files for the three populations.
Check the SNP IDs from Liam's `snp-ances-file.txt`, which has header: 
    SNP AFR EUR THAI (only SNP column is useful here, and is in format "chr:pos:non-effect_allele:effect_allele")
After harmonizing, align all SNPs by IDs in the `snp-ances-file.txt` file, 
    swapping alleles and flipping effect sizes as needed to match the ID, and
    dropping any SNPs that cannot be aligned (save them in a separate file).

FEMA header:  (MarkerName is rsID)
    CHR GENPOS MarkerName Allele1 Allele2 Freq1 FreqSE MinFreq MaxFreq
    Effect StdErr P-value Direction HetISq HetChiSq HetDf HetPVal TotalN

Thai stats (filtered) header: (ea_freq is the freq of a2)
    snp	chr	pos	a1	a2	id	ea_freq	info	beta	p	n	se	z	impute	maf

Output MAMA input header:
    SNP  CHR BP A1 A2 FREQ BETA SE P N
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import sys
import time
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, Iterator, Optional, Sequence, Tuple


COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}
OUTPUT_COLUMNS = ("SNP", "CHR", "BP", "A1", "A2", "FREQ", "BETA", "SE", "P", "N")
DROPPED_COLUMNS = ("SNP", "CHR", "BP", "REASON", "THAI_A1", "THAI_A2", "EUR_A1", "EUR_A2", "AFR_A1", "AFR_A2")


@dataclass(frozen=True)
class TargetSnp:
    snp_id: str
    chrom: str
    bp: int
    non_effect_allele: str
    effect_allele: str


@dataclass(frozen=True)
class StandardRecord:
    snp: str
    chrom: str
    bp: int
    effect_allele: str
    other_allele: str
    freq: float
    beta: float
    se: float
    p_value: float
    sample_size: float

    def to_mama_row(self, output_snp: str, ref_a1: str, ref_a2: str) -> Tuple[str, ...]:
        return (
            output_snp,
            self.chrom,
            str(self.bp),
            ref_a1,
            ref_a2,
            format_number(self.freq),
            format_number(self.beta),
            format_number(self.se),
            format_number(self.p_value),
            format_number(self.sample_size),
        )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Create harmonized MAMA input tables for Thai, EUR, and AFR summary statistics. "
            "EUR and AFR METAL inputs are aligned to the Thai effect allele."
        )
    )
    parser.add_argument("--thai", required=True, help="Thai summary statistics file (.tsv or .tsv.gz).")
    parser.add_argument("--eur", required=True, help="EUR position-annotated METAL output file.")
    parser.add_argument("--afr", required=True, help="AFR position-annotated METAL output file.")
    parser.add_argument(
        "--snp-ances-file",
        required=True,
        help="Liam's SNP ancestry file containing the SNP column formatted as chr:pos:non-effect:effect.",
    )
    parser.add_argument(
        "--out-prefix",
        default="mama_input",
        help="Prefix for output files. Outputs are written as <prefix>.THAI.tsv, <prefix>.EUR.tsv, and <prefix>.AFR.tsv.",
    )
    parser.add_argument(
        "--keep-palindromic",
        action="store_true",
        help="Keep A/T and C/G SNPs. By default these SNPs are dropped because strand orientation is ambiguous.",
    )
    parser.add_argument(
        "--dropped-out",
        default=None,
        help="Optional output path for SNPs dropped during alignment to snp-ances-file IDs (default: <prefix>.dropped.tsv).",
    )
    return parser


def open_text(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return open(path, "r", encoding="utf-8", newline="")


def normalize_chrom(chrom: str) -> str:
    value = chrom.strip()
    if value.lower().startswith("chr"):
        value = value[3:]
    return value


def normalize_allele(allele: str) -> str:
    return allele.strip().upper()


def is_palindromic(a1: str, a2: str) -> bool:
    return COMPLEMENT.get(a1) == a2


def valid_effect_alleles(a1: str, a2: str) -> bool:
    return a1 in COMPLEMENT and a2 in COMPLEMENT and a1 != a2


def parse_int(value: str) -> Optional[int]:
    text = value.strip()
    if not text or text.upper() == "NA":
        return None
    try:
        return int(float(text))
    except ValueError:
        return None


def parse_float(value: str) -> Optional[float]:
    text = value.strip()
    if not text or text.upper() == "NA":
        return None
    try:
        number = float(text)
    except ValueError:
        return None
    if math.isnan(number):
        return None
    return number


def format_number(value: float) -> str:
    if float(value).is_integer():
        return str(int(value))
    return f"{value:.10g}"


def format_duration(seconds: float) -> str:
    if seconds < 1:
        return f"{seconds * 1000:.0f}ms"
    if seconds < 60:
        return f"{seconds:.2f}s"
    minutes = int(seconds // 60)
    remain = seconds - (minutes * 60)
    return f"{minutes}m{remain:05.2f}s"


def log_progress(message: str, start_time: float, step_start: Optional[float] = None) -> float:
    now = time.perf_counter()
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    total_elapsed = format_duration(now - start_time)
    if step_start is None:
        step_part = ""
    else:
        step_part = f" | step {format_duration(now - step_start)}"
    print(f"[{timestamp}] {message} | total {total_elapsed}{step_part}", file=sys.stderr)
    return now


def resolve_snp_id(*values: str) -> Optional[str]:
    for value in values:
        text = value.strip()
        if text and text not in {".", "NA", "nan", "None"}:
            return text
    return None


def parse_target_snp_id(value: str) -> Optional[TargetSnp]:
    parts = value.strip().split(":")
    if len(parts) != 4:
        return None
    chrom = normalize_chrom(parts[0])
    bp = parse_int(parts[1])
    non_effect = normalize_allele(parts[2])
    effect = normalize_allele(parts[3])
    if not chrom or bp is None:
        return None
    if not valid_effect_alleles(effect, non_effect):
        return None
    return TargetSnp(snp_id=value.strip(), chrom=chrom, bp=bp, non_effect_allele=non_effect, effect_allele=effect)


def read_target_snps(path: str) -> Tuple[Sequence[TargetSnp], int]:
    targets = []
    invalid = 0
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "SNP" not in reader.fieldnames:
            raise ValueError("snp-ances-file must be tab-delimited and contain an SNP column")
        for row in reader:
            if row is None:
                continue
            parsed = parse_target_snp_id((row.get("SNP") or "").strip())
            if parsed is None:
                invalid += 1
                continue
            targets.append(parsed)
    return targets, invalid


def iter_delimited_rows(path: str, delimiter: Optional[str]) -> Iterator[dict]:
    with open_text(path) as handle:
        if delimiter is None:
            reader = csv.DictReader(handle, delimiter="\t")
        else:
            reader = csv.DictReader(handle, delimiter=delimiter, skipinitialspace=True)
        for row in reader:
            if row is None:
                continue
            yield {key.strip(): (value or "").strip() for key, value in row.items() if key is not None}


def read_thai_records(path: str) -> Dict[Tuple[str, int], StandardRecord]:
    records: Dict[Tuple[str, int], StandardRecord] = {}
    for row in iter_delimited_rows(path, "\t"):
        chrom = normalize_chrom(row.get("chr", ""))
        bp = parse_int(row.get("pos", ""))
        effect_allele = normalize_allele(row.get("a2", ""))
        other_allele = normalize_allele(row.get("a1", ""))
        freq = parse_float(row.get("ea_freq", ""))
        beta = parse_float(row.get("beta", ""))
        se = parse_float(row.get("se", ""))
        p_value = parse_float(row.get("p", ""))
        sample_size = parse_float(row.get("n", ""))
        snp = resolve_snp_id(row.get("id", ""), row.get("snp", "")) or ""

        if not chrom or bp is None:
            continue
        if not valid_effect_alleles(effect_allele, other_allele):
            continue
        if None in (freq, beta, se, p_value, sample_size):
            continue
        assert freq is not None and beta is not None and se is not None
        assert p_value is not None and sample_size is not None

        records[(chrom, bp)] = StandardRecord(
            snp=snp,
            chrom=chrom,
            bp=bp,
            effect_allele=effect_allele,
            other_allele=other_allele,
            freq=freq,
            beta=beta,
            se=se,
            p_value=p_value,
            sample_size=sample_size,
        )
    return records


def read_fema_records(path: str) -> Dict[Tuple[str, int], StandardRecord]:
    records: Dict[Tuple[str, int], StandardRecord] = {}
    for row in iter_delimited_rows(path, None):
        chrom = normalize_chrom(row.get("CHR", ""))
        bp = parse_int(row.get("GENPOS", ""))
        effect_allele = normalize_allele(row.get("Allele1", ""))
        other_allele = normalize_allele(row.get("Allele2", ""))
        freq = parse_float(row.get("Freq1", ""))
        beta = parse_float(row.get("Effect", ""))
        se = parse_float(row.get("StdErr", ""))
        p_value = parse_float(row.get("P-value", ""))
        sample_size = parse_float(row.get("TotalN", ""))
        snp = resolve_snp_id(row.get("MarkerName", "")) or ""

        if not chrom or bp is None:
            continue
        if not valid_effect_alleles(effect_allele, other_allele):
            continue
        if None in (freq, beta, se, p_value, sample_size):
            continue
        assert freq is not None and beta is not None and se is not None
        assert p_value is not None and sample_size is not None

        records[(chrom, bp)] = StandardRecord(
            snp=snp,
            chrom=chrom,
            bp=bp,
            effect_allele=effect_allele,
            other_allele=other_allele,
            freq=freq,
            beta=beta,
            se=se,
            p_value=p_value,
            sample_size=sample_size,
        )
    return records


def complement_allele(allele: str) -> str:
    return COMPLEMENT[allele]


def align_to_reference(record: StandardRecord, ref_a1: str, ref_a2: str) -> Optional[StandardRecord]:
    same = (record.effect_allele, record.other_allele)
    swapped = (record.other_allele, record.effect_allele)
    comp_same = (complement_allele(record.effect_allele), complement_allele(record.other_allele))
    comp_swapped = (complement_allele(record.other_allele), complement_allele(record.effect_allele))

    if same == (ref_a1, ref_a2) or comp_same == (ref_a1, ref_a2):
        return StandardRecord(
            snp=record.snp,
            chrom=record.chrom,
            bp=record.bp,
            effect_allele=ref_a1,
            other_allele=ref_a2,
            freq=record.freq,
            beta=record.beta,
            se=record.se,
            p_value=record.p_value,
            sample_size=record.sample_size,
        )

    if swapped == (ref_a1, ref_a2) or comp_swapped == (ref_a1, ref_a2):
        return StandardRecord(
            snp=record.snp,
            chrom=record.chrom,
            bp=record.bp,
            effect_allele=ref_a1,
            other_allele=ref_a2,
            freq=1.0 - record.freq,
            beta=-record.beta,
            se=record.se,
            p_value=record.p_value,
            sample_size=record.sample_size,
        )

    return None


def write_output(path: Path, rows: Iterable[Tuple[str, ...]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(OUTPUT_COLUMNS)
        writer.writerows(rows)


def write_dropped_output(path: Path, rows: Iterable[Tuple[str, ...]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(DROPPED_COLUMNS)
        writer.writerows(rows)


def build_harmonized_records(
    thai_records: Dict[Tuple[str, int], StandardRecord],
    eur_records: Dict[Tuple[str, int], StandardRecord],
    afr_records: Dict[Tuple[str, int], StandardRecord],
    keep_palindromic: bool,
) -> Tuple[Dict[Tuple[str, int], Tuple[StandardRecord, StandardRecord, StandardRecord]], dict]:
    harmonized: Dict[Tuple[str, int], Tuple[StandardRecord, StandardRecord, StandardRecord]] = {}
    stats = {
        "shared_sites": 0,
        "palindromic_dropped": 0,
        "allele_mismatch_dropped": 0,
        "harmonized": 0,
    }

    shared_keys = sorted(set(thai_records) & set(eur_records) & set(afr_records))
    stats["shared_sites"] = len(shared_keys)

    for key in shared_keys:
        thai = thai_records[key]
        eur = eur_records[key]
        afr = afr_records[key]
        ref_a1 = thai.effect_allele
        ref_a2 = thai.other_allele

        if not keep_palindromic and is_palindromic(ref_a1, ref_a2):
            stats["palindromic_dropped"] += 1
            continue

        harmonized_eur = align_to_reference(eur, ref_a1, ref_a2)
        harmonized_afr = align_to_reference(afr, ref_a1, ref_a2)
        if harmonized_eur is None or harmonized_afr is None:
            stats["allele_mismatch_dropped"] += 1
            continue

        thai = StandardRecord(
            snp=thai.snp,
            chrom=thai.chrom,
            bp=thai.bp,
            effect_allele=ref_a1,
            other_allele=ref_a2,
            freq=thai.freq,
            beta=thai.beta,
            se=thai.se,
            p_value=thai.p_value,
            sample_size=thai.sample_size,
        )
        harmonized_eur = StandardRecord(
            snp=harmonized_eur.snp,
            chrom=harmonized_eur.chrom,
            bp=harmonized_eur.bp,
            effect_allele=ref_a1,
            other_allele=ref_a2,
            freq=harmonized_eur.freq,
            beta=harmonized_eur.beta,
            se=harmonized_eur.se,
            p_value=harmonized_eur.p_value,
            sample_size=harmonized_eur.sample_size,
        )
        harmonized_afr = StandardRecord(
            snp=harmonized_afr.snp,
            chrom=harmonized_afr.chrom,
            bp=harmonized_afr.bp,
            effect_allele=ref_a1,
            other_allele=ref_a2,
            freq=harmonized_afr.freq,
            beta=harmonized_afr.beta,
            se=harmonized_afr.se,
            p_value=harmonized_afr.p_value,
            sample_size=harmonized_afr.sample_size,
        )

        harmonized[key] = (thai, harmonized_eur, harmonized_afr)
        stats["harmonized"] += 1

    return harmonized, stats


def build_target_aligned_rows(
    harmonized: Dict[Tuple[str, int], Tuple[StandardRecord, StandardRecord, StandardRecord]],
    targets: Sequence[TargetSnp],
) -> Tuple[Sequence[Tuple[str, ...]], Sequence[Tuple[str, ...]], Sequence[Tuple[str, ...]], Sequence[Tuple[str, ...]], dict]:
    thai_rows = []
    eur_rows = []
    afr_rows = []
    dropped_rows = []
    stats = {
        "target_total": len(targets),
        "target_missing_dropped": 0,
        "target_allele_mismatch_dropped": 0,
        "written": 0,
    }

    for target in targets:
        key = (target.chrom, target.bp)
        trio = harmonized.get(key)
        if trio is None:
            stats["target_missing_dropped"] += 1
            dropped_rows.append((target.snp_id, target.chrom, str(target.bp), "NOT_IN_HARMONIZED_SHARED_SET", "", "", "", "", "", ""))
            continue

        thai_rec, eur_rec, afr_rec = trio
        ref_a1 = target.effect_allele
        ref_a2 = target.non_effect_allele
        thai_aligned = align_to_reference(thai_rec, ref_a1, ref_a2)
        eur_aligned = align_to_reference(eur_rec, ref_a1, ref_a2)
        afr_aligned = align_to_reference(afr_rec, ref_a1, ref_a2)

        if thai_aligned is None or eur_aligned is None or afr_aligned is None:
            stats["target_allele_mismatch_dropped"] += 1
            dropped_rows.append(
                (
                    target.snp_id,
                    target.chrom,
                    str(target.bp),
                    "ALLELES_CANNOT_ALIGN_TO_TARGET_ID",
                    thai_rec.effect_allele,
                    thai_rec.other_allele,
                    eur_rec.effect_allele,
                    eur_rec.other_allele,
                    afr_rec.effect_allele,
                    afr_rec.other_allele,
                )
            )
            continue

        thai_rows.append(thai_aligned.to_mama_row(target.snp_id, ref_a1, ref_a2))
        eur_rows.append(eur_aligned.to_mama_row(target.snp_id, ref_a1, ref_a2))
        afr_rows.append(afr_aligned.to_mama_row(target.snp_id, ref_a1, ref_a2))
        stats["written"] += 1

    return thai_rows, eur_rows, afr_rows, dropped_rows, stats


def print_summary(parse_counts: dict, harmonized_stats: dict, target_stats: dict) -> None:
    print(f"Thai records retained after parsing: {parse_counts['thai']}", file=sys.stderr)
    print(f"EUR records retained after parsing: {parse_counts['eur']}", file=sys.stderr)
    print(f"AFR records retained after parsing: {parse_counts['afr']}", file=sys.stderr)
    print(f"Invalid target SNP IDs in snp-ances-file: {parse_counts['target_invalid']}", file=sys.stderr)
    print(f"Sites shared across all three populations: {harmonized_stats['shared_sites']}", file=sys.stderr)
    print(f"Dropped palindromic sites: {harmonized_stats['palindromic_dropped']}", file=sys.stderr)
    print(f"Dropped allele-mismatched sites in first harmonization: {harmonized_stats['allele_mismatch_dropped']}", file=sys.stderr)
    print(f"Sites retained after first harmonization: {harmonized_stats['harmonized']}", file=sys.stderr)
    print(f"Target SNP IDs considered: {target_stats['target_total']}", file=sys.stderr)
    print(f"Target IDs missing from harmonized set: {target_stats['target_missing_dropped']}", file=sys.stderr)
    print(f"Target IDs failing allele alignment: {target_stats['target_allele_mismatch_dropped']}", file=sys.stderr)
    print(f"Variants written to each output: {target_stats['written']}", file=sys.stderr)


def main(argv: Optional[Sequence[str]] = None) -> int:
    start_time = time.perf_counter()
    step_start = log_progress("Pipeline started", start_time)

    step_start = log_progress("Parsing command-line arguments", start_time, step_start)
    args = build_parser().parse_args(argv)

    step_start = log_progress("Reading Thai summary statistics", start_time, step_start)
    thai_records = read_thai_records(args.thai)
    step_start = log_progress(f"Thai records loaded: {len(thai_records)}", start_time, step_start)

    step_start = log_progress("Reading EUR METAL summary statistics", start_time, step_start)
    eur_records = read_fema_records(args.eur)
    step_start = log_progress(f"EUR records loaded: {len(eur_records)}", start_time, step_start)

    step_start = log_progress("Reading AFR METAL summary statistics", start_time, step_start)
    afr_records = read_fema_records(args.afr)
    step_start = log_progress(f"AFR records loaded: {len(afr_records)}", start_time, step_start)

    step_start = log_progress("Reading target SNP IDs from snp-ances-file", start_time, step_start)
    target_snps, invalid_target_count = read_target_snps(args.snp_ances_file)
    step_start = log_progress(
        f"Target SNP IDs loaded: {len(target_snps)} (invalid skipped: {invalid_target_count})",
        start_time,
        step_start,
    )

    step_start = log_progress("Harmonizing shared variants across THAI/EUR/AFR", start_time, step_start)
    harmonized, harmonized_stats = build_harmonized_records(
        thai_records,
        eur_records,
        afr_records,
        keep_palindromic=args.keep_palindromic,
    )
    step_start = log_progress(f"Variants retained after first harmonization: {harmonized_stats['harmonized']}", start_time, step_start)

    step_start = log_progress("Aligning harmonized variants to target SNP IDs", start_time, step_start)
    thai_rows, eur_rows, afr_rows, dropped_rows, target_stats = build_target_aligned_rows(
        harmonized,
        target_snps,
    )
    step_start = log_progress(f"Variants aligned to target IDs: {target_stats['written']}", start_time, step_start)

    out_prefix = Path(args.out_prefix)

    step_start = log_progress("Writing THAI output", start_time, step_start)
    write_output(out_prefix.with_suffix(".THAI.tsv"), thai_rows)

    step_start = log_progress("Writing EUR output", start_time, step_start)
    write_output(out_prefix.with_suffix(".EUR.tsv"), eur_rows)

    step_start = log_progress("Writing AFR output", start_time, step_start)
    write_output(out_prefix.with_suffix(".AFR.tsv"), afr_rows)

    dropped_out = Path(args.dropped_out) if args.dropped_out else out_prefix.with_suffix(".dropped.tsv")
    step_start = log_progress("Writing dropped SNP report", start_time, step_start)
    write_dropped_output(dropped_out, dropped_rows)
    step_start = log_progress(f"Dropped SNP report written: {dropped_out}", start_time, step_start)

    parse_counts = {
        "thai": len(thai_records),
        "eur": len(eur_records),
        "afr": len(afr_records),
        "target_invalid": invalid_target_count,
    }
    step_start = log_progress("Printing summary statistics", start_time, step_start)
    print_summary(parse_counts, harmonized_stats, target_stats)
    log_progress("Pipeline finished", start_time, step_start)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())