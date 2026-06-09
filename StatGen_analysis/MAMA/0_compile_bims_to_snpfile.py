#!/usr/bin/env python3
"""Build SNP ancestry presence table from PLINK .bim files.

The output format matches MAMA snp-ances-file.txt:

SNP\tAFR\tEUR\tTHAI
1:14464:A:T\t1\t0\t0
...
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import Dict, List, Set, Tuple


def parse_args() -> argparse.Namespace:
	parser = argparse.ArgumentParser(
		description=(
			"Read .bim files for each ancestry and chromosome, then create "
			"snp-ances-file.txt with ancestry presence indicators."
		)
	)
	parser.add_argument(
		"--panel-dir",
		default="MAMA/PLINK_panel",
		help=(
			"Directory containing ancestry subfolders with .bim files "
			"(default: MAMA/PLINK_panel)."
		),
	)
	parser.add_argument(
		"--ancestries",
		nargs="+",
		default=["afr", "eur", "thai"],
		help="Ancestry labels to process (default: afr eur thai).",
	)
	parser.add_argument(
		"--chrs",
		nargs="+",
		type=int,
		default=list(range(1, 23)),
		help="Chromosomes to process (default: 1..22).",
	)
	parser.add_argument(
		"--out",
		default="MAMA/input/snp-ances-file.txt",
		help="Output file path (default: MAMA/input/snp-ances-file.txt).",
	)
	return parser.parse_args()


def normalize_chr(raw_chr: str) -> str:
	s = str(raw_chr).strip()
	if s.lower().startswith("chr"):
		s = s[3:]
	return s


def parse_bim_line(line: str, source_file: str, line_num: int) -> Tuple[str, int, str, str]:
	parts = line.rstrip("\n").split()
	if len(parts) < 6:
		raise ValueError(
			f"Malformed .bim line in {source_file}:{line_num}. Expected 6 columns, got {len(parts)}"
		)
	chr_raw, _snp, _cm, bp_raw, a1, a2 = parts[:6]
	chr_clean = normalize_chr(chr_raw)
	try:
		bp = int(bp_raw)
	except ValueError as exc:
		raise ValueError(
			f"Non-integer BP in {source_file}:{line_num}: {bp_raw}"
		) from exc

	# Keep allele order as in .bim to preserve exact SNP identity used downstream.
	return (chr_clean, bp, a1.upper(), a2.upper())


def read_chr_snps(bim_path: str) -> Set[Tuple[str, int, str, str]]:
	"""Read one BIM file into a set of unique (chr, bp, a1, a2) tuples."""
	snps: Set[Tuple[str, int, str, str]] = set()
	with open(bim_path, "r", encoding="utf-8") as fh:
		for i, line in enumerate(fh, start=1):
			if not line.strip():
				continue
			snps.add(parse_bim_line(line, bim_path, i))
	return snps


def write_chr_rows(
	out_handle,
	chr_snps_by_ances: Dict[str, Set[Tuple[str, int, str, str]]],
	ancestry_order: List[str],
) -> int:
	"""Write one chromosome worth of rows and return written row count."""
	all_snps: Set[Tuple[str, int, str, str]] = set()
	for snps in chr_snps_by_ances.values():
		all_snps.update(snps)

	rows_written = 0
	for ch, bp, a1, a2 in sorted(all_snps, key=lambda x: (x[1], x[2], x[3])):
		snp_id = f"{ch}:{bp}:{a1}:{a2}"
		flags = ["1" if (ch, bp, a1, a2) in chr_snps_by_ances[a] else "0" for a in ancestry_order]
		out_handle.write("\t".join([snp_id] + flags) + "\n")
		rows_written += 1

	return rows_written


def main() -> int:
	args = parse_args()
	ancestry_order = [a.lower() for a in args.ancestries]
	os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
	rows_total = 0
	files_found = 0
	missing_files: List[str] = []

	with open(args.out, "w", encoding="utf-8") as out:
		out.write("\t".join(["SNP"] + [a.upper() for a in ancestry_order]) + "\n")

		for chrom in args.chrs:
			chr_snps_by_ances: Dict[str, Set[Tuple[str, int, str, str]]] = {}
			for ances in ancestry_order:
				bim_path = os.path.join(args.panel_dir, ances, f"{ances}_chr{chrom}.bim")
				if not os.path.exists(bim_path):
					missing_files.append(bim_path)
					chr_snps_by_ances[ances] = set()
					continue

				files_found += 1
				chr_snps_by_ances[ances] = read_chr_snps(bim_path)

			rows_total += write_chr_rows(out, chr_snps_by_ances, ancestry_order)

	if files_found == 0:
		print("Error: No .bim files found. Check --panel-dir and filenames.\n", file=sys.stderr)
		return 1

	if missing_files:
		print(f"Warning: Missing {len(missing_files)} expected .bim files.\n", file=sys.stderr)

	print(f"Wrote {rows_total} SNP rows from {files_found} .bim files to {args.out}\n")
	return 0


if __name__ == "__main__":
	raise SystemExit(main())
