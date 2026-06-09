#!/usr/bin/env python3
"""
Merge liftover coordinates into GWAS summary stats with integrated dbSNP validation
and filtering.

Supported GWAS headers include:
    snp  chr  pos  a1  a2  id  maf  info  beta  p  n  se  z  impute
    and MAMA-style output:
    SNP	CHR	BP	A1	A2	FREQ	BETA	SE	Z	P	N_EFF	N_ORIG

Liftover BED format (UCSC):
    target_chr  target_start  target_end  orig_chr:orig_pos-orig_pos  score

Always runs the full merge + filter pipeline:
  1. Merges liftover BED coordinates into GWAS rows (prepends <to>_chr, <to>_pos, <to>_id).
  2. Applies fast liftover sanity checks (missing match, non-chromosomal contig).
  3. Batch-queries dbSNP build 157 GRCh37 by rsID or lifted coordinate (bcftools),
      backfills missing SNP rsIDs when available, and filters coordinate mismatches.
  4. Writes passing variants (sorted by target coordinates) to --output-file.
  5. Writes filtered variants with reasons to --filtered-file.
  Optionally also writes the full merged (pre-filter) table via --merged-output-file.
"""

import argparse
import gzip
import os
import subprocess
import sys
import tempfile
import time
from multiprocessing import Pool

DBSNP_VCF_HG19 = "/lab-share/Hem-Sankaran-e2/Public/public_datasets/dbSNP/build_157/VCF/GCF_000001405.25.gz"

GRCH37_REFSEQ_PREFIXES = {
    "1": "NC_000001.",
    "2": "NC_000002.",
    "3": "NC_000003.",
    "4": "NC_000004.",
    "5": "NC_000005.",
    "6": "NC_000006.",
    "7": "NC_000007.",
    "8": "NC_000008.",
    "9": "NC_000009.",
    "10": "NC_000010.",
    "11": "NC_000011.",
    "12": "NC_000012.",
    "13": "NC_000013.",
    "14": "NC_000014.",
    "15": "NC_000015.",
    "16": "NC_000016.",
    "17": "NC_000017.",
    "18": "NC_000018.",
    "19": "NC_000019.",
    "20": "NC_000020.",
    "21": "NC_000021.",
    "22": "NC_000022.",
    "X": "NC_000023.",
    "Y": "NC_000024.",
    "M": "NC_012920.",
    "MT": "NC_012920.",
}

_WORKER_COORD_MAP = None
_WORKER_GWAS_SCHEMA = None


def normalize_chrom(chrom):
    return str(chrom).replace("chr", "")


def canonical_chrom_label(chrom):
    chrom = normalize_chrom(chrom).upper()
    if chrom == "MT":
        return "M"
    return chrom


def split_gwas_fields(line):
    return line.rstrip("\n").split()


def debug_log(enabled, message):
    if enabled:
        print(f"{time.ctime()}. DEBUG: {message}")


def normalize_header_name(name):
    return str(name).strip().upper()


def extract_first_rsid(value):
    for token in str(value).replace(";", ",").split(","):
        token = token.strip()
        if token.lower().startswith("rs"):
            return token
    return None


def detect_gwas_schema(raw_header_line):
    fields = split_gwas_fields(raw_header_line)
    lookup = {normalize_header_name(field): idx for idx, field in enumerate(fields)}

    def find_index(*names, required=True):
        for name in names:
            idx = lookup.get(name)
            if idx is not None:
                return idx
        if required:
            wanted = ", ".join(names)
            raise ValueError(f"Input GWAS header is missing required column(s): {wanted}")
        return None

    return {
        "fields": fields,
        "chrom_idx": find_index("CHR", "CHROM"),
        "pos_idx": find_index("POS", "BP"),
        "a1_idx": find_index("A1"),
        "a2_idx": find_index("A2"),
        "snp_idx": find_index("SNP", required=False),
        "id_idx": find_index("ID", required=False),
    }


def is_chromosomal(chrom):
    chrom = str(chrom).upper()
    if chrom in {"X", "Y", "M", "MT"}:
        return True
    try:
        num = int(chrom)
        return 1 <= num <= 22
    except ValueError:
        return False


def chrom_to_int(chrom):
    chrom = canonical_chrom_label(chrom)
    if chrom == "X":
        return 23
    if chrom == "Y":
        return 24
    if chrom == "M":
        return 25
    try:
        return int(chrom)
    except ValueError:
        return 999


def parse_liftover_line(line):
    if line.startswith("#") or line.startswith("track") or line.startswith("browser"):
        return None
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 4:
        return None

    name_field = fields[3]
    try:
        orig_chrom = normalize_chrom(name_field.split(":")[0])
        orig_pos = int(name_field.split(":")[1].split("-")[0])
    except (IndexError, ValueError):
        return None

    new_chrom = normalize_chrom(fields[0])
    try:
        new_pos = int(fields[2])
    except ValueError:
        return None

    return orig_chrom, orig_pos, new_chrom, new_pos


def load_liftover_map(liftover_bed):
    bed_open = gzip.open if liftover_bed.endswith(".gz") else open
    coord_map = {}
    n_dup = 0

    with bed_open(liftover_bed, "rt") as fh:
        for line in fh:
            parsed = parse_liftover_line(line)
            if not parsed:
                continue
            orig_chrom, orig_pos, new_chrom, new_pos = parsed
            key = (orig_chrom, orig_pos)
            if key in coord_map:
                n_dup += 1
                continue
            coord_map[key] = (new_chrom, new_pos)

    if n_dup:
        print(f"{time.ctime()}. Warning: {n_dup} duplicate positions in BED (kept first mapping).")
    print(f"{time.ctime()}. Loaded {len(coord_map):,} liftover mappings.")
    return coord_map


def init_worker(coord_map, gwas_schema):
    global _WORKER_COORD_MAP
    global _WORKER_GWAS_SCHEMA
    _WORKER_COORD_MAP = coord_map
    _WORKER_GWAS_SCHEMA = gwas_schema


def extract_candidate_rsid(fields, gwas_schema):
    snp_idx = gwas_schema.get("snp_idx")
    id_idx = gwas_schema.get("id_idx")

    if snp_idx is not None and len(fields) > snp_idx:
        rsid = extract_first_rsid(fields[snp_idx])
        if rsid:
            return rsid
    if id_idx is not None and len(fields) > id_idx:
        rsid = extract_first_rsid(fields[id_idx])
        if rsid:
            return rsid
    return "NA"


def should_add_rsid_column(gwas_schema):
    return gwas_schema.get("snp_idx") is not None and gwas_schema.get("id_idx") is None


def build_output_header_fields(raw_header_line, from_assembly, to_assembly, gwas_schema, add_rsid_column=False):
    fields = split_gwas_fields(raw_header_line)
    if len(fields) <= max(gwas_schema["chrom_idx"], gwas_schema["pos_idx"]):
        return fields

    fields[gwas_schema["chrom_idx"]] = f"{from_assembly}_chr"
    fields[gwas_schema["pos_idx"]] = f"{from_assembly}_pos"

    new_cols = [f"{to_assembly}_chr", f"{to_assembly}_pos", f"{to_assembly}_id"]
    if add_rsid_column:
        new_cols.append("rsID")
    return new_cols + fields


def build_output_row_fields(merged_fields, resolved_rsid, add_rsid_column=False):
    if not add_rsid_column:
        return list(merged_fields)
    return list(merged_fields[:3]) + [resolved_rsid] + list(merged_fields[3:])


def transform_data_line(line):
    """Transform one GWAS row using worker-global liftover map.

    Returns tuple:
      (count_total, count_mapped, merged_line, to_chr, to_pos, to_id, rsid, fast_filter_reason)
    """
    global _WORKER_COORD_MAP
    global _WORKER_GWAS_SCHEMA

    fields = split_gwas_fields(line)
    gwas_schema = _WORKER_GWAS_SCHEMA
    if gwas_schema is None:
        raise RuntimeError("Worker GWAS schema was not initialized.")

    min_fields = max(gwas_schema["chrom_idx"], gwas_schema["pos_idx"]) + 1
    if len(fields) < min_fields:
        return 0, 0, None, None, None, None, "NA", "Malformed_input_row"

    chrom = normalize_chrom(fields[gwas_schema["chrom_idx"]])
    try:
        pos = int(fields[gwas_schema["pos_idx"]])
    except ValueError:
        return 0, 0, None, None, None, None, "NA", "Invalid_input_position"

    rsid = extract_candidate_rsid(fields, gwas_schema)

    mapping = _WORKER_COORD_MAP.get((chrom, pos))
    if mapping:
        new_chrom, new_pos = mapping
        a1_idx = gwas_schema.get("a1_idx")
        a2_idx = gwas_schema.get("a2_idx")
        a1 = fields[a1_idx] if a1_idx is not None and len(fields) > a1_idx else "NA"
        a2 = fields[a2_idx] if a2_idx is not None and len(fields) > a2_idx else "NA"
        new_id = f"chr{new_chrom}:{new_pos}:{a1}:{a2}"
        merged_fields = [str(new_chrom), str(new_pos), new_id] + fields

        if not is_chromosomal(new_chrom):
            reason = f"Liftover_to_non_chromosomal_contig_{new_chrom}"
        elif "," in new_id or ";" in new_id:
            reason = "Multiple_liftover_matches"
        else:
            reason = None
        return 1, 1, merged_fields, str(new_chrom), int(new_pos), new_id, rsid, reason

    merged_fields = ["NA", "NA", "NA"] + fields
    return 1, 0, merged_fields, "NA", None, "NA", rsid, "No_liftover_match"


def select_best_mapping(existing, new):
    if existing is None:
        return new
    old_chr, old_pos = existing
    new_chr, new_pos = new

    old_rank = chrom_to_int(old_chr)
    new_rank = chrom_to_int(new_chr)
    if new_rank < old_rank:
        return new
    if new_rank == old_rank and int(new_pos) < int(old_pos):
        return new
    return existing


def select_preferred_rsid(existing, new):
    if existing is None:
        return new

    def sort_key(rsid):
        suffix = rsid[2:]
        if suffix.isdigit():
            return (0, int(suffix))
        return (1, rsid)

    return new if sort_key(new) < sort_key(existing) else existing


def fetch_dbsnp_contigs(dbsnp_vcf):
    cmd = f"bcftools view -h '{dbsnp_vcf}'"
    proc = subprocess.run(
        cmd,
        shell=True,
        check=True,
        text=True,
        capture_output=True,
    )

    contigs = []
    for line in proc.stdout.splitlines():
        if not line.startswith("##contig=<ID="):
            continue
        contig = line.split("##contig=<ID=", 1)[1].split(",", 1)[0].split(">", 1)[0].strip()
        if contig:
            contigs.append(contig)
    return contigs


def choose_dbsnp_contig(canonical_chrom, available_contigs):
    canonical_chrom = canonical_chrom_label(canonical_chrom)
    direct_candidates = [canonical_chrom, f"chr{canonical_chrom}"]
    if canonical_chrom == "M":
        direct_candidates.extend(["MT", "chrM", "chrMT"])

    for candidate in direct_candidates:
        if candidate in available_contigs:
            return candidate

    refseq_prefix = GRCH37_REFSEQ_PREFIXES.get(canonical_chrom)
    if refseq_prefix:
        for contig in available_contigs:
            if contig.startswith(refseq_prefix):
                return contig

    for contig in available_contigs:
        if canonical_chrom_label(contig) == canonical_chrom:
            return contig

    return canonical_chrom


def build_dbsnp_chrom_maps(dbsnp_vcf, debug=False):
    available_contigs = fetch_dbsnp_contigs(dbsnp_vcf)
    query_map = {}
    reverse_map = {}

    canonical_chroms = [str(chrom) for chrom in range(1, 23)] + ["X", "Y", "M"]
    for chrom in canonical_chroms:
        contig = choose_dbsnp_contig(chrom, available_contigs)
        query_map[chrom] = contig
        reverse_map[contig] = chrom

    debug_log(debug, f"dbSNP contig sample: {available_contigs[:5]}")
    debug_log(debug, f"dbSNP chromosome map sample: {[(chrom, query_map[chrom]) for chrom in ['1', '2', 'X', 'M']]}")
    return query_map, reverse_map


def convert_chrom_for_dbsnp(chrom, chrom_query_map):
    canonical = canonical_chrom_label(chrom)
    return chrom_query_map.get(canonical, chrom)


def canonicalize_dbsnp_chrom(chrom, chrom_reverse_map):
    if chrom in chrom_reverse_map:
        return chrom_reverse_map[chrom]
    normalized = canonical_chrom_label(chrom)
    return normalized if is_chromosomal(normalized) else None


def build_dbsnp_lookup(rsids, dbsnp_vcf, threads=1, debug=False, chrom_reverse_map=None):
    if not rsids:
        return {}

    sample_rsids = sorted(rsids)[:5]

    with tempfile.TemporaryDirectory(prefix="step1c_dbsnp_") as tmpdir:
        rsid_file = os.path.join(tmpdir, "rsids.txt")
        with open(rsid_file, "w", encoding="utf-8") as fh:
            fh.write("\n".join(sorted(rsids)))
            fh.write("\n")

        cmd = (
            f"bcftools view --threads {max(1, int(threads))} -i 'ID=@{rsid_file}' '{dbsnp_vcf}' -Ou "
            f"| bcftools query -f '%ID\\t%CHROM\\t%POS\\n'"
        )
        proc = subprocess.run(
            cmd,
            shell=True,
            check=True,
            text=True,
            capture_output=True,
        )

    debug_log(debug, f"rsID lookup queried {len(rsids):,} unique rsIDs.")
    debug_log(debug, f"rsID lookup sample rsIDs: {sample_rsids}")
    debug_log(debug, f"rsID lookup raw result lines: {len(proc.stdout.splitlines()):,}")
    if proc.stdout:
        debug_log(debug, f"rsID lookup raw stdout sample: {proc.stdout.splitlines()[:5]}")
    if proc.stderr:
        debug_log(debug, f"rsID lookup stderr sample: {proc.stderr.splitlines()[:5]}")

    lookup = {}
    for line in proc.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) != 3:
            continue
        id_field, chrom, pos = fields
        chrom = canonicalize_dbsnp_chrom(chrom, chrom_reverse_map or {})
        if chrom is None:
            continue
        try:
            pos_int = int(pos)
        except ValueError:
            continue
        for rsid in str(id_field).replace(";", ",").split(","):
            rsid = rsid.strip()
            if not rsid.lower().startswith("rs"):
                continue
            lookup[rsid] = select_best_mapping(lookup.get(rsid), (chrom, pos_int))

    if debug and not lookup and sample_rsids:
        debug_log(debug, "rsID lookup returned zero records for the sampled rsIDs.")

    return lookup


def debug_probe_dbsnp_by_region(coords, dbsnp_vcf, chrom_query_map=None, debug=False, label="dbSNP region probe"):
    if not debug or not coords:
        return

    sample_regions = ",".join(
        f"{convert_chrom_for_dbsnp(chrom, chrom_query_map or {})}:{pos}-{pos}"
        for chrom, pos in coords[:5]
    )
    cmd = f"bcftools query -r '{sample_regions}' -f '%CHROM\\t%POS\\t%ID\\n' '{dbsnp_vcf}'"
    proc = subprocess.run(
        cmd,
        shell=True,
        check=False,
        text=True,
        capture_output=True,
    )
    debug_log(debug, f"{label} regions={sample_regions}")
    debug_log(debug, f"{label} exit_code={proc.returncode}")
    if proc.stdout:
        debug_log(debug, f"{label} stdout sample: {proc.stdout.splitlines()[:5]}")
    if proc.stderr:
        debug_log(debug, f"{label} stderr sample: {proc.stderr.splitlines()[:5]}")


def debug_probe_dbsnp_by_window(coords, dbsnp_vcf, chrom_query_map=None, debug=False, label="dbSNP window probe"):
    if not debug or not coords:
        return

    sample_regions = ",".join(
        f"{convert_chrom_for_dbsnp(chrom, chrom_query_map or {})}:{max(1, pos - 1)}-{pos + 1}"
        for chrom, pos in coords[:5]
    )
    cmd = f"bcftools query -r '{sample_regions}' -f '%CHROM\\t%POS\\t%ID\\n' '{dbsnp_vcf}'"
    proc = subprocess.run(
        cmd,
        shell=True,
        check=False,
        text=True,
        capture_output=True,
    )
    debug_log(debug, f"{label} regions={sample_regions}")
    debug_log(debug, f"{label} exit_code={proc.returncode}")
    if proc.stdout:
        debug_log(debug, f"{label} stdout sample: {proc.stdout.splitlines()[:10]}")
    if proc.stderr:
        debug_log(debug, f"{label} stderr sample: {proc.stderr.splitlines()[:5]}")


def build_dbsnp_coord_lookup(coords, dbsnp_vcf, threads=1, debug=False, chrom_query_map=None, chrom_reverse_map=None):
    if not coords:
        return {}

    sorted_coords = sorted(coords, key=lambda item: (chrom_to_int(item[0]), item[1]))

    with tempfile.TemporaryDirectory(prefix="step1c_dbsnp_coords_") as tmpdir:
        regions_file = os.path.join(tmpdir, "regions.tsv")
        with open(regions_file, "w", encoding="utf-8") as fh:
            for chrom, pos in sorted_coords:
                dbsnp_chrom = convert_chrom_for_dbsnp(chrom, chrom_query_map or {})
                fh.write(f"{dbsnp_chrom}\t{pos}\t{pos}\n")

        cmd = (
            f"bcftools view --threads {max(1, int(threads))} -R '{regions_file}' '{dbsnp_vcf}' -Ou "
            f"| bcftools query -f '%CHROM\\t%POS\\t%ID\\n'"
        )
        proc = subprocess.run(
            cmd,
            shell=True,
            check=True,
            text=True,
            capture_output=True,
        )

    debug_log(debug, f"Coordinate lookup queried {len(sorted_coords):,} positions.")
    debug_log(debug, f"Coordinate lookup sample positions: {sorted_coords[:5]}")
    debug_log(debug, f"Coordinate lookup raw result lines: {len(proc.stdout.splitlines()):,}")
    if proc.stdout:
        debug_log(debug, f"Coordinate lookup raw stdout sample: {proc.stdout.splitlines()[:5]}")
    if proc.stderr:
        debug_log(debug, f"Coordinate lookup stderr sample: {proc.stderr.splitlines()[:5]}")

    lookup = {}
    for line in proc.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) != 3:
            continue
        chrom, pos, id_field = fields
        chrom = canonicalize_dbsnp_chrom(chrom, chrom_reverse_map or {})
        if chrom is None:
            continue
        try:
            pos_int = int(pos)
        except ValueError:
            continue

        rsid = extract_first_rsid(id_field)
        if rsid is None:
            continue
        key = (chrom, pos_int)
        lookup[key] = select_preferred_rsid(lookup.get(key), rsid)

    if debug and not lookup:
        debug_probe_dbsnp_by_region(
            sorted_coords[:5],
            dbsnp_vcf,
            chrom_query_map=chrom_query_map,
            debug=debug,
            label="Coordinate lookup probe",
        )
        debug_probe_dbsnp_by_window(
            sorted_coords[:5],
            dbsnp_vcf,
            chrom_query_map=chrom_query_map,
            debug=debug,
            label="Coordinate lookup off-by-one probe",
        )

        lookup = build_dbsnp_coord_lookup_by_chrom_scan(
            sorted_coords,
            dbsnp_vcf,
            threads=threads,
            debug=debug,
            chrom_query_map=chrom_query_map,
            chrom_reverse_map=chrom_reverse_map,
        )

    return lookup


def build_dbsnp_coord_lookup_by_chrom_scan(
    sorted_coords,
    dbsnp_vcf,
    threads=1,
    debug=False,
    chrom_query_map=None,
    chrom_reverse_map=None,
):
    if not sorted_coords:
        return {}

    coords_by_chrom = {}
    for chrom, pos in sorted_coords:
        canonical = canonical_chrom_label(chrom)
        coords_by_chrom.setdefault(canonical, set()).add(pos)

    lookup = {}
    for chrom in sorted(coords_by_chrom, key=chrom_to_int):
        dbsnp_chrom = convert_chrom_for_dbsnp(chrom, chrom_query_map or {})
        cmd = (
            f"bcftools view --threads {max(1, int(threads))} -r '{dbsnp_chrom}' '{dbsnp_vcf}' -Ou "
            f"| bcftools query -f '%CHROM\\t%POS\\t%ID\\n'"
        )
        proc = subprocess.run(
            cmd,
            shell=True,
            check=True,
            text=True,
            capture_output=True,
        )

        target_positions = coords_by_chrom[chrom]
        debug_log(
            debug,
            f"Chrom-scan lookup querying {dbsnp_chrom} for {len(target_positions):,} target positions; raw lines={len(proc.stdout.splitlines()):,}",
        )

        for line in proc.stdout.splitlines():
            fields = line.split("\t")
            if len(fields) != 3:
                continue
            raw_chrom, pos, id_field = fields
            canonical_hit_chrom = canonicalize_dbsnp_chrom(raw_chrom, chrom_reverse_map or {})
            if canonical_hit_chrom != chrom:
                continue
            try:
                pos_int = int(pos)
            except ValueError:
                continue
            if pos_int not in target_positions:
                continue
            rsid = extract_first_rsid(id_field)
            if rsid is None:
                continue
            key = (chrom, pos_int)
            lookup[key] = select_preferred_rsid(lookup.get(key), rsid)

        if debug and not any(key[0] == chrom for key in lookup):
            debug_log(debug, f"Chrom-scan lookup found no rsIDs on chromosome {chrom}.")

    debug_log(debug, f"Chrom-scan lookup total coordinate-based rsIDs: {len(lookup):,}")
    return lookup


def write_header(out_fh, raw_header_line, from_assembly, to_assembly, gwas_schema, add_rsid_column=False):
    fields = build_output_header_fields(
        raw_header_line,
        from_assembly,
        to_assembly,
        gwas_schema,
        add_rsid_column=add_rsid_column,
    )
    if fields:
        out_fh.write("\t".join(fields) + "\n")
    else:
        out_fh.write(raw_header_line)


def run_pipeline(
    gwas_file,
    liftover_bed,
    output_file,
    filtered_file,
    from_assembly,
    to_assembly,
    parallel=1,
    chunksize=5000,
    merged_output_file=None,
    dbsnp_vcf=DBSNP_VCF_HG19,
    debug=False,
):
    print(f"{time.ctime()}. Loading liftover BED into memory: {liftover_bed}")
    coord_map = load_liftover_map(liftover_bed)

    gwas_open = gzip.open if gwas_file.endswith(".gz") else open
    out_open = gzip.open if output_file.endswith(".gz") else open
    merge_out_open = gzip.open if (merged_output_file and merged_output_file.endswith(".gz")) else open
    filt_open = gzip.open if filtered_file.endswith(".gz") else open

    n_total = 0
    n_mapped = 0
    header_line = None
    gwas_schema = None
    add_rsid_column = False
    dbsnp_chrom_query_map = None
    dbsnp_chrom_reverse_map = None
    input_rsid_count = 0
    missing_rsid_examples = []

    passed_candidates = []
    filter_reasons = []

    with gwas_open(gwas_file, "rt") as gwas_fh:
        out_merge_fh = None
        if merged_output_file:
            out_merge_fh = merge_out_open(merged_output_file, "wt")

        try:
            for raw_line in gwas_fh:
                if raw_line.startswith("#"):
                    continue
                header_line = raw_line
                gwas_schema = detect_gwas_schema(header_line)
                add_rsid_column = should_add_rsid_column(gwas_schema)
                debug_log(debug, f"Detected GWAS columns ({len(gwas_schema['fields'])}): {gwas_schema['fields']}")
                debug_log(
                    debug,
                    (
                        "Schema indexes: "
                        f"chrom={gwas_schema['chrom_idx']}, pos={gwas_schema['pos_idx']}, "
                        f"a1={gwas_schema['a1_idx']}, a2={gwas_schema['a2_idx']}, "
                        f"snp={gwas_schema['snp_idx']}, id={gwas_schema['id_idx']}"
                    ),
                )
                if "P" not in {normalize_header_name(field) for field in gwas_schema["fields"]}:
                    debug_log(debug, "Input header does not contain a P column after parsing.")
                break

            if header_line is None:
                raise ValueError("No header/data found in input GWAS file.")

            dbsnp_chrom_query_map, dbsnp_chrom_reverse_map = build_dbsnp_chrom_maps(dbsnp_vcf, debug=debug)

            with out_open(output_file, "wt") as out_fh:
                # Prepare worker mode
                if parallel <= 1:
                    init_worker(coord_map, gwas_schema)

                    if out_merge_fh:
                        write_header(out_merge_fh, header_line, from_assembly, to_assembly, gwas_schema)

                    for raw_line in gwas_fh:
                        if raw_line.startswith("#"):
                            continue

                        c_total, c_mapped, merged_fields, to_chr, to_pos, _to_id, rsid, reason = transform_data_line(raw_line)
                        n_total += c_total
                        n_mapped += c_mapped
                        if merged_fields is None:
                            continue
                        if rsid != "NA":
                            input_rsid_count += 1
                        elif len(missing_rsid_examples) < 5:
                            raw_fields = split_gwas_fields(raw_line)
                            snp_idx = gwas_schema.get("snp_idx")
                            if snp_idx is not None and len(raw_fields) > snp_idx:
                                missing_rsid_examples.append(raw_fields[snp_idx])

                        if out_merge_fh:
                            out_merge_fh.write("\t".join(merged_fields) + "\n")
                        if reason:
                            filter_reasons.append((merged_fields, reason, rsid))
                        else:
                            passed_candidates.append((to_chr, to_pos, rsid, merged_fields))
                else:
                    print(
                        f"{time.ctime()}. Parallel merge enabled with {parallel} workers (chunksize={chunksize}).",
                        file=sys.stderr,
                    )
                    if out_merge_fh:
                        write_header(out_merge_fh, header_line, from_assembly, to_assembly, gwas_schema)

                    with Pool(processes=parallel, initializer=init_worker, initargs=(coord_map, gwas_schema)) as pool:
                        batch = []

                        def flush(lines):
                            nonlocal n_total, n_mapped, input_rsid_count
                            if not lines:
                                return
                            for c_total, c_mapped, merged_fields, to_chr, to_pos, _to_id, rsid, reason in pool.map(
                                transform_data_line, lines, chunksize=chunksize
                            ):
                                n_total += c_total
                                n_mapped += c_mapped
                                if merged_fields is None:
                                    continue
                                if rsid != "NA":
                                    input_rsid_count += 1
                                elif len(missing_rsid_examples) < 5:
                                    snp_idx = gwas_schema.get("snp_idx")
                                    merged_snp_idx = 3 + snp_idx if snp_idx is not None else None
                                    if merged_snp_idx is not None and merged_snp_idx < len(merged_fields):
                                        missing_rsid_examples.append(merged_fields[merged_snp_idx])
                                if out_merge_fh:
                                    out_merge_fh.write("\t".join(merged_fields) + "\n")
                                if reason:
                                    filter_reasons.append((merged_fields, reason, rsid))
                                else:
                                    passed_candidates.append((to_chr, to_pos, rsid, merged_fields))

                        for raw_line in gwas_fh:
                            if raw_line.startswith("#"):
                                continue
                            batch.append(raw_line)
                            if len(batch) >= chunksize * 4:
                                flush(batch)
                                batch = []

                        flush(batch)

                rsids = {r for _, _, r, _ in passed_candidates if r != "NA"}
                debug_log(debug, f"Rows carrying an input rsID before dbSNP lookup: {input_rsid_count:,}")
                if missing_rsid_examples:
                    debug_log(debug, f"Example non-rs SNP values: {missing_rsid_examples}")
                print(f"{time.ctime()}. Batch dbSNP query for {len(rsids):,} unique rsIDs...")
                dbsnp_lookup = build_dbsnp_lookup(
                    rsids,
                    dbsnp_vcf,
                    threads=parallel,
                    debug=debug,
                    chrom_reverse_map=dbsnp_chrom_reverse_map,
                )
                print(f"{time.ctime()}. Retrieved {len(dbsnp_lookup):,} chromosomal dbSNP mappings.")

                coord_candidates = {(to_chr, to_pos) for to_chr, to_pos, rsid, _ in passed_candidates if rsid == "NA"}
                if coord_candidates:
                    print(
                        f"{time.ctime()}. Batch dbSNP coordinate query for {len(coord_candidates):,} lifted positions without rsID..."
                    )
                    dbsnp_coord_lookup = build_dbsnp_coord_lookup(
                        coord_candidates,
                        dbsnp_vcf,
                        threads=parallel,
                        debug=debug,
                        chrom_query_map=dbsnp_chrom_query_map,
                        chrom_reverse_map=dbsnp_chrom_reverse_map,
                    )
                    print(f"{time.ctime()}. Retrieved {len(dbsnp_coord_lookup):,} coordinate-based rsIDs.")
                else:
                    dbsnp_coord_lookup = {}

                passed_variants = []
                for to_chr, to_pos, rsid, merged_fields in passed_candidates:
                    resolved_rsid = rsid
                    if rsid != "NA" and rsid in dbsnp_lookup:
                        db_chr, db_pos = dbsnp_lookup[rsid]
                        if db_chr != to_chr or db_pos != to_pos:
                            reason = (
                                f"dbSNP_coordinate_mismatch_dbSNP_{db_chr}_{db_pos}_"
                                f"vs_liftover_{to_chr}_{to_pos}"
                            )
                            filter_reasons.append((merged_fields, reason, resolved_rsid))
                            continue
                    elif rsid == "NA":
                        resolved_rsid = dbsnp_coord_lookup.get((to_chr, to_pos), "NA")

                    output_fields = build_output_row_fields(
                        merged_fields,
                        resolved_rsid,
                        add_rsid_column=add_rsid_column,
                    )
                    passed_variants.append((to_chr, to_pos, "\t".join(output_fields)))

                passed_variants.sort(key=lambda x: (chrom_to_int(x[0]), x[1]))

                # Write passing variants
                write_header(
                    out_fh,
                    header_line,
                    from_assembly,
                    to_assembly,
                    gwas_schema,
                    add_rsid_column=add_rsid_column,
                )
                if debug:
                    output_header = build_output_header_fields(
                        header_line,
                        from_assembly,
                        to_assembly,
                        gwas_schema,
                        add_rsid_column=add_rsid_column,
                    )
                    debug_log(debug, f"Output header columns ({len(output_header)}): {output_header}\n header_line={header_line.strip()}") 
                    debug_log(debug, f"Output header contains P: {'P' in {normalize_header_name(field) for field in output_header}}")
                    debug_log(debug, f"Output header line: {' '.join(output_header)}")
                    if passed_variants:
                        first_output_fields = passed_variants[0][2].split("\t")
                        debug_log(debug, f"First output row column count: {len(first_output_fields)}")
                        debug_log(debug, f"First output row sample: {first_output_fields[:15]}")
                for _chr, _pos, line in passed_variants:
                    out_fh.write(line + "\n")

        finally:
            if out_merge_fh is not None:
                out_merge_fh.close()

    # Write filtered file
    with filt_open(filtered_file, "wt") as filt_fh:
        if header_line is not None:
            hdr = build_output_header_fields(
                header_line,
                from_assembly,
                to_assembly,
                gwas_schema,
                add_rsid_column=add_rsid_column,
            )
            if hdr:
                filt_fh.write("\t".join(hdr) + "\tFILTER_REASON\n")
        for fields, reason, rsid in filter_reasons:
            output_fields = build_output_row_fields(fields, rsid, add_rsid_column=add_rsid_column)
            filt_fh.write("\t".join(output_fields) + f"\t{reason}\n")

    pct = 100.0 * n_mapped / n_total if n_total else 0.0
    print(f"{time.ctime()}. Mapped {n_mapped:,}/{n_total:,} variants ({pct:.1f}%).")
    n_pass = len(passed_variants)
    n_filt = len(filter_reasons)
    print(f"{time.ctime()}. Filtering complete: passed={n_pass:,}, filtered={n_filt:,}")


def main():
    parser = argparse.ArgumentParser(
        description="Merge liftover coordinates and filter with dbSNP validation.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
  %(prog)s -g gwas.tsv.gz -b lift.bed -o pass.tsv.gz -f filtered.tsv.gz \\
           -from hg38 -to hg19 --parallel 4

  Optionally keep the full merged (pre-filter) table:
  %(prog)s -g gwas.tsv.gz -b lift.bed -o pass.tsv.gz -f filtered.tsv.gz \\
           --merged-output-file merged_all.tsv.gz -from hg38 -to hg19
""",
    )
    parser.add_argument("-g", "--gwas-file", required=True, help="Input GWAS summary stats file (can be gzipped)")
    parser.add_argument("-b", "--liftover-bed", required=True, help="Liftover BED file from UCSC (can be gzipped)")
    parser.add_argument("-o", "--output-file", required=True, help="Output file for passing variants (sorted by target coordinates).")
    parser.add_argument("-from", "--from-assembly", required=True, help="Source assembly label (e.g. hg38)")
    parser.add_argument("-to", "--to-assembly", required=True, help="Target assembly label (e.g. hg19)")

    parser.add_argument("-f", "--filtered-file", required=True, help="Output file for filtered variants with FILTER_REASON column.")
    parser.add_argument("--merged-output-file", default=None, help="Optional: also write the full merged (pre-filter) table here.")
    parser.add_argument("-d", "--dbsnp-vcf", default=DBSNP_VCF_HG19, help="dbSNP build 157 GRCh37 VCF path")

    parser.add_argument("--parallel", type=int, default=None, help="Worker count (default: SLURM_CPUS_PER_TASK or 1)")
    parser.add_argument("--chunksize", type=int, default=5000, help="Parallel map chunksize (default: 5000)")
    parser.add_argument("--debug", action="store_true", help="Print extra diagnostics about schema parsing and dbSNP lookup.")

    args = parser.parse_args()

    if args.parallel is None:
        slurm_cpus = os.environ.get("SLURM_CPUS_PER_TASK")
        args.parallel = int(slurm_cpus) if slurm_cpus else 1

    t0 = time.time()
    try:
        run_pipeline(
            gwas_file=args.gwas_file,
            liftover_bed=args.liftover_bed,
            output_file=args.output_file,
            filtered_file=args.filtered_file,
            from_assembly=args.from_assembly,
            to_assembly=args.to_assembly,
            parallel=max(1, args.parallel),
            chunksize=max(100, args.chunksize),
            merged_output_file=args.merged_output_file,
            dbsnp_vcf=args.dbsnp_vcf,
            debug=args.debug,
        )
        dt = time.time() - t0
        print(f"{time.ctime()}. Done in {int(dt // 60)}m {int(dt % 60)}s")
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
