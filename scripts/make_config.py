#!/usr/bin/env python3
"""
make_config.py — generate and validate a Snakemake config YAML.

Usage examples:

# 1) Create a config from scratch
python make_config.py init \
  --output config.yaml \
  --samples HG001_subset \
  --start-from ubam \
  --input-dir /ifs/data/research/unique/projects/test_config \
  --output-dir /ifs/data/research/unique/projects/test_config \
  --reference /ifs/data/research/unique/leena/references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
  --reference-name GRCh38 \
  --reference-TE /ifs/data/research/unique/repeat-catalogs/TEs/teref.ont.human.fa \
  --tr-catalog /ifs/data/research/unique/repeat-catalogs/TRs/GRCh38/pathogenic/STRchive-disease-loci.v2.2.1.GRCh38.longTR.bed \
  --te-catalog /ifs/data/research/unique/repeat-catalogs/TEs/GRCh38/TE_classes/GRCh38_TEs_retroposon.bed \
  --min-read-quality 7 \
  --min-read-length 500 \
  --flanking-length-bp 250 \
  --extension-repeat-consensus 1000

# 2) Validate an existing config (prints a summary or errors)
python make_config.py validate config.yaml

# 3) Start from a template, edit, then validate
python make_config.py template --output config.yaml
python make_config.py validate config.yaml
"""

import argparse
import sys
import os
import re
from pathlib import Path
from typing import Any, Dict, List, Tuple, Optional

try:
    import yaml
except ImportError as e:
    sys.stderr.write("ERROR: PyYAML not installed. Run: pip install pyyaml\n")
    raise

ALLOWED_START_FROM = {"pod5", "ubam", "fastq", "bam"}
ALLOWED_REFERENCE_NAME = {"chm13v2", "GRCh38"}
ALLOWED_TYPE_OF_TE = {
    "all","all_cpg",
    "DNA","DNA_cpg",
    "LINE","LINE_cpg",
    "helitron","helitron_cpg",
    "SINE","SINE_cpg",
    "LTR","LTR_cpg",
    "retroposon","retroposon_cpg",
}
ALLOWED_TYPE_OF_TR = {"genome-wide", "pathogenic", "forensic"}

def _file_exists(p: str) -> bool:
    return p and Path(p).is_file()

def _dir_exists(p: str) -> bool:
    return p and Path(p).is_dir()

def _missing_reference_indexes(reference: str) -> List[str]:
    """
    Return a list of expected index files that are missing for the given reference FASTA.

    Adjust this list depending on which tools you actually use
    (e.g., minimap2 .mmi, samtools .fai only, etc.).
    """
    if not reference:
        return []

    # Base path is the full FASTA path; append suffixes for different index files.
    # Example: /path/genome.fa -> /path/genome.fa.fai, /path/genome.fa.bwt, ...
    expected_suffixes = [
        ".fai"]  # samtools/htslib FASTA index

    missing = []
    for suf in expected_suffixes:
        idx_path = reference + suf
        if not _file_exists(idx_path):
            missing.append(idx_path)
    return missing

def _as_int(name: str, v: Any, min_val: Optional[int] = None) -> Tuple[Optional[int], Optional[str]]:
    try:
        iv = int(v)
        if min_val is not None and iv < min_val:
            return None, f"'{name}' must be ≥ {min_val}, got {iv}"
        return iv, None
    except Exception:
        return None, f"'{name}' must be an integer, got {v!r}"

def _comma_or_list_to_list(values: Any) -> List[str]:
    # Accept list, space-separated values, or comma-separated string
    if values is None:
        return []
    if isinstance(values, list):
        items = values
    else:
        s = str(values)
        if "," in s:
            items = [x.strip() for x in s.split(",")]
        else:
            items = [x.strip() for x in s.split()]
    return [x for x in items if x]

def _infer_te_type_from_path(te_catalog: str) -> Optional[str]:
    """
    Infer TE type from filename.
    Expected pattern: ... TE_XXX.bed or TEs_XXX.bed; capture XXX (may include _cpg).
    """
    if not te_catalog:
        return None
    name = Path(te_catalog).name
    # Try 'TE_' first
    m = re.search(r"(?:^|[_-])TE_([^./]+)\.bed$", name)
    if not m:
        # Support 'TEs_'
        m = re.search(r"(?:^|[_-])TEs_([^./]+)\.bed$", name)
    if m:
        return m.group(1)
    return None

def _infer_tr_type_from_path(tr_catalog: str) -> Optional[str]:
    """
    Infer TR type by scanning the full path for one of the allowed strings.
    """
    if not tr_catalog:
        return None
    s = str(tr_catalog).lower()
    for k in ("genome-wide", "pathogenic", "forensic"):
        if k in s:
            return k
    return None

def _validate_config(cfg: Dict[str, Any]) -> List[str]:
    errors: List[str] = []

    # samples
    samples = cfg.get("samples")
    if not isinstance(samples, list) or not samples:
        errors.append("Field 'samples' must be a non-empty list of sample IDs.")
    else:
        bad = [s for s in samples if not isinstance(s, str) or not s.strip()]
        if bad:
            errors.append(f"'samples' contains empty/invalid entries: {bad}")

    # start_from
    start_from = cfg.get("start_from")
    if start_from not in ALLOWED_START_FROM:
        errors.append(f"'start_from' must be one of {sorted(ALLOWED_START_FROM)}, got {start_from!r}")

    # input_dir / output_dir
    input_dir = cfg.get("input_dir")
    output_dir = cfg.get("output_dir")
    if not _dir_exists(input_dir):
        errors.append(f"'input_dir' does not exist or is not a directory: {input_dir!r}")
    if not _dir_exists(output_dir):
        errors.append(f"'output_dir' does not exist or is not a directory: {output_dir!r}")

    # reference files
    reference = cfg.get("reference")
    reference_name = cfg.get("reference_name")
    reference_TE = cfg.get("reference_TE")

    if not _file_exists(reference):
        errors.append(f"'reference' file not found: {reference!r}")
    else:
        missing_idx = _missing_reference_indexes(reference)
        if missing_idx:
            print(
                "WARNING: some reference index files are missing:\n  " +
                "\n  ".join(missing_idx),
                file=sys.stderr,
            )

    if reference_name not in ALLOWED_REFERENCE_NAME:
        errors.append(
            f"'reference_name' must be one of {sorted(ALLOWED_REFERENCE_NAME)}, "
            f"got {reference_name!r}"
        )

    if not _file_exists(reference_TE):
        errors.append(f"'reference_TE' file not found: {reference_TE!r}")    


    # catalogs (optional but validated if provided)
    tr_catalog = cfg.get("tr_catalog")
    te_catalog = cfg.get("te_catalog")
    if tr_catalog and not _file_exists(tr_catalog):
        errors.append(f"'tr_catalog' file not found: {tr_catalog!r}")
    if te_catalog and not _file_exists(te_catalog):
        errors.append(f"'te_catalog' file not found: {te_catalog!r}")

    # types (must be present and allowed)
    type_of_te = cfg.get("type_of_te")
    if type_of_te not in ALLOWED_TYPE_OF_TE:
        errors.append(f"'type_of_te' must be one of {sorted(ALLOWED_TYPE_OF_TE)}, got {type_of_te!r}")

    type_of_tr = cfg.get("type_of_tr")
    if type_of_tr not in ALLOWED_TYPE_OF_TR:
        errors.append(f"'type_of_tr' must be one of {sorted(ALLOWED_TYPE_OF_TR)}, got {type_of_tr!r}")

    # fastq_filtering
    fq = cfg.get("fastq_filtering", {})
    if not isinstance(fq, dict):
        errors.append("'fastq_filtering' must be a mapping with 'min_read_quality' and 'min_read_length'.")
    else:
        q, e1 = _as_int("fastq_filtering.min_read_quality", fq.get("min_read_quality"), min_val=0)
        l, e2 = _as_int("fastq_filtering.min_read_length", fq.get("min_read_length"), min_val=0)
        if e1: errors.append(e1)
        if e2: errors.append(e2)

    # flanking_length_bp / extension_repeat_consensus
    flanking, e3 = _as_int("flanking_length_bp", cfg.get("flanking_length_bp"), min_val=0)
    ext, e4 = _as_int("extension_repeat_consensus", cfg.get("extension_repeat_consensus"), min_val=0)
    if e3: errors.append(e3)
    if e4: errors.append(e4)

    return errors

def _load_yaml(path: Path) -> Dict[str, Any]:
    with path.open("r") as fh:
        data = yaml.safe_load(fh) or {}
    if not isinstance(data, dict):
        raise ValueError("YAML root must be a mapping.")
    return data

def _dump_yaml(path: Path, data: Dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        yaml.safe_dump(data, fh, sort_keys=False)

def _template_config() -> Dict[str, Any]:
    return {
        "samples": ["HG001_subset"],
        "start_from": "ubam",  # Options: pod5, ubam, fastq, bam
        "input_dir": "/path/to/input_dir",
        "output_dir": "/path/to/output_dir",
        "reference": "/path/to/GRCh38.fasta",
        "reference_name": "GRCh38",  # or chm13v2
        "reference_TE": "/path/to/teref.ont.human.fa",
        "fastq_filtering": {
            "min_read_quality": 7,
            "min_read_length": 500,
        },
        "tr_catalog": "/path/to/STRchive-disease-loci.v2.2.1.GRCh38.longTR.bed",
        "te_catalog": "/path/to/GRCh38_TEs_retroposon.bed",
        # auto-filled during init; kept here to show expected fields:
        "type_of_te": "retroposon",
        "type_of_tr": "pathogenic",
        "flanking_length_bp": 250,
        "extension_repeat_consensus": 1000,
    }

def cmd_template(args: argparse.Namespace) -> int:
    data = _template_config()
    out = Path(args.output)
    _dump_yaml(out, data)
    print(f"Template written: {out}")
    return 0

def cmd_validate(args: argparse.Namespace) -> int:
    path = Path(args.config)
    if not path.exists():
        print(f"ERROR: Config not found: {path}", file=sys.stderr)
        return 2
    cfg = _load_yaml(path)
    errors = _validate_config(cfg)
    if errors:
        print("CONFIG INVALID:")
        for e in errors:
            print(f" - {e}")
        return 1
    else:
        print("CONFIG OK ✅")
        # Optional: echo a brief normalized view of critical paths
        for k in ("input_dir","output_dir","reference","reference_TE","tr_catalog","te_catalog"):
            if k in cfg and cfg[k]:
                print(f"{k}: {cfg[k]}")
        return 0

def cmd_init(args: argparse.Namespace) -> int:
    # Build dict from args
    cfg: Dict[str, Any] = {}
    cfg["samples"] = args.samples or _comma_or_list_to_list(args.samples_str)
    if args.samples_file:
        with open(args.samples_file, "r") as fh:
            file_samples = [line.strip() for line in fh if line.strip() and not line.startswith("#")]
        cfg["samples"].extend(file_samples)
    cfg["samples"] = sorted(set(cfg["samples"]))

    cfg["start_from"] = args.start_from
    cfg["input_dir"] = args.input_dir
    cfg["output_dir"] = args.output_dir
    cfg["reference"] = args.reference
    cfg["reference_name"] = args.reference_name
    cfg["reference_TE"] = args.reference_TE
    cfg["fastq_filtering"] = {
        "min_read_quality": args.min_read_quality,
        "min_read_length": args.min_read_length,
    }
    cfg["tr_catalog"] = args.tr_catalog
    cfg["te_catalog"] = args.te_catalog

    # Infer types
    inferred_te = _infer_te_type_from_path(args.te_catalog)
    if not inferred_te:
        print("ERROR: Could not infer 'type_of_te' from TE catalog filename. "
              "Expected pattern like '*TE_<TYPE>.bed' or '*TEs_<TYPE>.bed'.", file=sys.stderr)
        return 1
    if inferred_te not in ALLOWED_TYPE_OF_TE:
        print(f"ERROR: Inferred 'type_of_te'='{inferred_te}' is not one of allowed values: "
              f"{sorted(ALLOWED_TYPE_OF_TE)}", file=sys.stderr)
        return 1
    cfg["type_of_te"] = inferred_te

    inferred_tr = _infer_tr_type_from_path(args.tr_catalog)
    if not inferred_tr:
        print("ERROR: Could not infer 'type_of_tr' from TR catalog path. "
              "None of the patterns ['genome-wide','pathogenic','forensic'] were found in the path.", file=sys.stderr)
        return 1
    cfg["type_of_tr"] = inferred_tr

    cfg["flanking_length_bp"] = args.flanking_length_bp
    cfg["extension_repeat_consensus"] = args.extension_repeat_consensus

    errors = _validate_config(cfg)
    if errors:
        print("Refusing to write invalid config. Found:", file=sys.stderr)
        for e in errors:
            print(f" - {e}", file=sys.stderr)
        return 1

    out = Path(args.output)
    _dump_yaml(out, cfg)
    print(f"Wrote valid config to: {out}")
    return 0

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Create and validate Snakemake config YAML.")
    sub = p.add_subparsers(dest="cmd", required=True)

    # template
    p_t = sub.add_parser("template", help="Write a filled-out template you can edit.")
    p_t.add_argument("--output", required=True, help="Path to write the template YAML.")
    p_t.set_defaults(func=cmd_template)

    # validate
    p_v = sub.add_parser("validate", help="Validate an existing config YAML.")
    p_v.add_argument("config", help="Path to config YAML.")
    p_v.set_defaults(func=cmd_validate)

    # init
    p_i = sub.add_parser("init", help="Create a config from CLI arguments.")
    p_i.add_argument("--output", required=True, help="Where to write the YAML.")

    # samples options
    g_s = p_i.add_argument_group("Samples")
    g_s.add_argument("--samples", nargs="+", default=[],
                     help="Sample IDs (space-separated).")
    g_s.add_argument("--samples-file", help="File with one sample ID per line (supports # comments).")
    g_s.add_argument("--samples-str", default="", help="Samples as a comma- or space-separated string.")

    # core inputs
    g_c = p_i.add_argument_group("Core paths and references")
    g_c.add_argument("--start-from", choices=sorted(ALLOWED_START_FROM), required=True)
    g_c.add_argument("--input-dir", required=True)
    g_c.add_argument("--output-dir", required=True)
    g_c.add_argument("--reference", required=True)
    g_c.add_argument("--reference-name", choices=sorted(ALLOWED_REFERENCE_NAME), required=True)
    g_c.add_argument("--reference-TE", required=True, dest="reference_TE")

    # catalogs
    g_cat = p_i.add_argument_group("Catalogs")
    g_cat.add_argument("--tr-catalog", required=True, dest="tr_catalog")
    g_cat.add_argument("--te-catalog", required=True, dest="te_catalog")

    # filtering + params
    g_f = p_i.add_argument_group("Read filtering")
    g_f.add_argument("--min-read-quality", type=int, default=7)
    g_f.add_argument("--min-read-length", type=int, default=500)

    g_p = p_i.add_argument_group("Analysis parameters")
    g_p.add_argument("--flanking-length-bp", type=int, default=250, dest="flanking_length_bp")
    g_p.add_argument("--extension-repeat-consensus", type=int, default=1000, dest="extension_repeat_consensus")

    p_i.set_defaults(func=cmd_init)

    return p

def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)

if __name__ == "__main__":
    raise SystemExit(main())

