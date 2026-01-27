#!/usr/bin/env python3
"""
make_config.py — generate and validate a Snakemake config YAML.

Usage examples:

# 1) Create a config from scratch
python make_config.py init \
  --output <config-name>.yaml \				# Name of the config file to create
  --samples <sample1> <sample2> .. \			# list of sampleIDs (space or comma-separated)
  --start-from <pod5|ubam|fastq|bam> \       		# Start from a specific input type (choose one)
  --input-dir <path-to-project-input-dir> \		# Path to the project input directory
  --output-dir <path-to-project-output-dir> \           # Path to the project output directory
  --reference <path to reference fasta> \		# path to human reference genome fasta file (should be chr naming)
  --reference-name <GRCh38|chm13v2> \			# Genome build (e.g., GRCh38 or chm13v2), default = GRCh38
  --tr-catalog <path-to-tr-catalog> \			# Optional: Path to the TR catalog (if not default GRCh38 genome-wide adotto catalogue is used)
  --te-catalog <path-to-te-catalog> \			# Optional: Path to the TR catalog (if not default GRCh38 genome-wide TE (all) catalog is used)
  --min-read-quality <value> \				# Minimum read quality (default: 7)
  --min-read-length <value>  \				# Minimum read length (default: 500)
  --flanking-length-bp <value> \			# Flanking length in base pairs (default: 250)
  --extension-repeat-consensus <value> \		# Repeat consensus extension length (default: 1000)

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
    "all",
    "DNA",
    "LINE",
    "helitron",
    "SINE",
    "LTR",
    "retroposon",
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
    Infer TE type from filename
    """
    if not te_catalog:
        return None
    name = Path(te_catalog).name
    
    m = re.search(r"_TEs_([^./]+)\.bed$", name)
    if m:
        return m.group(1)

    m = re.search(r"_TE_([^./]+)\.bed$", name)
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

def _bundled_cpg_str_bed(db_root: str, build: str) -> str:
    """
    Return the default CpG-containing STR BED for the bundled echoDB.
    """
    if build == "GRCh38":
        return str(Path(db_root) / "TRs/GRCh38/genome-wide-str-cpg/adotto_longTR_STR_cpgmotif.bed")
    if build == "chm13v2":
        raise ValueError("No bundled CpG STR BED available for chm13v2 in echoDB_v1.")
    raise ValueError(f"Unknown build for CpG STR BED: {build!r}")


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

    reference_TE = cfg.get("reference_TE")
    if not reference_TE or not _file_exists(reference_TE):
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

    # tr_methylation (CpG STR filtering)
    tm = cfg.get("tr_methylation", {})
    if tm and not isinstance(tm, dict):
        errors.append("'tr_methylation' must be a mapping.")
    elif isinstance(tm, dict):
        filt = bool(tm.get("filter_to_cpg_str", False))
        bed = tm.get("cpg_str_bed", None)
        if filt:
            if not bed:
                errors.append("tr_methylation.filter_to_cpg_str is true but tr_methylation.cpg_str_bed is missing.")
            elif not _file_exists(bed):
                errors.append(f"'tr_methylation.cpg_str_bed' file not found: {bed!r}")

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
    db_root = args.db_root if getattr(args, "use_bundled_db", False) else "resources/echoDB_v1"
    cfg["reference_TE"] = str(Path(db_root) / "TEs/teref.ont.human.fa")
    cfg["fastq_filtering"] = {
        "min_read_quality": args.min_read_quality,
        "min_read_length": args.min_read_length,
    }

    # --- Catalog selection: bundled DB mode OR explicit paths ---
    cfg["tr_catalog"] = args.tr_catalog
    cfg["te_catalog"] = args.te_catalog

    if args.use_bundled_db: 
        # Bundled DB mode: don't require tr_catalog or te_catalog

        db_root = args.db_root
        build = args.reference_name # "GRCh38" or "chm13v2"

        # TR catalog: default is Adotto for genome-wide TR genotyping
        if build == "GRCh38":
            if args.tr_type == "genome-wide":
                cfg["tr_catalog"] = str(Path(db_root) / "TRs/GRCh38/genome-wide/adotto_longTR.bed")
            elif args.tr_type == "pathogenic":
                cfg["tr_catalog"] = str(Path(db_root) / "TRs/GRCh38/pathogenic/STRchive-disease-loci.v2.2.1.GRCh38.longTR.bed")
            elif args.tr_type == "forensic":
                cfg["tr_catalog"] = str(Path(db_root) / "TRs/GRCh38/forensic/STRbase_GRCh38_STRloci_longTR.bed")
            else:
                print(f"ERROR: Unsupported --tr-type {args.tr_type!r} for GRCh38.", file=sys.stderr)
                return 1
        elif build == "chm13v2":
            if args.tr_type == "pathogenic":
                cfg["tr_catalog"] = str(Path(db_root) / "TRs/T2T-CHM13v2/pathogenic/STRchive-disease-loci.v2.2.1.T2T-CHM13.longTR.bed")
            else:
                print("ERROR: For chm13v2, only pathogenic TR catalog is available in echoDB_v1.", file=sys.stderr)
                return 1
        else:
            print(f"ERROR: Unknown reference_name/build: {build!r}", file=sys.stderr)
            return 1

        # TE catalog:
        te_base = args.te_type or "all"
        if build == "GRCh38":
            if te_base == "all":
                cfg["te_catalog"] = str(Path(db_root) / "TEs/GRCh38/GRCh38_TEs_all.bed")
            else:
                cfg["te_catalog"] = str(Path(db_root) / f"TEs/GRCh38/TE_classes/GRCh38_TEs_{te_base}.bed")
        elif build == "chm13v2":
            if te_base == "all":
                cfg["te_catalog"] = str(Path(db_root) / "TEs/T2T-CHM13v2/T2T-CHM13_TEs_all.bed")
            else:
                cfg["te_catalog"] = str(Path(db_root) / f"TEs/T2T-CHM13v2/TE_classes/T2T-CHM13_TEs_{te_base}.bed")

        # Set types directly (no inference needed)
        cfg["type_of_tr"] = args.tr_type
        cfg["type_of_te"] = args.te_type or "all"

        # Record bundled-db intent in config (so Snakemake can default later)
        cfg["catalog_defaults"] = {
            "enabled": True,
            "db_root": db_root,
            "build": build,
        }

    else:
        # Manual mode: require explicit paths and keep legacy inference
        if not args.tr_catalog or not args.te_catalog:
            print("ERROR: Provide --tr-catalog and --te-catalog (or use --use-bundled-db).", file=sys.stderr)
            return 1

        # Infer types
        inferred_te = _infer_te_type_from_path(args.te_catalog)
        if not inferred_te:
            print("ERROR: Could not infer 'type_of_te' from TE catalog filename.", file=sys.stderr)
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

    # --- TR methylation CpG STR filtering (defaults depend on bundled vs manual mode) ---

    cfg["tr_methylation"] = {}
    cfg_cpg_bed = args.cpg_str_bed  # may be None

    if args.use_bundled_db:
        # Bundled DB: default ON unless user explicitly disables
        cpg_filter = True if args.cpg_filter is None else bool(args.cpg_filter)

        if cpg_filter:
            if cfg_cpg_bed is None:
                try:
                    cfg_cpg_bed = _bundled_cpg_str_bed(args.db_root, args.reference_name)
                except Exception as e:
                    print(f"ERROR: CpG filtering requested but no bundled CpG STR BED available: {e}", file=sys.stderr)
                    return 1
            cfg["tr_methylation"]["filter_to_cpg_str"] = True
            cfg["tr_methylation"]["cpg_str_bed"] = cfg_cpg_bed
        else:
            cfg["tr_methylation"]["filter_to_cpg_str"] = False
            cfg["tr_methylation"]["cpg_str_bed"] = None

    else:
        # Manual catalogs: default OFF to avoid mismatches unless explicitly enabled
        cpg_filter = False if args.cpg_filter is None else bool(args.cpg_filter)

        if cpg_filter:
            if cfg_cpg_bed is None:
                print("ERROR: --cpg-filter enabled with custom TR catalogs, but --cpg-str-bed was not provided. "
                      "Provide a CpG STR BED matched to your TR catalog, or disable filtering.", file=sys.stderr)
                return 1
            cfg["tr_methylation"]["filter_to_cpg_str"] = True
            cfg["tr_methylation"]["cpg_str_bed"] = cfg_cpg_bed
        else:
            cfg["tr_methylation"]["filter_to_cpg_str"] = False
            cfg["tr_methylation"]["cpg_str_bed"] = None

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
    p_t.add_argument("--output", required=True, help="Path and filename to write the template YAML, i.e config.yaml")
    p_t.set_defaults(func=cmd_template)

    # validate
    p_v = sub.add_parser("validate", help="Validate an existing config YAML.")
    p_v.add_argument("config", help="Path to config YAML.")
    p_v.set_defaults(func=cmd_validate)

    # init
    p_i = sub.add_parser("init", help="Create a config from CLI arguments.")
    p_i.add_argument("--output", required=True, help="Path and filename to write config YAML.")

    # samples options
    g_s = p_i.add_argument_group("Samples")
    g_s.add_argument("--samples", nargs="+", default=[],
                     help="Sample IDs (space-separated).")
    g_s.add_argument("--samples-file", help="File with one sample ID per line (supports # comments).")
    g_s.add_argument("--samples-str", default="", help="Samples as a comma- or space-separated string.")

    # core inputs
    g_c = p_i.add_argument_group("Core paths and references")
    g_c.add_argument("--start-from", choices=sorted(ALLOWED_START_FROM), required=True, help="Starting point for processing <pod5|ubam|fastq|bam>. If pre-basecalled, needs to be with methylation-aware basecalling model" )
    g_c.add_argument("--input-dir", required=True, help="Project directory for input data")
    g_c.add_argument("--output-dir", required=True, help="Project directory for output data")
    g_c.add_argument("--reference", required=True, help="Path to the human reference genome FASTA file")
    g_c.add_argument("--reference-name", choices=sorted(ALLOWED_REFERENCE_NAME), default="GRCh38", help="Reference genome build <GRCh38|T2T-CHM13v2>(default: GRCh38)")

    # catalogs
    g_cat = p_i.add_argument_group("Catalogs")
    g_cat.add_argument("--tr-catalog", required=False, dest="tr_catalog",
                       help="Path to TR catalog BED (optional if using --use-bundled-db).")
    g_cat.add_argument("--te-catalog", required=False, dest="te_catalog",
                       help="Path to TE catalog BED (optional if using --use-bundled-db).")
    g_cat.add_argument("--use-bundled-db", dest="use_bundled_db", action="store_true",
                       help="Use the ECHO repeat database installed in resources (no need to pass catalog paths).")
    g_cat.add_argument("--custom-db", dest="use_bundled_db", action="store_false",
                       help="Do not use echoDB; require explicit --tr-catalog and --te-catalog.")
    g_cat.set_defaults(use_bundled_db=True)
    g_cat.add_argument("--db-root", default="resources/echoDB_v1",
                       help="Root folder of the bundled ECHO repeat databases (default: resources/echoDB_v1).")
    g_cat.add_argument("--tr-type", choices=sorted(ALLOWED_TYPE_OF_TR), default="genome-wide",
                       help="ECHO TR catalog type to use with --use-bundled-db (genome-wide/pathogenic/forensic).")
    g_cat.add_argument("--te-type", choices=sorted(ALLOWED_TYPE_OF_TE), default=all,
                       help="ECHO TE catalog type to use with --use-bundled-db (all/LINE/SINE/...).")
    g_cat.add_argument("--cpg-filter", dest="cpg_filter", action="store_true",
                       help="Enable filtering to CpG-containing STR loci before TR methylation profiling.")
    g_cat.add_argument("--no-cpg-filter", dest="cpg_filter", action="store_false",
                       help="Disable CpG STR filtering before TR methylation profiling.")
    g_cat.set_defaults(cpg_filter=None)  # None means choose default based on bundled/custom mode
    g_cat.add_argument("--cpg-str-bed", dest="cpg_str_bed", default=None,
                       help="BED file of CpG-containing STR loci to filter to (required if --cpg-filter with custom TR catalogs).")

    # filtering + params
    g_f = p_i.add_argument_group("Read filtering")
    g_f.add_argument("--min-read-quality", type=int, default=7, help="Minimum read quality (default: 7)")
    g_f.add_argument("--min-read-length", type=int, default=500, help="Minimum read length (default: 500 bp)")

    g_p = p_i.add_argument_group("Analysis parameters")
    g_p.add_argument("--flanking-length-bp", type=int, default=250, dest="flanking_length_bp", help="Flanking length (bp) of repeats used for re-alignment (default: 250)")
    g_p.add_argument("--extension-repeat-consensus", type=int, default=1000, dest="extension_repeat_consensus", help="Length (bp) of TE consensus extension (default: 1000)")

    p_i.set_defaults(func=cmd_init)

    return p

def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)

if __name__ == "__main__":
    raise SystemExit(main())

