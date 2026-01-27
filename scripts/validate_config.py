# validate_config.py
from pathlib import Path

def _file_exists(p):
    return p and Path(p).is_file()

def _dir_exists(p):
    return p and Path(p).is_dir()

def validate_config(cfg):
    errors = []

    for k in ["samples", "input_dir", "output_dir", "reference"]:
        if k not in cfg:
            errors.append(f"Missing required field: {k}")

    if not _dir_exists(cfg["input_dir"]):
        errors.append(f"input_dir not found: {cfg['input_dir']}")

    if not _file_exists(cfg["reference"]):
        errors.append(f"reference FASTA not found: {cfg['reference']}")

    for f in ["tr_catalog", "te_catalog", "reference_TE"]:
        if not _file_exists(cfg.get(f)):
            errors.append(f"{f} file not found: {cfg.get(f)}")

    tm = cfg.get("tr_methylation", {})
    if tm.get("filter_to_cpg_str"):
        if not _file_exists(tm.get("cpg_str_bed")):
            errors.append("CpG STR BED missing or invalid")

    return errors

