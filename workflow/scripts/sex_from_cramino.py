#!/usr/bin/env python3
import re
import os
import sys
import argparse
import pandas as pd
from typing import Dict, List

def parse_cramino_normalized_counts(path: str) -> Dict[str, float]:
    counts = {}
    in_section = False
    with open(path, 'r', encoding='utf-8') as f:
        for line in f:
            s = line.strip()
            if not in_section:
                if s.startswith('# Normalized read count per chromosome'):
                    in_section = True
                continue
            if not s:
                continue
            if not (s.startswith('chr') or s.startswith('chrUn_') or s.startswith('chrEBV')):
                break
            parts = re.split(r'\s+', s)
            if len(parts) >= 2:
                try:
                    counts[parts[0]] = float(parts[1])
                except ValueError:
                    pass
    return counts

def median_autosome_baseline(counts: Dict[str, float]) -> float:
    autosomes = [f"chr{i}" for i in range(1, 23)]
    vals = [counts[c] for c in autosomes if c in counts]
    if not vals:
        raise ValueError("No autosome counts found.")
    return float(pd.Series(vals).median())

def sample_id_from_path(path: str) -> str:
    base = os.path.basename(path)
    m = re.match(r'([A-Za-z0-9\-]+)_', base)
    if m:
        return m.group(1)
    return os.path.splitext(base)[0]

def classify_karyotype(x_ratio: float, y_ratio: float) -> str:
    y_present = y_ratio >= 0.15
    y_high = y_ratio >= 0.75
    if (0.30 <= x_ratio <= 0.70) and y_present and not y_high:
        return "XY"
    if (0.80 <= x_ratio <= 1.20) and not y_present:
        return "XX"
    if (1.30 <= x_ratio <= 1.70) and y_present and not y_high:
        return "XXY"
    if (0.30 <= x_ratio <= 0.70) and y_high:
        return "XYY"
    if (1.80 <= x_ratio <= 2.20) and not y_present:
        return "XXX"
    if (0.30 <= x_ratio <= 0.70) and not y_present:
        return "XO"
    return "Ambiguous"

def haploid_chroms_from_karyotype(k: str):
    """Return list of haploid chromosomes""" 
    k = k.upper()
    haploids = ["chrM"] # chrM is always haploid
    
    if k == "XX":
        haploids += ["chrY"] # In case of background mapping to chrY
        return haploids
    if k == "XY":
        haploids += ["chrX", "chrY"]
        return haploids
    if k == "XXY":
        haploids += ["chrY"]
        return haploids
    if k == "XYY":
        haploids += ["chrX"]
        return haploids
    if k == "XO":
        haploids += ["chrX"]
        return haploids
    if k == "XXX":
        return haploids
    return haploids

def analyze_files(paths: List[str]) -> pd.DataFrame:
    rows = []
    for p in paths:
        counts = parse_cramino_normalized_counts(p)
        base = sample_id_from_path(p)
        baseline = median_autosome_baseline(counts)
        x = counts.get("chrX", float('nan'))
        y = counts.get("chrY", float('nan'))
        x_ratio = x / baseline
        y_ratio = y / baseline
        k = classify_karyotype(x_ratio, y_ratio)
        haploids = haploid_chroms_from_karyotype(k)
        rows.append({
            "sample_id": base,
            "autosome_median": baseline,
            "chrX_norm": x,
            "chrY_norm": y,
            "X_over_autosomes": x_ratio,
            "Y_over_autosomes": y_ratio,
            "predicted_karyotype": k,
            "haploid_chromosomes": ",".join(haploids) if haploids else "-",
        })
    return pd.DataFrame(rows)

def main():
    ap = argparse.ArgumentParser(description="Predict sex/karyotype from cramino outputs and recommend haploid chromosomes.")
    ap.add_argument("inputs", nargs="+", help="Paths to cramino output .txt files")
    ap.add_argument("-o", "--out", default="cramino_sex_calls.csv", help="Output CSV path")
    args = ap.parse_args()
    df = analyze_files(args.inputs)
    df.to_csv(args.out, index=False)
    print(df.to_string(index=False))

if __name__ == "__main__":
    main()
