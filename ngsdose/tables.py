"""Per-sample summary rows and TSV/JSON helpers shared by the command line and the report."""
from __future__ import annotations

import csv
import gzip
import json
import sys
from pathlib import Path

import numpy as np


class Encoder(json.JSONEncoder):
    def default(self, o):
        if isinstance(o, (np.floating, np.integer)):
            return o.item()
        if isinstance(o, np.ndarray):
            return o.tolist()
        return super().default(o)


def dump(obj, path):
    data = json.dumps(obj, cls=Encoder, allow_nan=True)
    if str(path) == "-":
        sys.stdout.write(data + "\n")
    elif str(path).endswith(".gz"):
        with gzip.open(path, "wt") as fh:
            fh.write(data)
    else:
        Path(path).write_text(data)


def load_result(path) -> dict:
    with (gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)) as fh:
        return json.load(fh)


def summary_row(r: dict) -> dict:
    row = dict(sample=r["sample"], mode=r["mode"], engine=f"{r.get('engine_version', '?')}+{(r.get('engine_build') or 'unknown')[:7]}",
               depth=round(r["depth_equiv"], 3), read_length=r["read_length"],
               insert_median=r["insert_median"], gc_L=r["gc_L"], ctrl_dup_frac=round(r["ctrl_dup_frac"], 4),
               gc_rel_35=r["gc_rel"].get("35"), gc_rel_65=r["gc_rel"].get("65"), gc_curve_max_se=round(r["gc_curve_max_se"], 4),
               ctrl_region_sd=None if not r["control_qc"] else round(r["control_qc"]["region_log_mad_sd"], 4),
               flagged_chromosomes=None if not r["control_qc"] else ",".join(r["control_qc"]["flagged_chromosomes"]))
    for label, t in r.get("truth_regions", {}).items():
        # known-truth sets are scored against their answer; dosage sets (chrM, chrEBV) are copies per cell
        row[f"{label}.copies" if t.get("role") == "dosage" else f"truth.{label}"] = round(t["cn"], 4)
    if r.get("eof_marker"):
        row["eof_marker"] = r["eof_marker"]
    for name, c in r["classes"].items():
        if c["kind"] == "positional":
            row[f"{name}.cn_single"] = round(c["cn"], 3)
            row[f"{name}.cn_anchor"] = round(c["cn_anchor"], 2)
            row[f"{name}.cn_all"] = round(c["cn_all"], 2)
            row[f"{name}.cn_median"] = round(c["cn_median"], 2)
            row[f"{name}.window_log_sd"] = round(c["window_log_sd"], 4)
            row[f"{name}.cn_all_flat"] = round(c["cn_all_flat"], 2)
            row[f"{name}.dup_flag_frac"] = round(c["dup_flag_frac"], 4)
            for fn, fv in c["features"].items():
                row[f"{name}.{fn}"] = round(fv["cn"], 2)
                row[f"{name}.{fn}.flat"] = round(fv["cn_flat"], 2)
        else:
            row[f"{name}.mass_Mb"] = round(c["mass_Mb"], 4)
    return row


def write_table(rows: list[dict], path):
    cols: list[str] = []
    for r in rows:
        cols += [k for k in r if k not in cols]
    fh = sys.stdout if str(path) == "-" else open(path, "w", newline="")
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n", restval="NA")
    w.writeheader()
    for r in rows:
        w.writerow({k: ("NA" if v is None or (isinstance(v, float) and not np.isfinite(v)) else v) for k, v in r.items()})
    if fh is not sys.stdout:
        fh.close()


def read_table(path) -> list[dict]:
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def num(row: dict, col: str) -> float:
    """A cell as a float; NaN for NA, empty or missing."""
    v = row.get(col)
    if v is None or v == "" or v == "NA":
        return float("nan")
    try:
        return float(v)
    except (TypeError, ValueError):
        return float("nan")
