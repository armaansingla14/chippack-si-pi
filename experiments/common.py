"""Shared helpers for experiment scripts."""
from __future__ import annotations

import csv
import os
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RESULTS = os.path.join(ROOT, "results")
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)


def write_csv(name, rows):
    os.makedirs(RESULTS, exist_ok=True)
    with open(os.path.join(RESULTS, name), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


def md_table(rows, cols=None, fmt=None):
    cols = cols or list(rows[0].keys())
    fmt = fmt or {}
    out = ["| " + " | ".join(cols) + " |", "|" + "|".join("---" for _ in cols) + "|"]
    for r in rows:
        cells = []
        for c in cols:
            v = r[c]
            if isinstance(v, float):
                v = format(v, fmt.get(c, ".4g"))
            elif v is None:
                v = "-"
            cells.append(str(v))
        out.append("| " + " | ".join(cells) + " |")
    return "\n".join(out)


def write_md(name, text):
    os.makedirs(RESULTS, exist_ok=True)
    with open(os.path.join(RESULTS, name), "w") as f:
        f.write(text.rstrip() + "\n")
