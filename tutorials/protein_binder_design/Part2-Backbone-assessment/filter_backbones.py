#!/usr/bin/env python3
"""
filter_backbones.py

Filter a metrics CSV (one row per structure) by user-defined thresholds and
copy the corresponding structure files to an output directory.

Typical usage:
  python filter_backbones.py \
    --metrics metrics.csv \
    --pdb-root /path/to/pdbs \
    --outdir Passing_Backbones \
    --min rg=9.5 --max rg=12.0 \
    --min pct_helix_total=70 \
    --max pct_C=15 \
    --min nterm_sc=0.9 --max nterm_sc=1.4

Threshold syntax:
  --min <col>=<value>   keep rows where col >= value
  --max <col>=<value>   keep rows where col <= value

You can provide multiple --min / --max flags.

By default, the script uses:
  - `pdb_path` if present and exists
  - otherwise `pdb_root / pdb_file`

It writes:
  - copies of passing structure files to --outdir
  - optionally, a filtered CSV with only passing rows.
"""

import argparse
import csv
import os
import shutil
import sys
from typing import Dict, List, Tuple, Optional


def parse_threshold_kv(s: str) -> Tuple[str, float]:
    """Parse 'col=value' into (col, float(value))."""
    if "=" not in s:
        raise argparse.ArgumentTypeError(f"Threshold must be in the form col=value, got: {s}")
    col, val = s.split("=", 1)
    col = col.strip()
    try:
        fval = float(val.strip())
    except ValueError as e:
        raise argparse.ArgumentTypeError(f"Value must be numeric in {s}") from e
    if not col:
        raise argparse.ArgumentTypeError(f"Column name is empty in {s}")
    return col, fval


def safe_float(x: str) -> Optional[float]:
    """Convert string to float; return None for empty/NA/non-numeric."""
    if x is None:
        return None
    x = str(x).strip()
    if x == "" or x.upper() in {"NA", "NAN", "NONE"}:
        return None
    try:
        return float(x)
    except ValueError:
        return None


def row_passes_thresholds(
    row: Dict[str, str],
    mins: Dict[str, float],
    maxs: Dict[str, float],
    require_ok: bool = False,
) -> Tuple[bool, str]:
    """
    Return (passes, reason_if_failed).
    If a needed value is missing/non-numeric, row fails.
    """
    if require_ok:
        status = (row.get("status") or "").strip().upper()
        if status != "OK":
            return False, f"status != OK ({status or 'missing'})"

    for col, minv in mins.items():
        v = safe_float(row.get(col))
        if v is None:
            return False, f"missing/non-numeric {col}"
        if v < minv:
            return False, f"{col}={v} < min {minv}"

    for col, maxv in maxs.items():
        v = safe_float(row.get(col))
        if v is None:
            return False, f"missing/non-numeric {col}"
        if v > maxv:
            return False, f"{col}={v} > max {maxv}"

    return True, ""


def resolve_source_path(row: Dict[str, str], pdb_root: str) -> Optional[str]:
    """
    Prefer absolute `pdb_path` if it exists. Otherwise use pdb_root/pdb_file.
    """
    pdb_path = (row.get("pdb_path") or "").strip()
    if pdb_path and os.path.isfile(pdb_path):
        return pdb_path

    pdb_file = (row.get("pdb_file") or "").strip()
    if pdb_file and pdb_root:
        candidate = os.path.join(pdb_root, pdb_file)
        if os.path.isfile(candidate):
            return candidate

    return None


def main() -> int:
    parser = argparse.ArgumentParser(
        prog="filter_backbones.py",
        description="Filter a metrics CSV by thresholds and copy passing structure files to a new folder.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--metrics", required=True, help="Path to the metrics CSV.")
    parser.add_argument(
        "--pdb-root",
        default="",
        help="Folder containing the structure files. Used if pdb_path in CSV is missing or not valid.",
    )
    parser.add_argument("--outdir", default="Passing_Backbones", help="Output directory for copied files.")
    parser.add_argument(
        "--min",
        action="append",
        default=[],
        metavar="COL=VAL",
        help="Minimum threshold (keep rows where COL >= VAL). Can be provided multiple times.",
    )
    parser.add_argument(
        "--max",
        action="append",
        default=[],
        metavar="COL=VAL",
        help="Maximum threshold (keep rows where COL <= VAL). Can be provided multiple times.",
    )
    parser.add_argument(
        "--require-ok",
        action="store_true",
        help="Only accept rows with status == OK.",
    )
    parser.add_argument(
        "--copy-mode",
        choices=["copy", "symlink"],
        default="copy",
        help="Copy files or create symlinks in the output directory.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Do not copy anything; just print what would pass/fail.",
    )
    parser.add_argument(
        "--filtered-csv",
        default="",
        help="Optional path to write a CSV containing only passing rows.",
    )

    args = parser.parse_args()

    # Parse thresholds into dicts
    mins: Dict[str, float] = {}
    maxs: Dict[str, float] = {}

    for item in args.min:
        col, val = parse_threshold_kv(item)
        mins[col] = val
    for item in args.max:
        col, val = parse_threshold_kv(item)
        maxs[col] = val

    # Ensure outdir exists (unless dry-run)
    if not args.dry_run:
        os.makedirs(args.outdir, exist_ok=True)

    passed_rows: List[Dict[str, str]] = []
    total = 0
    passed = 0
    missing_files = 0

    with open(args.metrics, newline="") as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames:
            print("ERROR: metrics CSV has no header.", file=sys.stderr)
            return 2

        for row in reader:
            total += 1
            ok, reason = row_passes_thresholds(row, mins, maxs, require_ok=args.require_ok)
            if not ok:
                if args.dry_run:
                    name = row.get("pdb_file") or row.get("pdb_path") or f"row{total}"
                    print(f"FAIL  {name}: {reason}")
                continue

            src = resolve_source_path(row, args.pdb_root)
            if src is None:
                missing_files += 1
                if args.dry_run:
                    name = row.get("pdb_file") or row.get("pdb_path") or f"row{total}"
                    print(f"PASS (but missing file) {name}")
                continue

            dst = os.path.join(args.outdir, os.path.basename(src))

            if args.dry_run:
                print(f"PASS  {os.path.basename(src)}  ->  {dst}")
            else:
                # Avoid overwriting silently if same name appears twice
                if os.path.exists(dst):
                    # Add a suffix
                    base, ext = os.path.splitext(os.path.basename(src))
                    i = 2
                    while True:
                        candidate = os.path.join(args.outdir, f"{base}__{i}{ext}")
                        if not os.path.exists(candidate):
                            dst = candidate
                            break
                        i += 1

                if args.copy_mode == "copy":
                    shutil.copy2(src, dst)
                else:
                    os.symlink(os.path.abspath(src), dst)

            passed += 1
            passed_rows.append(row)

    # Optional: write filtered CSV
    if args.filtered_csv:
        if args.dry_run:
            print(f"(dry-run) Would write filtered CSV to: {args.filtered_csv}")
        else:
            with open(args.filtered_csv, "w", newline="") as out_f:
                writer = csv.DictWriter(out_f, fieldnames=passed_rows[0].keys() if passed_rows else [])
                if passed_rows:
                    writer.writeheader()
                    writer.writerows(passed_rows)

    print("\nSummary")
    print(f"  Total rows read:      {total}")
    print(f"  Passed thresholds:    {len(passed_rows)}")
    print(f"  Files copied/symlink: {passed}")
    print(f"  Missing files:        {missing_files}")
    print(f"  Output directory:     {args.outdir}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())