#!/usr/bin/env python3
"""Analyze backbone metrics for a directory of PDB files.

Computes, per PDB:
1) N-term and C-term normalized end scores against an interface point
2) Radius of gyration (Rg) for a chain + atom selection
3) STRIDE secondary-structure percentages

Example:
python analyze_backbones.py \
  --pdb-dir ./pdbs \
  --binder-chain A \
  --interface-point 0 0 -16.17 \
  --stride /path/to/stride \
  --atoms alpha-carbons \
  --output metrics.csv \
  --threads 8
"""

from __future__ import annotations

import argparse
import csv
import math
import shutil
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np


STRIDE_CLASSES: Tuple[str, ...] = ("H", "E", "C", "T", "G", "I", "B")
ATOM_SELECTIONS: Dict[str, Optional[set[str]]] = {
    "alpha-carbons": {"CA"},
    "backbone": {"N", "CA", "C", "O"},
    "backbone-no-carbonyl": {"N", "CA", "C"},
    "all": None,  # All ATOM records in the selected chain
}

OUTPUT_COLUMNS: Tuple[str, ...] = (
    "pdb_file",
    "pdb_path",
    "n_atoms_rg",
    "rg",
    "n_ca_binder",
    "nterm_sc",
    "cterm_sc",
    "pct_H",
    "pct_E",
    "pct_C",
    "pct_T",
    "pct_G",
    "pct_I",
    "pct_B",
    "pct_helix_total",
    "pct_extended_total",
    "stride_total_residues",
    "status",
    "error_message",
)


class StrideFailureError(RuntimeError):
    """Raised when STRIDE fails and --fail-on-missing-stride is enabled."""


def parse_pdb_atoms(pdb_path: Path) -> List[Tuple[str, str, str, float, float, float]]:
    """Parse ATOM records from a PDB file.

    Returns tuples:
    (atom_name, chain_id, residue_index, x, y, z)
    """
    atoms: List[Tuple[str, str, str, float, float, float]] = []
    with pdb_path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith("ATOM"):
                continue
            if len(line) < 54:
                continue

            atom_name = line[12:16].strip()
            chain_id = line[21].strip()
            residue_index = line[22:26].strip()

            try:
                x_coord = float(line[30:38])
                y_coord = float(line[38:46])
                z_coord = float(line[46:54])
            except ValueError:
                parts = line.split()
                if len(parts) < 9:
                    continue
                try:
                    atom_name = parts[2]
                    chain_id = parts[4]
                    x_coord = float(parts[6])
                    y_coord = float(parts[7])
                    z_coord = float(parts[8])
                except (ValueError, IndexError):
                    continue

            atoms.append((atom_name, chain_id, residue_index, x_coord, y_coord, z_coord))
    return atoms


def compute_end_scores(
    binder_ca_coords: np.ndarray, interface_point: np.ndarray
) -> Tuple[int, Optional[float], Optional[float], Optional[str]]:
    """Compute normalized N-term and C-term scores."""
    n_ca = int(binder_ca_coords.shape[0])
    if n_ca == 0:
        return 0, None, None, "No binder-chain CA atoms found for end scores."

    deltas = binder_ca_coords - interface_point
    squared_distances = np.sum(deltas * deltas, axis=1)
    distances = np.sqrt(squared_distances)
    rmsd_to_point = float(math.sqrt(float(np.mean(squared_distances))))

    if rmsd_to_point == 0.0:
        return n_ca, None, None, "RMS distance to interface point is 0. End scores are undefined."

    nterm_sc = float(distances[0] / rmsd_to_point)
    cterm_sc = float(distances[-1] / rmsd_to_point)
    return n_ca, nterm_sc, cterm_sc, None


def compute_rg(rg_coords: np.ndarray) -> Tuple[int, Optional[float], Optional[str]]:
    """Compute radius of gyration using geometric center."""
    n_atoms = int(rg_coords.shape[0])
    if n_atoms == 0:
        return 0, None, "No atoms found for requested Rg chain/atom selection."

    center = np.mean(rg_coords, axis=0)
    squared_distances = np.sum((rg_coords - center) ** 2, axis=1)
    rg = float(math.sqrt(float(np.mean(squared_distances))))
    return n_atoms, rg, None


def empty_stride_metrics() -> Dict[str, Optional[float]]:
    metrics: Dict[str, Optional[float]] = {f"pct_{ss_class}": None for ss_class in STRIDE_CLASSES}
    metrics["pct_helix_total"] = None
    metrics["pct_extended_total"] = None
    metrics["stride_total_residues"] = None
    return metrics


def run_stride_and_parse(
    stride_executable: str, pdb_path: Path, timeout_seconds: int = 30
) -> Tuple[Dict[str, Optional[float]], Optional[str]]:
    """Run STRIDE and parse secondary-structure percentages."""
    command = [stride_executable, "-f", str(pdb_path)]
    try:
        completed = subprocess.run(
            command,
            capture_output=True,
            text=True,
            timeout=timeout_seconds,
            check=False,
        )
    except FileNotFoundError as exc:
        raise RuntimeError(f"STRIDE executable not found: {stride_executable}") from exc
    except subprocess.TimeoutExpired as exc:
        raise RuntimeError(f"STRIDE timed out after {timeout_seconds}s.") from exc

    if completed.returncode != 0:
        stderr = completed.stderr.strip()
        stdout = completed.stdout.strip()
        details = stderr or stdout or f"return code {completed.returncode}"
        raise RuntimeError(f"STRIDE failed ({details}).")

    counts = {ss_class: 0 for ss_class in STRIDE_CLASSES}
    unknown_classes: Dict[str, int] = {}
    total_residues = 0

    for line in completed.stdout.splitlines():
        if not line.startswith("ASG"):
            continue
        columns = line.split()
        if len(columns) <= 5:
            continue
        ss_class = columns[5].upper()
        total_residues += 1
        if ss_class in counts:
            counts[ss_class] += 1
        else:
            unknown_classes[ss_class] = unknown_classes.get(ss_class, 0) + 1

    if total_residues == 0:
        metrics = empty_stride_metrics()
        metrics["stride_total_residues"] = 0
        return metrics, "STRIDE returned zero ASG residue assignments."

    metrics = {
        f"pct_{ss_class}": 100.0 * counts[ss_class] / total_residues for ss_class in STRIDE_CLASSES
    }
    metrics["pct_helix_total"] = metrics["pct_H"] + metrics["pct_G"] + metrics["pct_I"]
    metrics["pct_extended_total"] = metrics["pct_E"] + metrics["pct_B"]
    metrics["stride_total_residues"] = total_residues

    warning = None
    if unknown_classes:
        unknown_repr = ",".join(f"{key}:{value}" for key, value in sorted(unknown_classes.items()))
        warning = f"Unknown STRIDE classes encountered ({unknown_repr})."
    return metrics, warning


def _collect_chain_coords(
    atoms: Iterable[Tuple[str, str, str, float, float, float]],
    chain_id: str,
    atom_filter: Optional[set[str]],
) -> np.ndarray:
    coords = []
    for atom_name, current_chain, _residue_idx, x_coord, y_coord, z_coord in atoms:
        if current_chain != chain_id:
            continue
        if atom_filter is not None and atom_name not in atom_filter:
            continue
        coords.append((x_coord, y_coord, z_coord))

    if not coords:
        return np.empty((0, 3), dtype=float)
    return np.asarray(coords, dtype=float)


def _make_base_row(pdb_path: Path) -> Dict[str, object]:
    row: Dict[str, object] = {column: None for column in OUTPUT_COLUMNS}
    row["pdb_file"] = pdb_path.name
    row["pdb_path"] = str(pdb_path.resolve())
    row["status"] = "OK"
    row["error_message"] = ""
    return row


def analyze_pdb(
    pdb_path_str: str,
    binder_chain: str,
    rg_chain: str,
    atom_selection: str,
    interface_point_xyz: Sequence[float],
    stride_executable: str,
    stride_timeout: int,
    fail_on_missing_stride: bool,
) -> Dict[str, object]:
    """Analyze one PDB and return a row dict for output."""
    pdb_path = Path(pdb_path_str)
    row = _make_base_row(pdb_path)
    warnings: List[str] = []

    try:
        atoms = parse_pdb_atoms(pdb_path)

        binder_ca = _collect_chain_coords(atoms, binder_chain, {"CA"})
        n_ca_binder, nterm_sc, cterm_sc, end_warning = compute_end_scores(
            binder_ca, np.asarray(interface_point_xyz, dtype=float)
        )
        row["n_ca_binder"] = n_ca_binder
        row["nterm_sc"] = nterm_sc
        row["cterm_sc"] = cterm_sc
        if end_warning:
            warnings.append(end_warning)

        atom_filter = ATOM_SELECTIONS[atom_selection]
        rg_coords = _collect_chain_coords(atoms, rg_chain, atom_filter)
        n_atoms_rg, rg, rg_warning = compute_rg(rg_coords)
        row["n_atoms_rg"] = n_atoms_rg
        row["rg"] = rg
        if rg_warning:
            warnings.append(rg_warning)

        try:
            stride_metrics, stride_warning = run_stride_and_parse(
                stride_executable, pdb_path, timeout_seconds=stride_timeout
            )
            row.update(stride_metrics)
            if stride_warning:
                warnings.append(stride_warning)
        except Exception as exc:
            if fail_on_missing_stride:
                raise StrideFailureError(f"{pdb_path}: {exc}") from exc
            row.update(empty_stride_metrics())
            warnings.append(f"STRIDE failed: {exc}")

        if warnings:
            row["status"] = "WARNING"
            row["error_message"] = "; ".join(warnings)
        else:
            row["status"] = "OK"
            row["error_message"] = ""
        return row

    except StrideFailureError:
        raise
    except Exception as exc:
        row.update(empty_stride_metrics())
        row["status"] = "ERROR"
        row["error_message"] = str(exc)
        return row


def discover_pdb_files(pdb_dir: Path, glob_pattern: str, recursive: bool) -> List[Path]:
    if recursive:
        pdb_files = [path for path in pdb_dir.rglob(glob_pattern) if path.is_file()]
    else:
        pdb_files = [path for path in pdb_dir.glob(glob_pattern) if path.is_file()]
    return sorted(path.resolve() for path in pdb_files)


def _serialize_cell(value: object) -> str:
    if value is None:
        return "NA"
    if isinstance(value, float):
        if math.isnan(value) or math.isinf(value):
            return "NA"
        return f"{value:.6f}"
    return str(value)


def write_output(rows: List[Dict[str, object]], output_path: Path, delimiter_name: str) -> None:
    delimiter = "," if delimiter_name == "csv" else "\t"
    rows.sort(key=lambda row: (str(row.get("pdb_file", "")), str(row.get("pdb_path", ""))))

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(OUTPUT_COLUMNS), delimiter=delimiter)
        writer.writeheader()
        for row in rows:
            writer.writerow({column: _serialize_cell(row.get(column)) for column in OUTPUT_COLUMNS})


def resolve_stride_executable(stride_arg: str) -> str:
    stride_path = Path(stride_arg)
    if stride_path.is_file():
        return str(stride_path.resolve())

    found = shutil.which(stride_arg)
    if found:
        return found

    raise FileNotFoundError(
        f"Could not resolve STRIDE executable from '{stride_arg}'. "
        "Pass an absolute path or a command available in PATH."
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Analyze a directory of PDB backbones and output per-PDB metrics "
            "(end scores, Rg, STRIDE secondary structure percentages)."
        ),
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=(
            "Example usage:\n"
            "  python analyze_backbones.py \\\n"
            "    --pdb-dir ./pdbs \\\n"
            "    --binder-chain A \\\n"
            "    --interface-point 0 0 -16.17 \\\n"
            "    --stride /path/to/stride \\\n"
            "    --atoms alpha-carbons \\\n"
            "    --output metrics.csv \\\n"
            "    --threads 8\n"
        ),
    )

    parser.add_argument("--pdb-dir", required=True, type=Path, help="Directory containing PDB files.")
    parser.add_argument(
        "--binder-chain",
        required=True,
        type=str,
        help="Chain ID used for NtermSc/CtermSc computation.",
    )
    parser.add_argument(
        "--interface-point",
        required=True,
        nargs=3,
        type=float,
        metavar=("X", "Y", "Z"),
        help="Interface point coordinates (X Y Z).",
    )
    parser.add_argument(
        "--stride",
        required=True,
        type=str,
        help="Path/name for STRIDE executable (used as: stride -f <pdb>).",
    )

    parser.add_argument(
        "--rg-chain",
        default=None,
        type=str,
        help="Chain ID for Rg computation (default: same as --binder-chain).",
    )
    parser.add_argument(
        "--atoms",
        default="alpha-carbons",
        choices=tuple(ATOM_SELECTIONS.keys()),
        help=(
            "Atom selection for Rg: "
            "alpha-carbons=CA; backbone=N,CA,C,O; "
            "backbone-no-carbonyl=N,CA,C; all=all ATOM records."
        ),
    )
    parser.add_argument("--output", default="metrics.csv", type=Path, help="Output CSV/TSV file path.")
    parser.add_argument(
        "--delimiter",
        default="csv",
        choices=("csv", "tsv"),
        help="Output delimiter format.",
    )
    parser.add_argument(
        "--recursive",
        action="store_true",
        help="Recursively search for PDB files in subdirectories.",
    )
    parser.add_argument(
        "--glob",
        default="*.pdb",
        type=str,
        help="Glob pattern for PDB discovery (default: *.pdb).",
    )
    parser.add_argument(
        "--fail-on-missing-stride",
        action="store_true",
        help="Fail immediately if STRIDE fails for any PDB.",
    )
    parser.add_argument(
        "--threads",
        default=1,
        type=int,
        help="Number of worker processes (ProcessPoolExecutor).",
    )
    parser.add_argument(
        "--stride-timeout",
        default=30,
        type=int,
        help="Timeout in seconds for each STRIDE call (default: 30).",
    )
    return parser


def _print_progress(done: int, total: int) -> None:
    if done % 25 == 0 or done == total:
        print(f"[{done}/{total}] analyzed", file=sys.stderr)


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    if args.threads < 1:
        parser.error("--threads must be >= 1")
    if args.stride_timeout < 1:
        parser.error("--stride-timeout must be >= 1")
    if not args.pdb_dir.exists() or not args.pdb_dir.is_dir():
        parser.error(f"--pdb-dir does not exist or is not a directory: {args.pdb_dir}")

    try:
        stride_executable = resolve_stride_executable(args.stride)
    except FileNotFoundError as exc:
        parser.error(str(exc))
        return 2

    rg_chain = args.rg_chain if args.rg_chain is not None else args.binder_chain
    pdb_files = discover_pdb_files(args.pdb_dir, args.glob, args.recursive)
    if not pdb_files:
        parser.error(
            f"No PDB files found under {args.pdb_dir} using pattern '{args.glob}' "
            f"(recursive={args.recursive})."
        )

    total = len(pdb_files)
    print(f"Discovered {total} PDB files.", file=sys.stderr)

    task_args = [
        (
            str(pdb_path),
            args.binder_chain,
            rg_chain,
            args.atoms,
            tuple(args.interface_point),
            stride_executable,
            args.stride_timeout,
            args.fail_on_missing_stride,
        )
        for pdb_path in pdb_files
    ]

    rows: List[Dict[str, object]] = []
    try:
        if args.threads == 1:
            for index, task in enumerate(task_args, start=1):
                rows.append(analyze_pdb(*task))
                _print_progress(index, total)
        else:
            try:
                with ProcessPoolExecutor(max_workers=args.threads) as executor:
                    futures = [executor.submit(analyze_pdb, *task) for task in task_args]
                    for index, future in enumerate(as_completed(futures), start=1):
                        rows.append(future.result())
                        _print_progress(index, total)
            except PermissionError as exc:
                print(
                    "WARNING: ProcessPoolExecutor is unavailable in this environment "
                    f"({exc}). Falling back to serial execution.",
                    file=sys.stderr,
                )
                rows = []
                for index, task in enumerate(task_args, start=1):
                    rows.append(analyze_pdb(*task))
                    _print_progress(index, total)
    except StrideFailureError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    write_output(rows, args.output.resolve(), args.delimiter)
    print(f"Wrote {len(rows)} rows to {args.output.resolve()}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
