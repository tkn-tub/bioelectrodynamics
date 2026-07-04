#!/usr/bin/env python3
"""
Compute ProDy ANM vibrational frequencies without plotting.

This script is intended for remote HPC use. It downloads or reads protein
structure files, computes ANM eigenvalues with ProDy, converts them to
frequencies using the scale factor used in the original notebook:

    frequency_THz = 0.225 * sqrt(eigenvalue)

Outputs:
  - results/frequencies_long.csv
  - results/summary.csv
  - results/frequencies.npz
  - results/run_manifest.json

The default structure list is taken from the uploaded Colab notebook. Edit
RCSB_STRUCTURES and ALPHAFOLD_STRUCTURES below if you want to change it.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
import traceback
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
from prody import ANM, parsePDB


# RCSB structures used in the notebook.
# The first four were used in the combined VDOS plot.
# 1TBU and 2THI also appear later in the notebook for individual calculations/animations.
RCSB_STRUCTURES = ["2AXQ", "1YAA", "3FKY", "3IG5"]

# AlphaFold structures downloaded from the GitHub repository in the notebook.
ALPHAFOLD_STRUCTURES = {
    "YPR021C": "AF-Q12482-F1-model_v6.pdb",
    "YLR174W": "AF-P41939-F1-model_v6.pdb",
    "YGR019W": "AF-P17649-F1-model_v6.pdb",
    "YAL062W": "AF-P39708-F1-model_v6.pdb",
    "YMR113W": "AF-Q12676-F1-model_v6.pdb",
    "YLR017W": "AF-Q07938-F1-model_v6.pdb",
}

ALPHAFOLD_BASE_URL = (
    "https://raw.githubusercontent.com/tkn-tub/bioelectrodynamics/main/code/database/"
)
RCSB_BASE_URL = "https://files.rcsb.org/download/"


@dataclass(frozen=True)
class StructureJob:
    label: str
    filename: str
    url: str
    source_type: str


def parse_n_modes(value: str) -> int | None:
    value_clean = value.strip().lower()
    if value_clean in {"all", "none", "full"}:
        return None
    try:
        n_modes = int(value_clean)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            "--n-modes must be an integer or 'all'"
        ) from exc
    if n_modes <= 0:
        raise argparse.ArgumentTypeError("--n-modes must be positive or 'all'")
    return n_modes


def build_jobs(include_rcsb: bool = True, include_alphafold: bool = True) -> list[StructureJob]:
    jobs: list[StructureJob] = []

    if include_rcsb:
        for pdb_id in RCSB_STRUCTURES:
            pdb_id = pdb_id.upper()
            jobs.append(
                StructureJob(
                    label=pdb_id,
                    filename=f"{pdb_id}.pdb",
                    url=f"{RCSB_BASE_URL}{pdb_id}.pdb",
                    source_type="RCSB",
                )
            )

    if include_alphafold:
        for label, filename in ALPHAFOLD_STRUCTURES.items():
            jobs.append(
                StructureJob(
                    label=label,
                    filename=filename,
                    url=f"{ALPHAFOLD_BASE_URL}{filename}",
                    source_type="AlphaFold/GitHub",
                )
            )

    return jobs


def download_file(url: str, destination: Path, overwrite: bool = False) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and not overwrite:
        print(f"[download] exists: {destination}")
        return

    print(f"[download] {url} -> {destination}")
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "prody-hpc-frequency-script/1.0"},
    )
    with urllib.request.urlopen(request, timeout=120) as response:
        data = response.read()
    destination.write_bytes(data)


def download_structures(
    jobs: Iterable[StructureJob],
    structures_dir: Path,
    overwrite: bool = False,
) -> None:
    structures_dir.mkdir(parents=True, exist_ok=True)
    for job in jobs:
        download_file(job.url, structures_dir / job.filename, overwrite=overwrite)


def compute_one_structure(
    job: StructureJob,
    structure_path: Path,
    n_modes: int | None,
    frequency_scale: float,
) -> tuple[np.ndarray, np.ndarray, dict]:
    start = time.time()
    print(f"\n[compute] Processing {job.label} from {structure_path}")

    if not structure_path.exists():
        raise FileNotFoundError(
            f"Missing structure file: {structure_path}. "
            "Run this script with --download-only first, or run without --no-download."
        )

    structure = parsePDB(str(structure_path))
    if structure is None:
        raise RuntimeError(f"ProDy could not parse {structure_path}")

    calphas = structure.select("protein and name CA")
    if calphas is None or len(calphas) == 0:
        raise RuntimeError(f"No protein C-alpha atoms found for {job.label}")

    anm = ANM(job.label)
    anm.buildHessian(calphas)
    anm.calcModes(n_modes=n_modes)

    eigvals = np.asarray(anm.getEigvals(), dtype=float)
    negative_eigvals = int(np.sum(eigvals < 0.0))

    # Very small negative eigenvalues can occur from numerical roundoff.
    # Clipping avoids NaN values from sqrt while preserving the count in summary.csv.
    frequencies = frequency_scale * np.sqrt(np.clip(eigvals, 0.0, None))

    elapsed = time.time() - start
    summary = {
        "label": job.label,
        "source_type": job.source_type,
        "source_file": str(structure_path),
        "n_calpha": int(len(calphas)),
        "n_modes": int(len(frequencies)),
        "min_frequency_thz": float(np.min(frequencies)) if len(frequencies) else "",
        "max_frequency_thz": float(np.max(frequencies)) if len(frequencies) else "",
        "mean_frequency_thz": float(np.mean(frequencies)) if len(frequencies) else "",
        "modes_microwave_lt_0_1_thz": int(np.sum(frequencies < 0.1)),
        "modes_thz_0_1_to_3_0_thz": int(np.sum((frequencies >= 0.1) & (frequencies <= 3.0))),
        "modes_beyond_3_0_thz": int(np.sum(frequencies > 3.0)),
        "negative_eigenvalues_clipped": negative_eigvals,
        "elapsed_seconds": round(elapsed, 3),
        "status": "ok",
        "error": "",
    }

    print(
        f"[compute] {job.label}: {summary['n_calpha']} C-alpha atoms, "
        f"{summary['n_modes']} modes, "
        f"{summary['min_frequency_thz']:.6g}-{summary['max_frequency_thz']:.6g} THz"
    )

    return eigvals, frequencies, summary


def write_outputs(
    output_dir: Path,
    summaries: list[dict],
    eigenvalues_by_label: dict[str, np.ndarray],
    frequencies_by_label: dict[str, np.ndarray],
    manifest: dict,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    summary_path = output_dir / "summary.csv"
    long_path = output_dir / "frequencies_long.csv"
    npz_path = output_dir / "frequencies.npz"
    manifest_path = output_dir / "run_manifest.json"

    summary_fields = [
        "label",
        "source_type",
        "source_file",
        "n_calpha",
        "n_modes",
        "min_frequency_thz",
        "max_frequency_thz",
        "mean_frequency_thz",
        "modes_microwave_lt_0_1_thz",
        "modes_thz_0_1_to_3_0_thz",
        "modes_beyond_3_0_thz",
        "negative_eigenvalues_clipped",
        "elapsed_seconds",
        "status",
        "error",
    ]

    with summary_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_fields)
        writer.writeheader()
        for row in summaries:
            writer.writerow({field: row.get(field, "") for field in summary_fields})

    with long_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["label", "mode_index", "eigenvalue", "frequency_thz"])
        for label, freqs in frequencies_by_label.items():
            eigvals = eigenvalues_by_label[label]
            for mode_index, (eigval, freq) in enumerate(zip(eigvals, freqs), start=1):
                writer.writerow([label, mode_index, f"{eigval:.17g}", f"{freq:.17g}"])

    # Store each protein's frequency array as a separate array in one compact binary file.
    np.savez_compressed(npz_path, **frequencies_by_label)

    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    print("\n[done] Wrote outputs:")
    print(f"  {summary_path}")
    print(f"  {long_path}")
    print(f"  {npz_path}")
    print(f"  {manifest_path}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Compute ProDy ANM frequencies on an HPC node without plotting."
    )
    parser.add_argument(
        "--structures-dir",
        type=Path,
        default=Path("structures"),
        help="Directory containing/downloading PDB files. Default: structures",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("../database/results_motif_1"),
        help="Output directory for CSV/NPZ files. Default: ../database/results_1",
    )
    parser.add_argument(
        "--n-modes",
        type=parse_n_modes,
        default=None,
        help="Number of ANM modes to calculate, or 'all'. Default: all",
    )
    parser.add_argument(
        "--frequency-scale",
        type=float,
        default=0.225,
        help="Scale factor for frequency_THz = scale * sqrt(eigenvalue). Default: 0.225",
    )
    parser.add_argument(
        "--no-download",
        action="store_true",
        help="Do not download files; require all structures to exist in --structures-dir.",
    )
    parser.add_argument(
        "--download-only",
        action="store_true",
        help="Download structure files and exit without calculations.",
    )
    parser.add_argument(
        "--overwrite-downloads",
        action="store_true",
        help="Re-download structure files even if they already exist.",
    )
    parser.add_argument(
        "--only-rcsb",
        action="store_true",
        help="Run only RCSB/PDB structures from RCSB_STRUCTURES.",
    )
    parser.add_argument(
        "--only-alphafold",
        action="store_true",
        help="Run only AlphaFold/GitHub structures from ALPHAFOLD_STRUCTURES.",
    )
    parser.add_argument(
        "--labels",
        nargs="+",
        default=None,
        help=(
            "Optional list of structure labels to run, e.g. --labels 3PYM "
            "or --labels 3FKY 3PYM YDL215C. Matching is case-insensitive."
        ),
    )

    args = parser.parse_args(argv)

    if args.only_rcsb and args.only_alphafold:
        parser.error("Use only one of --only-rcsb or --only-alphafold.")

    include_rcsb = not args.only_alphafold
    include_alphafold = not args.only_rcsb

    jobs = build_jobs(include_rcsb=include_rcsb, include_alphafold=include_alphafold)

    if args.labels:
        requested_labels = {label.upper() for label in args.labels}
        jobs = [job for job in jobs if job.label.upper() in requested_labels]
        found_labels = {job.label.upper() for job in jobs}
        missing_labels = sorted(requested_labels - found_labels)
        if missing_labels:
            parser.error(f"Unknown label(s): {', '.join(missing_labels)}")
        if not jobs:
            parser.error("No jobs selected after applying --labels.")

    manifest = {
        "script": Path(__file__).name,
        "timestamp_unix": time.time(),
        "n_modes_requested": "all" if args.n_modes is None else args.n_modes,
        "frequency_scale": args.frequency_scale,
        "structures_dir": str(args.structures_dir),
        "outdir": str(args.outdir),
        "jobs": [job.__dict__ for job in jobs],
        "python_version": sys.version,
        "numpy_version": np.__version__,
    }

    if not args.no_download:
        download_structures(
            jobs,
            structures_dir=args.structures_dir,
            overwrite=args.overwrite_downloads,
        )

    if args.download_only:
        print("\n[done] Download-only mode completed.")
        return 0

    summaries: list[dict] = []
    eigenvalues_by_label: dict[str, np.ndarray] = {}
    frequencies_by_label: dict[str, np.ndarray] = {}

    for job in jobs:
        structure_path = args.structures_dir / job.filename
        try:
            eigvals, freqs, summary = compute_one_structure(
                job=job,
                structure_path=structure_path,
                n_modes=args.n_modes,
                frequency_scale=args.frequency_scale,
            )
            eigenvalues_by_label[job.label] = eigvals
            frequencies_by_label[job.label] = freqs
            summaries.append(summary)
        except Exception as exc:
            print(f"\n[error] {job.label}: {exc}", file=sys.stderr)
            traceback.print_exc()
            summaries.append(
                {
                    "label": job.label,
                    "source_type": job.source_type,
                    "source_file": str(structure_path),
                    "status": "failed",
                    "error": repr(exc),
                }
            )

    write_outputs(
        output_dir=args.outdir,
        summaries=summaries,
        eigenvalues_by_label=eigenvalues_by_label,
        frequencies_by_label=frequencies_by_label,
        manifest=manifest,
    )

    n_ok = sum(1 for row in summaries if row.get("status") == "ok")
    n_failed = sum(1 for row in summaries if row.get("status") == "failed")
    print(f"\n[summary] Successful structures: {n_ok}; failed structures: {n_failed}")

    return 0 if n_ok > 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
