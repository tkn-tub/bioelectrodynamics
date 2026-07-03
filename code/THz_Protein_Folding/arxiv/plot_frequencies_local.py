#!/usr/bin/env python3
"""
Plot ProDy ANM frequency results on your local computer.

Input:
  results/frequencies_long.csv produced by compute_prody_frequencies.py

Outputs:
  local_plots/vdos_overlay.png
  local_plots/histograms/<label>_histogram.png
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def read_frequencies(path: Path) -> dict[str, np.ndarray]:
    grouped: dict[str, list[float]] = defaultdict(list)
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"label", "frequency_thz"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing required columns in {path}: {sorted(missing)}")

        for row in reader:
            grouped[row["label"]].append(float(row["frequency_thz"]))

    return {label: np.asarray(values, dtype=float) for label, values in grouped.items()}


def plot_overlay(
    grouped: dict[str, np.ndarray],
    outpath: Path,
    xlim: tuple[float, float] = (0.0, 2.0),
    bins: int = 120,
) -> None:
    outpath.parent.mkdir(parents=True, exist_ok=True)

    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(111)

    edges = np.linspace(xlim[0], xlim[1], bins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])

    for label, freqs in grouped.items():
        in_range = freqs[(freqs >= xlim[0]) & (freqs <= xlim[1])]
        if len(in_range) == 0:
            continue
        density, _ = np.histogram(in_range, bins=edges, density=True)
        ax.plot(centers, density, label=f"{label} ({len(freqs)} modes)")

    ax.axvspan(0.1, 3.0, alpha=0.1, label="Standard THz range")
    ax.axvline(0.1, linestyle="--", linewidth=1)
    ax.set_title("Vibrational Density of States (VDOS)")
    ax.set_xlabel("Frequency (THz)")
    ax.set_ylabel("Density of modes")
    ax.set_xlim(*xlim)
    ax.grid(axis="y", alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(outpath, dpi=300)
    plt.close(fig)


def plot_histograms(grouped: dict[str, np.ndarray], outdir: Path, bins: int = 40) -> None:
    outdir.mkdir(parents=True, exist_ok=True)

    for label, freqs in grouped.items():
        fig = plt.figure(figsize=(10, 6))
        ax = fig.add_subplot(111)
        ax.hist(freqs, bins=bins, edgecolor="black", alpha=0.7)
        ax.axvline(0.1, linestyle="--", label="100 GHz = 0.1 THz")
        ax.axvline(3.0, linestyle="--", label="3.0 THz")
        ax.set_title(f"Full Vibrational Frequency Distribution for {label}")
        ax.set_xlabel("Frequency (THz)")
        ax.set_ylabel("Number of modes")
        ax.grid(axis="y", linestyle=":", alpha=0.6)
        ax.legend()
        fig.tight_layout()
        fig.savefig(outdir / f"{label}_histogram.png", dpi=300)
        plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description="Plot frequencies calculated on HPC.")
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("results") / "frequencies_long.csv",
        help="Path to frequencies_long.csv from the HPC run.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("local_plots"),
        help="Directory for PNG plots.",
    )
    parser.add_argument(
        "--xmax",
        type=float,
        default=2.0,
        help="Maximum x value for overlay VDOS plot. Default: 2.0 THz",
    )
    args = parser.parse_args()

    grouped = read_frequencies(args.input)
    if not grouped:
        raise RuntimeError(f"No frequencies found in {args.input}")

    plot_overlay(grouped, args.outdir / "vdos_overlay.png", xlim=(0.0, args.xmax))
    plot_histograms(grouped, args.outdir / "histograms")

    print(f"Wrote {args.outdir / 'vdos_overlay.png'}")
    print(f"Wrote individual histograms in {args.outdir / 'histograms'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
