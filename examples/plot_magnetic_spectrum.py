#!/usr/bin/env python3
"""Plot the converged magnetic hybrid-Wannier subband spectrum."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def main(
    spectrum_path: Path,
    output_stem: Path,
    energy_min: float | None,
    energy_max: float | None,
) -> None:
    spectrum = np.genfromtxt(spectrum_path, delimiter=",", names=True)
    required = {"k1", "k2", "magnetic_band", "energy_meV"}
    available = set(spectrum.dtype.names or ())
    if not required.issubset(available):
        raise RuntimeError(
            f"Magnetic-spectrum table is missing columns: {required - available}"
        )

    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 8.2,
            "axes.labelsize": 8.2,
            "xtick.labelsize": 7.6,
            "ytick.labelsize": 7.6,
            "axes.linewidth": 0.75,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    figure, axis = plt.subplots(figsize=(3.45, 2.45))
    colors = plt.get_cmap("viridis")(
        np.linspace(0.18, 0.82, len(np.unique(spectrum["k2"])))
    )

    for color, k2 in zip(colors, np.unique(spectrum["k2"]), strict=True):
        k2_rows = spectrum[np.isclose(spectrum["k2"], k2)]
        first_curve = True
        for band in np.unique(k2_rows["magnetic_band"]):
            band_rows = k2_rows[k2_rows["magnetic_band"] == band]
            order = np.argsort(band_rows["k1"])
            axis.plot(
                band_rows["k1"][order],
                band_rows["energy_meV"][order],
                color=color,
                linewidth=0.42,
                alpha=0.86,
                label=rf"$k_2={k2:.3f}$" if first_curve else None,
            )
            first_curve = False

    axis.set_xlabel(r"magnetic momentum $k_1$")
    axis.set_ylabel("Energy (meV)")
    axis.set_xlim(0.0, 1.0)
    if energy_min is not None or energy_max is not None:
        current_min, current_max = axis.get_ylim()
        axis.set_ylim(
            current_min if energy_min is None else energy_min,
            current_max if energy_max is None else energy_max,
        )
    axis.tick_params(direction="in", length=2.5, width=0.7)
    axis.legend(frameon=False, fontsize=7.0, loc="best")
    figure.tight_layout(pad=0.45)

    output_stem.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_stem.with_suffix(".pdf"), dpi=500)
    figure.savefig(output_stem.with_suffix(".png"), dpi=500)
    plt.close(figure)
    print(output_stem.with_suffix(".pdf").resolve())


if __name__ == "__main__":
    repository_root = Path(__file__).resolve().parents[1]
    default_directory = (
        repository_root / "results" / "tmbg_pq1_20_lg11_smax5"
    )
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--spectrum",
        type=Path,
        default=default_directory / "magnetic_spectrum.csv",
    )
    parser.add_argument("--output-stem", type=Path, default=None)
    parser.add_argument("--energy-min", type=float, default=15.0)
    parser.add_argument("--energy-max", type=float, default=40.0)
    parser.add_argument(
        "--full-range",
        action="store_true",
        help="Plot the full exported energy range instead of the target-band window.",
    )
    arguments = parser.parse_args()
    output = arguments.output_stem or arguments.spectrum.with_suffix("")
    energy_min = None if arguments.full_range else arguments.energy_min
    energy_max = None if arguments.full_range else arguments.energy_max
    main(arguments.spectrum, output, energy_min, energy_max)
