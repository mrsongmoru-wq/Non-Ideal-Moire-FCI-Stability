#!/usr/bin/env python3
"""Plot the magnetic spectrum exported by the LL expansion."""

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
    required = {"k1", "k2", "central_band_1based", "energy_meV"}
    available = set(spectrum.dtype.names or ())
    if not required.issubset(available):
        raise RuntimeError(
            f"LL spectrum table is missing columns: {required - available}"
        )

    figure, axis = plt.subplots(figsize=(3.6, 2.6))
    k2_values = np.unique(spectrum["k2"])
    colors = plt.get_cmap("viridis")(np.linspace(0.15, 0.85, len(k2_values)))
    for color, k2 in zip(colors, k2_values):
        rows_at_k2 = spectrum[np.isclose(spectrum["k2"], k2)]
        for band in np.unique(rows_at_k2["central_band_1based"]):
            band_rows = rows_at_k2[
                rows_at_k2["central_band_1based"] == band
            ]
            order = np.argsort(band_rows["k1"])
            axis.plot(
                band_rows["k1"][order],
                band_rows["energy_meV"][order],
                color=color,
                linewidth=0.45,
                alpha=0.85,
            )

    axis.set_xlabel(r"magnetic momentum $k_1$")
    axis.set_ylabel("Energy (meV)")
    axis.set_xlim(0.0, 1.0)
    if energy_min is not None or energy_max is not None:
        current_min, current_max = axis.get_ylim()
        axis.set_ylim(
            current_min if energy_min is None else energy_min,
            current_max if energy_max is None else energy_max,
        )
    axis.tick_params(direction="in", top=True, right=True)
    figure.tight_layout(pad=0.5)
    output_stem.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_stem.with_suffix(".pdf"), bbox_inches="tight")
    figure.savefig(output_stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(figure)
    print(output_stem.with_suffix(".pdf").resolve())
    print(output_stem.with_suffix(".png").resolve())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--spectrum",
        type=Path,
        default=Path("results/ll_reference/magnetic_spectrum.csv"),
    )
    parser.add_argument("--output-stem", type=Path)
    parser.add_argument("--energy-min", type=float)
    parser.add_argument("--energy-max", type=float)
    arguments = parser.parse_args()
    main(
        arguments.spectrum,
        arguments.output_stem or arguments.spectrum.with_suffix(""),
        arguments.energy_min,
        arguments.energy_max,
    )
