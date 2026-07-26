#!/usr/bin/env python3
"""Plot Berry curvature, Quantum Metric, and trace excess from LL output."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def main(input_path: Path, output_stem: Path) -> None:
    data = np.load(input_path)
    nk1, nk2 = int(data["nk1"]), int(data["nk2"])
    mbz_area = float(data["plaquette_area"]) * nk1 * nk2
    magnetic_length_squared = 2.0 * np.pi / mbz_area
    chern = float(data["chern"])

    curvature = np.asarray(data["berry_curvature"], dtype=float).T
    metric_trace = np.asarray(data["metric_trace"], dtype=float).T
    fields = (
        (
            chern * curvature / magnetic_length_squared,
            r"$C\Omega/l_B^2$",
            "coolwarm",
            (0.98, 1.02),
        ),
        (
            metric_trace / magnetic_length_squared,
            r"$\mathrm{Tr}\,g/l_B^2$",
            "viridis",
            (0.98, 1.06),
        ),
        (
            (metric_trace - np.abs(curvature)) / magnetic_length_squared,
            r"$(\mathrm{Tr}\,g-|\Omega|)/l_B^2$",
            "magma",
            (0.0, 0.04),
        ),
    )

    plt.rcParams.update(
        {
            "font.size": 9,
            "axes.labelsize": 9,
            "axes.titlesize": 9,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    figure, axes = plt.subplots(
        1,
        3,
        figsize=(8.4, 2.8),
        layout="constrained",
    )
    for index, (axis, (field, title, colormap, limits)) in enumerate(
        zip(axes, fields)
    ):
        image = axis.pcolormesh(
            np.linspace(0.0, 1.0, nk1 + 1),
            np.linspace(0.0, 1.0, nk2 + 1),
            field,
            shading="flat",
            cmap=colormap,
            vmin=limits[0],
            vmax=limits[1],
            rasterized=True,
        )
        axis.set_title(f"{chr(97 + index)}   {title}", loc="left")
        axis.set_xlabel(r"$k_x/Q_x$")
        axis.set_xticks([0.0, 0.5, 1.0])
        axis.set_yticks([0.0, 0.5, 1.0])
        if index == 0:
            axis.set_ylabel(r"$q_{\mathrm{LL}}k_y/Q_y$")
        else:
            axis.set_yticklabels([])
        axis.tick_params(direction="in", top=True, right=True)
        figure.colorbar(image, ax=axis, pad=0.025, fraction=0.05)

    figure.suptitle(
        rf"LL reference: $N_{{\rm LL}}={int(data['nll'])}$, "
        rf"{nk1}$\times${nk2}, $C={chern:.0f}$, "
        rf"$\eta={float(data['eta']):.4g}$",
        fontsize=9,
    )
    output_stem.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_stem.with_suffix(".pdf"), bbox_inches="tight")
    figure.savefig(output_stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(figure)
    print(output_stem.with_suffix(".pdf").resolve())
    print(output_stem.with_suffix(".png").resolve())


if __name__ == "__main__":
    repository_root = Path(__file__).resolve().parents[2]
    default_input = (
        repository_root
        / "reference"
        / "tmbg_pq1_20_ll_nll400_40x16"
        / "quantum_geometry.npz"
    )
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=default_input)
    parser.add_argument("--output-stem", type=Path)
    arguments = parser.parse_args()
    main(
        arguments.input,
        arguments.output_stem or arguments.input.with_suffix(""),
    )
