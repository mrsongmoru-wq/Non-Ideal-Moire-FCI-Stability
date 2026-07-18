#!/usr/bin/env python3
"""Plot the Berry curvature and quantum metric from the paper example."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def grid(values: np.ndarray, field: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    k1 = np.unique(values["k1"])
    k2 = np.unique(values["k2"])
    z = values[field].reshape(len(k2), len(k1))
    return k1, k2, z


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "csv",
        type=Path,
        nargs="?",
        default=Path("results/paper_q20_converged/ideal_component_geometry.csv"),
    )
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()

    values = np.genfromtxt(args.csv, delimiter=",", names=True)
    fields = (
        ("berry_curvature", r"$\Omega(\mathbf{k})$"),
        ("intrinsic_trace", r"$\mathrm{Tr}\,g_{\rm intrinsic}(\mathbf{k})$"),
        ("trace_excess", r"$\mathrm{Tr}\,g_{\rm residual}(\mathbf{k})$"),
    )
    figure, axes = plt.subplots(1, 3, figsize=(11.2, 3.35), constrained_layout=True)
    for axis, (field, title) in zip(axes, fields):
        k1, k2, z = grid(values, field)
        image = axis.pcolormesh(k1, k2, z, shading="nearest", cmap="viridis")
        axis.set_xlabel(r"$k_1$")
        axis.set_ylabel(r"$k_2$")
        axis.set_title(title)
        figure.colorbar(image, ax=axis)

    output = args.output or args.csv.with_name("ideal_component_geometry.png")
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=240)
    print(output.resolve())


if __name__ == "__main__":
    main()
