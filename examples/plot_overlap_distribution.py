#!/usr/bin/env python3
"""Plot target-link and correction-frame overlap distributions."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "csv",
        type=Path,
        nargs="?",
        default=Path("results/paper_q20_converged/overlap_distribution.csv"),
    )
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()
    values = np.genfromtxt(args.csv, delimiter=",", names=True)

    figure, axes = plt.subplots(1, 3, figsize=(10.8, 3.25), constrained_layout=True)

    for field, label in (
        ("raw_link1_abs", r"$|M_{00}^{(1)}|$"),
        ("raw_link2_abs", r"$|M_{00}^{(2)}|$"),
    ):
        axes[0].hist(values[field], bins=24, density=True, histtype="step", lw=1.8, label=label)
    axes[0].set_xlabel("raw target-link modulus")

    for field, label in (
        ("covariant_link1_abs", r"$1-|U_{00}^{(1)}|^2$"),
        ("covariant_link2_abs", r"$1-|U_{00}^{(2)}|^2$"),
    ):
        distance = np.maximum(0.0, 1.0 - values[field] ** 2)
        axes[1].hist(distance, bins=24, density=True, histtype="step", lw=1.8, label=label)
    axes[1].set_xlabel("covariant projector distance")
    axes[1].ticklabel_format(axis="x", style="sci", scilimits=(-2, 2))

    for field, label in (
        ("frame_singular1", r"$s_{\min}^{(1)}$"),
        ("frame_singular2", r"$s_{\min}^{(2)}$"),
        ("frame_singular_difference", r"$s_{\min}^{(2-1)}$"),
    ):
        axes[2].hist(values[field], bins=24, density=True, histtype="step", lw=1.8, label=label)
    axes[2].set_xlabel("minimum frame singular value")

    for axis in axes:
        axis.set_ylabel("density")
        axis.legend(frameon=False, fontsize=8)

    output = args.output or args.csv.with_name("overlap_distribution.png")
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=240)
    print(output.resolve())


if __name__ == "__main__":
    main()
