#!/usr/bin/env python3
"""Reproduce the path-and-hexagon ideal-component distribution of Fig. 4."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
from matplotlib.path import Path as MplPath
from matplotlib.ticker import FormatStrFormatter


def read_metadata(path: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            result[key.strip()] = value.strip()
    return result


def periodic_bilinear(grid: np.ndarray, k1, k2):
    n1, n2 = grid.shape
    x = (np.asarray(k1) % 1.0) * n1
    y = (np.asarray(k2) % 1.0) * n2
    i0 = np.floor(x).astype(int) % n1
    j0 = np.floor(y).astype(int) % n2
    i1, j1 = (i0 + 1) % n1, (j0 + 1) % n2
    tx, ty = x - np.floor(x), y - np.floor(y)
    return (
        (1 - tx) * (1 - ty) * grid[i0, j0]
        + tx * (1 - ty) * grid[i1, j0]
        + (1 - tx) * ty * grid[i0, j1]
        + tx * ty * grid[i1, j1]
    )


def colored_line(axis, x, y, values, cmap, norm):
    points = np.column_stack((x, y)).reshape(-1, 1, 2)
    segments = np.concatenate((points[:-1], points[1:]), axis=1)
    collection = LineCollection(
        segments,
        cmap=cmap,
        norm=norm,
        linewidth=2.0,
        capstyle="round",
        zorder=5,
    )
    collection.set_array(0.5 * (values[:-1] + values[1:]))
    axis.add_collection(collection)


def wigner_seitz_hexagon(g1: np.ndarray, g2: np.ndarray) -> np.ndarray:
    reciprocals = [
        m * g1 + n * g2
        for m in range(-2, 3)
        for n in range(-2, 3)
        if m != 0 or n != 0
    ]
    vertices: list[np.ndarray] = []
    for index, first in enumerate(reciprocals):
        for second in reciprocals[index + 1 :]:
            matrix = np.vstack((first, second))
            if abs(np.linalg.det(matrix)) < 1e-12:
                continue
            point = np.linalg.solve(
                matrix,
                np.array((np.dot(first, first), np.dot(second, second))) / 2,
            )
            if all(
                np.dot(point, vector) <= np.dot(vector, vector) / 2 + 1e-10
                for vector in reciprocals
            ) and not any(np.linalg.norm(point - old) < 1e-9 for old in vertices):
                vertices.append(point)
    polygon = np.asarray(vertices)
    return polygon[np.argsort(np.arctan2(polygon[:, 1], polygon[:, 0]))]


def fold_to_hexagon(reduced, gamma_reduced, basis, polygon):
    path = MplPath(polygon)
    best, best_norm = None, np.inf
    for m in range(-2, 3):
        for n in range(-2, 3):
            cartesian = basis @ (reduced - gamma_reduced - np.array((m, n)))
            distance = np.linalg.norm(cartesian)
            if path.contains_point(cartesian, radius=1e-10):
                return cartesian
            if distance < best_norm:
                best, best_norm = cartesian, distance
    return best


def build_hexagonal_map(grid, g1, g2, gamma_reduced, resolution=520):
    basis = np.column_stack((g1, g2))
    polygon = wigner_seitz_hexagon(g1, g2)
    x_axis = np.linspace(polygon[:, 0].min(), polygon[:, 0].max(), resolution)
    y_axis = np.linspace(polygon[:, 1].min(), polygon[:, 1].max(), resolution)
    x_grid, y_grid = np.meshgrid(x_axis, y_axis)
    points = np.column_stack((x_grid.ravel(), y_grid.ravel()))
    inside = MplPath(polygon).contains_points(points, radius=1e-11)
    reduced = gamma_reduced[:, None] + np.linalg.solve(basis, points.T)
    values = periodic_bilinear(grid, reduced[0], reduced[1]).reshape(x_grid.shape)
    values[~inside.reshape(x_grid.shape)] = np.nan
    return x_grid, y_grid, values, polygon, basis


def main(data_dir: Path, output_stem: Path) -> None:
    path_metadata = read_metadata(data_dir / "zero_field_path_metadata.txt")
    projection_metadata = read_metadata(data_dir / "metadata.txt")
    g1 = np.fromstring(path_metadata["g1"], sep=",")
    g2 = np.fromstring(path_metadata["g2"], sep=",")
    gamma_reduced = np.array((2 / 3, 1 / 3))

    projection = np.genfromtxt(
        data_dir / "upper_component_grid.csv", delimiter=",", names=True
    )
    if "target_weight_raw" not in (projection.dtype.names or ()):
        raise RuntimeError("Fig. 4 table does not contain the auditable raw weight")
    dimension = round(np.sqrt(projection.size))
    raw_weight = np.asarray(projection["target_weight"]).reshape(dimension, dimension).T
    relative_weight = raw_weight / np.max(raw_weight)
    vmin, vmax = float(np.min(relative_weight)), 1.0
    norm, cmap = Normalize(vmin=vmin, vmax=vmax), plt.get_cmap("Reds")

    path = np.genfromtxt(
        data_dir / "zero_field_path_kpoints.csv", delimiter=",", names=True
    )
    energies = np.loadtxt(
        data_dir / "zero_field_path_energies.csv", delimiter=",", skiprows=1
    )
    target_local = int(path_metadata["target_local_band"]) - 1
    node_indices = [
        int(value) for value in path_metadata["node_indices_zero_based"].split(",")
    ]
    path_weight = periodic_bilinear(
        relative_weight, path["k1_wang"], path["k2_wang"]
    )

    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 8.2,
            "axes.labelsize": 8.2,
            "xtick.labelsize": 7.6,
            "ytick.labelsize": 7.6,
            "axes.linewidth": 0.75,
            "mathtext.fontset": "stixsans",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    figure = plt.figure(figsize=(3.55, 1.9))
    layout = figure.add_gridspec(
        1,
        2,
        width_ratios=(1.12, 1.0),
        wspace=0.18,
        left=0.105,
        right=0.81,
        bottom=0.19,
        top=0.96,
    )
    band_axis = figure.add_subplot(layout[0, 0])
    hex_axis = figure.add_subplot(layout[0, 1])
    colorbar_axis = figure.add_axes((0.835, 0.255, 0.018, 0.57))

    coordinate = np.asarray(path["k_distance"])
    for band in range(energies.shape[1]):
        if band != target_local:
            band_axis.plot(
                coordinate,
                energies[:, band],
                color="#aeb6bf",
                linewidth=0.55,
                alpha=0.8,
            )
    colored_line(
        band_axis,
        coordinate,
        energies[:, target_local],
        path_weight,
        cmap,
        norm,
    )
    node_positions = coordinate[node_indices]
    band_axis.set_xticks(node_positions, [r"$K$", r"$\Gamma$", r"$M$", r"$K'$" ])
    for position in node_positions[1:-1]:
        band_axis.axvline(position, color="#d5d8dc", linewidth=0.55, zorder=0)
    band_axis.set_xlim(coordinate[0], coordinate[-1])
    band_axis.set_ylim(12, 30)
    band_axis.set_yticks((15, 20, 25, 30))
    band_axis.set_ylabel("Energy  meV", labelpad=2)
    band_axis.tick_params(direction="in", length=2.5, width=0.7)
    band_axis.text(0.025, 0.965, "a", transform=band_axis.transAxes, va="top", fontweight="bold")

    x_grid, y_grid, hex_weight, polygon, basis = build_hexagonal_map(
        relative_weight, g1, g2, gamma_reduced
    )
    image = hex_axis.pcolormesh(
        x_grid,
        y_grid,
        hex_weight,
        shading="auto",
        cmap=cmap,
        norm=norm,
        rasterized=True,
    )
    high_symmetry = {
        r"$\Gamma$": gamma_reduced,
        r"$K$": np.array((1 / 3, 2 / 3)),
        r"$M$": np.array((1 / 6, 1 / 3)),
        r"$K'$": np.array((0.0, 0.0)),
    }
    scale = np.linalg.norm(g1)
    for label, reduced in high_symmetry.items():
        point = fold_to_hexagon(reduced, gamma_reduced, basis, polygon)
        hex_axis.scatter(*point, s=5.5, color="#202124", linewidth=0, zorder=4)
        text_point = point + np.array((0.0, 0.045 * scale)) if np.linalg.norm(point) < 1e-12 else 0.84 * point
        hex_axis.text(*text_point, label, fontsize=7.6, ha="center", va="center")
    for spine in hex_axis.spines.values():
        spine.set_visible(False)
    hex_axis.set_aspect("equal")
    hex_axis.set_xlim(polygon[:, 0].min(), polygon[:, 0].max())
    hex_axis.set_ylim(polygon[:, 1].min(), polygon[:, 1].max())
    hex_axis.set_xticks([])
    hex_axis.set_yticks([])
    hex_axis.text(0.015, 0.965, "b", transform=hex_axis.transAxes, va="top", fontweight="bold")

    colorbar = figure.colorbar(image, cax=colorbar_axis)
    colorbar.set_ticks((vmin, 0.5 * (vmin + vmax), vmax))
    colorbar.ax.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    colorbar.set_label("Ideal component", labelpad=2.5, fontsize=7.8)
    colorbar.ax.tick_params(length=2.2, width=0.7, labelsize=7.2)

    output_stem.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_stem.with_suffix(".pdf"), dpi=500)
    figure.savefig(output_stem.with_suffix(".png"), dpi=500)
    plt.close(figure)

    source_weights = np.genfromtxt(
        data_dir / "upper_zero_field_band_sources.csv", delimiter=",", names=True
    )
    integrated_target = source_weights["weight"][target_local]
    if not np.isclose(raw_weight.sum(), integrated_target, atol=1e-12):
        raise RuntimeError("Grid and integrated target-band weights disagree")
    if float(projection_metadata["target_fraction_within_central_pair"]) <= 0.5:
        raise RuntimeError("The separated branch is not dominated by the target band")
    if float(projection_metadata["target_c3_relative_error"]) > 1e-12:
        raise RuntimeError("Exported Fig. 4 weight is not C3 covariant")
    print(output_stem.with_suffix(".pdf").resolve())
    print(
        "target_fraction_within_central_pair=",
        projection_metadata["target_fraction_within_central_pair"],
    )
    print(
        "target_c3_raw/restored=",
        projection_metadata["target_c3_raw_relative_error"],
        projection_metadata["target_c3_relative_error"],
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    root = Path(__file__).resolve().parents[1]
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=root / "results" / "paper_q20_converged" / "fig4",
    )
    parser.add_argument("--output-stem", type=Path, default=None)
    arguments = parser.parse_args()
    stem = arguments.output_stem or arguments.data_dir / "ideal_component_projection"
    main(arguments.data_dir, stem)
