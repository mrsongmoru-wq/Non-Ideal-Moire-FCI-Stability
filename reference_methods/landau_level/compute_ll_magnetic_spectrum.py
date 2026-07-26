#!/usr/bin/env python3
"""Export the fixed-flux tMBG magnetic spectrum in a Landau-level basis."""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
from joblib import Parallel, delayed
from scipy import linalg

from tmbg_landau_level import TMBGLandauLevelModel


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/ll_reference/magnetic_spectrum.csv"),
    )
    parser.add_argument("--p", type=int, default=1)
    parser.add_argument("--q", type=int, default=10)
    parser.add_argument("--nll", type=int, default=200)
    parser.add_argument("--nk1", type=int, default=20)
    parser.add_argument("--nk2", type=int, default=8)
    parser.add_argument(
        "--central-window",
        type=int,
        default=60,
        help="Number of levels retained on either side of the spectrum center.",
    )
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--theta", type=float, default=1.04)
    parser.add_argument("--kappa", type=float, default=0.70)
    parser.add_argument("--hbn", type=float, default=-30.0)
    parser.add_argument("--displacement", type=float, default=50.0)
    parser.add_argument(
        "--quick",
        action="store_true",
        help="Run a small installation smoke test; its values are not converged.",
    )
    args = parser.parse_args()
    if args.quick:
        args.nll = 20
        args.nk1 = 4
        args.nk2 = 2
        args.central_window = 10
        args.jobs = 1
    return args


def main() -> None:
    args = parse_args()
    model = TMBGLandauLevelModel(
        valley=1,
        signB=1,
        theta_d=args.theta,
        kappa=args.kappa,
        hBN_pot=args.hbn,
        V_pot=args.displacement,
        LayerhBN=0,
    )
    dimension = 6 * (args.nll + 1) * args.p
    center = dimension // 2
    first = center - args.central_window
    last = center + args.central_window - 1
    if first < 0 or last >= dimension:
        raise ValueError(
            f"Central window [{first}, {last}] exceeds basis dimension {dimension}"
        )

    def solve(i1: int, i2: int):
        hamiltonian = model.get_hamiltonian(
            i1 / args.nk1,
            i2 / args.nk2,
            0.0,
            0.0,
            args.p,
            args.q,
            args.nll,
        )
        values = linalg.eigh(
            hamiltonian,
            subset_by_index=[first, last],
            eigvals_only=True,
            driver="evr",
            check_finite=False,
        )
        return i1, i2, values

    started = time.time()
    points = [(i1, i2) for i1 in range(args.nk1) for i2 in range(args.nk2)]
    if args.jobs == 1:
        results = [solve(*point) for point in points]
    else:
        results = Parallel(n_jobs=args.jobs, prefer="processes")(
            delayed(solve)(*point) for point in points
        )

    rows = []
    for i1, i2, values in results:
        for offset, energy in enumerate(values):
            rows.append(
                (
                    i1 / args.nk1,
                    i2 / args.nk2,
                    offset + 1,
                    first + offset + 1,
                    energy,
                )
            )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(
        args.output,
        np.asarray(rows),
        delimiter=",",
        header="k1,k2,central_band_1based,full_basis_band_1based,energy_meV",
        comments="",
        fmt=["%.12g", "%.12g", "%d", "%d", "%.12g"],
    )
    metadata = args.output.with_suffix(".txt")
    magnetic_length, _, field_tesla = model.get_magnetic_strength(args.p, args.q)
    metadata.write_text(
        "\n".join(
            [
                f"output={args.output}",
                f"LL_input_flux={args.p}/{args.q}",
                f"physical_moire_flux={args.p}/{2 * args.q}",
                f"NLL={args.nll}",
                f"mesh={args.nk1}x{args.nk2}",
                f"central_bands={2 * args.central_window}",
                f"basis_dimension={dimension}",
                f"magnetic_length_angstrom={magnetic_length:.12g}",
                f"magnetic_field_tesla={field_tesla:.12g}",
                f"quick_smoke_test={args.quick}",
                f"elapsed_seconds={time.time() - started:.3f}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    print(args.output.resolve())
    print(metadata.resolve())


if __name__ == "__main__":
    main()
