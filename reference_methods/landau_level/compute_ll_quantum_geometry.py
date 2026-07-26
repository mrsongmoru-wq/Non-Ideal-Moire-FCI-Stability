#!/usr/bin/env python3
"""Compute tMBG magnetic-subband geometry in a Landau-level basis.

The band labels follow the frozen 120-band spectrum
used by the paper example: band 1 is the first of the 120 central levels.
Both rank-one and multiband projectors are supported.  The calculation retains
the full layer/sublattice/Landau-orbital basis and includes the analytic
Landau-orbital transfer matrices, Landau-gauge phase, background Berry phase,
and physical Cartesian trace condition.
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
from joblib import Parallel, delayed
from scipy import linalg


DEFAULT_OUTPUT = Path("results/ll_reference/quantum_geometry.npz")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--band-start", type=int, default=82,
                        help="First selected central-spectrum band (one based).")
    parser.add_argument("--nbands", type=int, default=1)
    parser.add_argument("--nk1", type=int, default=20)
    parser.add_argument("--nk2", type=int, default=8)
    parser.add_argument("--nll", type=int, default=200)
    parser.add_argument("--jobs", type=int, default=1,
                        help="Concurrent eigensystems; BLAS threads should be one.")
    parser.add_argument("--parallel-prefer", choices=("processes", "threads"),
                        default="processes")
    parser.add_argument("--geometry-backend", choices=("library", "block"),
                        default="block",
                        help="block is algebraically identical for p=1 and uses less memory; library supports general p.")
    parser.add_argument("--central-window", type=int, default=60,
                        help="Half-width used by the frozen 120-band spectrum.")
    parser.add_argument("--p", type=int, default=1,
                        help="LL-code flux numerator (physical flux is p/(2q)).")
    parser.add_argument("--q", type=int, default=10,
                        help="LL-code flux denominator (physical flux is p/(2q)).")
    parser.add_argument("--theta", type=float, default=1.04)
    parser.add_argument("--kappa", type=float, default=0.7)
    parser.add_argument("--hbn", type=float, default=-30.0)
    parser.add_argument("--displacement", type=float, default=50.0)
    parser.add_argument("--state-cache", type=Path)
    parser.add_argument("--reuse-state-cache", action="store_true")
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
        args.central_window = 30
        args.band_start = 30
        args.nbands = 1
        args.jobs = 1
    if args.output.suffix != ".npz":
        parser.error("--output must end in .npz")
    return args


def selected_eigensystem(model, args: argparse.Namespace):
    rank = args.nbands
    hdim = 6 * (args.nll + 1) * args.p
    center = 3 * (args.nll + 1) * args.p
    first_central_zero = args.band_start - 1
    global_first = center - args.central_window + first_central_zero
    global_last = global_first + rank - 1
    if global_first < 0 or global_last >= hdim:
        raise ValueError(
            "Selected central-spectrum interval lies outside the LL basis: "
            f"global interval [{global_first}, {global_last}], dimension {hdim}"
        )

    if args.reuse_state_cache:
        if args.state_cache is None:
            raise ValueError("--reuse-state-cache requires --state-cache")
        cached = np.load(args.state_cache)
        states = cached["states"]
        energies = cached["energies"]
        expected = (args.nk1, args.nk2, hdim, rank)
        if states.shape != expected:
            raise ValueError(f"State-cache shape {states.shape} != {expected}")
        return energies, states

    states = np.empty((args.nk1, args.nk2, hdim, rank), dtype=np.complex128)
    energies = np.empty((args.nk1, args.nk2, rank), dtype=np.float64)
    total = args.nk1 * args.nk2
    started = time.time()
    count = 0
    def solve_one(i1: int, i2: int):
        hamiltonian = model.get_hamiltonian(
            i1 / args.nk1,
            i2 / args.nk2,
            0.0,
            0.0,
            args.p,
            args.q,
            args.nll,
        )
        values, vectors = linalg.eigh(
            hamiltonian,
            subset_by_index=[global_first, global_last],
            driver="evr",
            check_finite=False,
        )
        return i1, i2, values, vectors

    points = [(i1, i2) for i1 in range(args.nk1) for i2 in range(args.nk2)]
    batch_size = max(args.jobs, total // 10)
    for batch_start in range(0, total, batch_size):
        batch = points[batch_start:batch_start + batch_size]
        if args.jobs == 1:
            results = [solve_one(i1, i2) for i1, i2 in batch]
        else:
            results = Parallel(n_jobs=args.jobs, prefer=args.parallel_prefer)(
                delayed(solve_one)(i1, i2) for i1, i2 in batch
            )
        for i1, i2, values, vectors in results:
            energies[i1, i2, :] = values
            states[i1, i2, :, :] = vectors
            count += 1
        print(
            f"eigensystems {count}/{total}; elapsed={time.time()-started:.1f}s",
            flush=True,
        )

    if args.state_cache is not None:
        args.state_cache.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            args.state_cache,
            states=states,
            energies=energies,
            band_start=np.int32(args.band_start),
            nbands=np.int32(args.nbands),
            nk1=np.int32(args.nk1),
            nk2=np.int32(args.nk2),
            nll=np.int32(args.nll),
            p=np.int32(args.p),
            q=np.int32(args.q),
        )
    return energies, states


def ll_transfer_block(Lmn, nll: int, z: complex) -> np.ndarray:
    """Return one LL form-factor block, matching get_transfer_matrix exactly."""
    block = np.zeros((nll + 1, nll + 1), dtype=np.complex128)
    # The library intentionally uses range(NLL), leaving the cutoff row and
    # column zero.  Keeping this detail is essential for exact validation.
    for n in range(nll):
        for m in range(nll):
            block[m, n] = Lmn(z, n, m)
    return block


def block_ll_geometry(model, states: np.ndarray, args: argparse.Namespace, Lmn):
    """Memory-efficient p=1 form of get_chern_number_non_ab."""
    if args.p != 1:
        raise ValueError("The block backend currently requires LL numerator p=1")
    rank = args.nbands
    magnetic_length, _, _ = model.get_magnetic_strength(args.p, args.q)
    dkx = model.Qx / args.nk1
    dky = model.Qy / args.q / args.nk2
    zx = dkx * magnetic_length / np.sqrt(2.0)
    zy = 1j * dky * magnetic_length / np.sqrt(2.0)
    zxy = (-dkx + 1j * dky) * magnetic_length / np.sqrt(2.0)
    tx = ll_transfer_block(Lmn, args.nll, zx)
    ty = ll_transfer_block(Lmn, args.nll, zy)
    txy = ll_transfer_block(Lmn, args.nll, zxy)

    def apply(block: np.ndarray, vectors: np.ndarray) -> np.ndarray:
        shaped = vectors.reshape(3, 2, args.nll + 1, rank)
        result = np.einsum("mn,abnr->abmr", block, shaped, optimize=True)
        return result.reshape(6 * (args.nll + 1), rank)

    berry_curvature = np.zeros((args.nk1, args.nk2), dtype=np.float64)
    metric_trace = np.zeros_like(berry_curvature)
    determinant_condition = np.zeros_like(berry_curvature)
    wilson_spread_terms = []
    positive_flux = 0.0
    negative_flux = 0.0
    total_flux = 0.0
    background = np.exp(-1j * 2.0 * np.pi / args.p / args.nk1 / args.nk2)

    for i1 in range(args.nk1):
        for i2 in range(args.nk2):
            u = states[i1, i2]
            up1 = states[(i1 + 1) % args.nk1, i2]
            up2 = states[i1, (i2 + 1) % args.nk2]
            up12 = states[(i1 + 1) % args.nk1, (i2 + 1) % args.nk2]

            vx_x = u.conj().T @ apply(tx.conj().T, up1)
            vy_x = up1.conj().T @ apply(ty, up12)
            vx_y = up2.conj().T @ apply(tx.conj().T, up12)
            vy_y = u.conj().T @ apply(ty, up2)
            wilson = (
                vx_x @ vy_x @ np.linalg.inv(vx_y) @ np.linalg.inv(vy_y)
                * background
            )
            flux = float(np.imag(np.log(np.linalg.det(wilson))))
            total_flux += flux
            positive_flux += max(0.0, flux)
            negative_flux += min(0.0, flux)
            berry_curvature[i1, i2] = flux / (dkx * dky)
            wilson_spread_terms.append(
                np.sum(np.imag(np.log(np.linalg.eigvals(wilson))) ** 2)
                / (4.0 * np.pi**2)
            )

            phase = np.exp(1j * magnetic_length**2 * dky * i1 * dkx)
            mx = up1.conj().T @ apply(tx, u)
            my = (u.conj().T @ apply(ty, up2)) * phase
            mixed = (up1.conj().T @ apply(txy, up2)) * phase
            norm_x = np.trace(mx @ mx.conj().T).real
            norm_y = np.trace(my @ my.conj().T).real
            norm_mixed = np.trace(mixed @ mixed.conj().T).real
            gxx = (rank - norm_x) / dkx**2
            gyy = (rank - norm_y) / dky**2
            gxy = (rank - norm_x - norm_y + norm_mixed) / (2.0 * dkx * dky)
            metric_trace[i1, i2] = gxx + gyy
            determinant_condition[i1, i2] = (
                gxx * gyy - gxy**2 - abs(berry_curvature[i1, i2]) ** 2 / 4.0
            )

    trace_violation = metric_trace - np.abs(berry_curvature)
    plaquette_area = dkx * dky
    return (
        total_flux / (2.0 * np.pi),
        berry_curvature,
        positive_flux / (2.0 * np.pi),
        negative_flux / (2.0 * np.pi),
        float(np.sqrt(np.sum(wilson_spread_terms))),
        metric_trace,
        float(np.sum(trace_violation) * plaquette_area),
        determinant_condition,
        float(np.sum(determinant_condition)),
    )


def main() -> None:
    args = parse_args()
    from tmbg_landau_level import (
        Lmn,
        TMBGLandauLevelModel,
        compute_non_abelian_geometry,
    )

    model = TMBGLandauLevelModel(
        valley=1,
        signB=1,
        theta_d=args.theta,
        kappa=args.kappa,
        hBN_pot=args.hbn,
        V_pot=args.displacement,
        LayerhBN=0,
    )
    energies, states = selected_eigensystem(model, args)
    # The selected-state array starts at local band zero.  All geometric
    # operations below are the original LL-library implementation.
    model.eigsvec2D = states
    geometry_started = time.time()
    geometry_function = (
        (lambda: block_ll_geometry(model, states, args, Lmn))
        if args.geometry_backend == "block"
        else (lambda: compute_non_abelian_geometry(
            model,
            0,
            args.nbands,
            args.nk1,
            args.nk2,
            args.p,
            args.q,
            args.nll,
        ))
    )
    (
        chern,
        berry_curvature,
        positive_flux,
        negative_flux,
        wilson_spread,
        metric_trace,
        eta,
        determinant_condition,
        integrated_determinant_condition,
    ) = geometry_function()
    dkx = model.Qx / args.nk1
    dky = model.Qy / args.q / args.nk2
    plaquette_area = dkx * dky
    integrated_trace = float(np.sum(metric_trace) * plaquette_area)
    integrated_abs_curvature = float(
        np.sum(np.abs(berry_curvature)) * plaquette_area
    )
    trace_violation = metric_trace - np.abs(berry_curvature)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.output,
        energies=energies,
        berry_curvature=berry_curvature,
        metric_trace=metric_trace,
        trace_violation=trace_violation,
        determinant_condition=determinant_condition,
        chern=np.float64(chern),
        eta=np.float64(eta),
        integrated_trace=np.float64(integrated_trace),
        integrated_abs_curvature=np.float64(integrated_abs_curvature),
        positive_flux=np.float64(positive_flux),
        negative_flux=np.float64(negative_flux),
        wilson_spread=np.float64(wilson_spread),
        integrated_determinant_condition=np.float64(
            integrated_determinant_condition
        ),
        plaquette_area=np.float64(plaquette_area),
        band_start=np.int32(args.band_start),
        nbands=np.int32(args.nbands),
        nk1=np.int32(args.nk1),
        nk2=np.int32(args.nk2),
        nll=np.int32(args.nll),
        p=np.int32(args.p),
        q=np.int32(args.q),
        physical_p=np.int32(args.p),
        physical_q=np.int32(2 * args.q),
    )
    summary_path = args.output.with_suffix(".txt")
    summary_path.write_text(
        "\n".join(
            [
                f"output={args.output}",
                f"LL_input_flux={args.p}/{args.q}",
                f"physical_moire_flux={args.p}/{2*args.q}",
                f"band_interval_1based={args.band_start}:{args.band_start+args.nbands-1}",
                f"projector_rank={args.nbands}",
                f"mesh={args.nk1}x{args.nk2}",
                f"NLL={args.nll}",
                f"geometry_backend={args.geometry_backend}",
                f"quick_smoke_test={args.quick}",
                f"mean_energy_meV={float(np.mean(energies)):.12g}",
                f"chern={chern:.12g}",
                f"eta={eta:.12g}",
                f"integrated_trace={integrated_trace:.12g}",
                f"integrated_abs_curvature={integrated_abs_curvature:.12g}",
                f"trace_violation_min={float(np.min(trace_violation)):.12g}",
                f"trace_violation_max={float(np.max(trace_violation)):.12g}",
                f"geometry_seconds={time.time()-geometry_started:.3f}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    print(summary_path.read_text(encoding="utf-8"), flush=True)


if __name__ == "__main__":
    main()
