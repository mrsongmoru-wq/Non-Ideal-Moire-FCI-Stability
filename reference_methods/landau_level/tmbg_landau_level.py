"""Minimal Landau-level expansion of the tMBG continuum Hamiltonian.

The module contains only the pieces used by the repository's linear LL
comparison workflow: the Hamiltonian, analytic Landau-orbital form factors,
and the full-matrix non-Abelian geometry backend.  The production p=1
calculation uses the algebraically equivalent block backend in
``compute_ll_quantum_geometry.py`` to reduce memory use.
"""

from __future__ import annotations

import ctypes

import numpy as np
from numba import njit
from numba.extending import get_cython_function_address


_DOUBLE = ctypes.c_double
_LONG = ctypes.c_long
_laguerre_address = get_cython_function_address(
    "scipy.special.cython_special",
    "__pyx_fuse_1_1eval_genlaguerre",
)
_laguerre_type = ctypes.CFUNCTYPE(_DOUBLE, _LONG, _DOUBLE, _DOUBLE)
_eval_genlaguerre = _laguerre_type(_laguerre_address)

_gammaln_address = get_cython_function_address(
    "scipy.special.cython_special",
    "gammaln",
)
_gammaln_type = ctypes.CFUNCTYPE(_DOUBLE, _DOUBLE)
_gammaln = _gammaln_type(_gammaln_address)

PLANCK_CONSTANT = 6.626e-34
ELEMENTARY_CHARGE = 1.601e-19
FLUX_QUANTUM = PLANCK_CONSTANT / ELEMENTARY_CHARGE


@njit
def Lmn(z: complex, m: int, n: int) -> complex:
    """Landau-orbital form factor in the convention used by the LL model."""
    z_squared = np.abs(z) ** 2
    if m < n:
        return (
            (-1) ** (n - m)
            * np.exp(
                (_gammaln(m + 1) - _gammaln(n + 1)) / 2
                + np.log(z) * (n - m)
                - z_squared / 2
            )
            * _eval_genlaguerre(m, n - m, z_squared)
        )
    return (
        np.exp(
            (_gammaln(n + 1) - _gammaln(m + 1)) / 2
            + np.log(np.conj(z)) * (m - n)
            - z_squared / 2
        )
        * _eval_genlaguerre(n, m - n, z_squared)
    )


@njit
def _positive_valley_hamiltonian(
    k1,
    k2,
    kphi1,
    kphi2,
    p,
    q,
    nll,
    t1,
    hbarv3,
    hbarv4,
    kappa,
    w,
    eps,
    magnetic_length,
    k_theta,
    delta,
    qx,
    qy,
    displacement,
    hbn_potential,
    hbn_layers,
    field_sign,
):
    layer_potential = np.array([-1.0, 0.0, 1.0])
    w1 = kappa * w
    w2 = w
    phase = -field_sign * 2.0 * np.pi / 3.0
    phase_factor = np.exp(1j * phase)
    tunneling0 = np.array([[w1, w2], [w2, w1]], dtype=np.complex128)
    tunneling1 = (
        np.array(
            [[w1, w2 * phase_factor], [w2 / phase_factor, w1]],
            dtype=np.complex128,
        )
        * np.exp(1j * 2.0 * np.pi * kphi1)
    )
    tunneling2 = (
        np.array(
            [[w1, w2 / phase_factor], [w2 * phase_factor, w1]],
            dtype=np.complex128,
        )
        * np.exp(1j * 2.0 * np.pi * kphi2)
    )
    tunneling = (tunneling1, tunneling2)

    dimension = 6 * (nll + 1) * p
    indices = np.arange(dimension).reshape((3, 2, nll + 1, p))
    hamiltonian = np.zeros((dimension, dimension), np.complex128)
    cone_shift = np.array([1.0, 0.0, 0.0])

    for magnetic_index in range(p):
        for layer in range(3):
            for n in range(nll):
                hamiltonian[
                    indices[layer, 0, n, magnetic_index],
                    indices[layer, 1, n + 1, magnetic_index],
                ] = (
                    -field_sign
                    * np.sqrt(2.0)
                    * eps
                    / k_theta
                    / magnetic_length
                    * np.sqrt(n + 1)
                )
                hamiltonian[
                    indices[layer, 0, n, magnetic_index],
                    indices[layer, 1, n, magnetic_index],
                ] = -field_sign * 1j * eps * cone_shift[layer]
                hamiltonian[
                    indices[layer, 0, n, magnetic_index],
                    indices[layer, 0, n, magnetic_index],
                ] = (
                    hbn_potential * hbn_layers[layer] / 4.0
                    + displacement * layer_potential[layer] / 4.0
                )
                hamiltonian[
                    indices[layer, 1, n, magnetic_index],
                    indices[layer, 1, n, magnetic_index],
                ] = (
                    -hbn_potential * hbn_layers[layer] / 4.0
                    + displacement * layer_potential[layer] / 4.0
                )
            hamiltonian[
                indices[layer, 0, nll, magnetic_index],
                indices[layer, 0, nll, magnetic_index],
            ] = -5.0e5

    for magnetic_index in range(p):
        for n in range(nll):
            hamiltonian[
                indices[2, 0, n + 1, magnetic_index],
                indices[1, 0, n, magnetic_index],
            ] = (
                field_sign
                * np.sqrt(2.0)
                * hbarv3
                / magnetic_length
                * np.sqrt(n + 1)
            )
            hamiltonian[
                indices[2, 1, n + 1, magnetic_index],
                indices[1, 1, n, magnetic_index],
            ] = (
                field_sign
                * np.sqrt(2.0)
                * hbarv3
                / magnetic_length
                * np.sqrt(n + 1)
            )
            hamiltonian[
                indices[2, 1, n, magnetic_index],
                indices[1, 0, n + 1, magnetic_index],
            ] = (
                field_sign
                * np.sqrt(2.0)
                * hbarv4
                / magnetic_length
                * np.sqrt(n + 1)
            )
            hamiltonian[
                indices[2, 0, n, magnetic_index],
                indices[1, 1, n, magnetic_index],
            ] = t1

    kx = k1 * qx
    ky = k2 * qy / q
    reciprocal_x = np.array([1, -1])
    reciprocal_y = np.array([1, 1])
    form_factor_argument = magnetic_length / np.sqrt(2.0) * np.array(
        [
            field_sign * reciprocal_x[index] * qx
            + 1j * reciprocal_y[index] * qy
            for index in range(2)
        ]
    )
    for magnetic_index in range(p):
        translated_index = np.array(
            [
                (magnetic_index - reciprocal_x[index]) % p
                for index in range(2)
            ]
        )
        for sublattice1 in range(2):
            for sublattice2 in range(2):
                for m in range(nll + 1):
                    hamiltonian[
                        indices[0, sublattice1, m, magnetic_index],
                        indices[1, sublattice2, m, magnetic_index],
                    ] += tunneling0[sublattice1, sublattice2]
                    for n in range(nll + 1):
                        for index in range(2):
                            hamiltonian[
                                indices[
                                    0,
                                    sublattice1,
                                    m,
                                    translated_index[index],
                                ],
                                indices[1, sublattice2, n, magnetic_index],
                            ] += (
                                tunneling[index][sublattice1, sublattice2]
                                * Lmn(form_factor_argument[index], m, n)
                                * np.exp(
                                    field_sign
                                    * 1j
                                    * ky
                                    * reciprocal_x[index]
                                    * delta
                                )
                                * np.exp(
                                    -field_sign
                                    * 1j
                                    * reciprocal_y[index]
                                    * kx
                                    * qy
                                    * magnetic_length**2
                                )
                                * np.exp(
                                    -field_sign
                                    * 1j
                                    * (
                                        reciprocal_y[index] * magnetic_index
                                        - reciprocal_x[index]
                                        * reciprocal_y[index]
                                        / 2.0
                                    )
                                    * 2.0
                                    * np.pi
                                    * q
                                    / p
                                )
                            )
    return hamiltonian + hamiltonian.T.conj()


@njit
def _negative_valley_hamiltonian(
    k1,
    k2,
    kphi1,
    kphi2,
    p,
    q,
    nll,
    t1,
    hbarv3,
    hbarv4,
    kappa,
    w,
    eps,
    magnetic_length,
    k_theta,
    delta,
    qx,
    qy,
    displacement,
    hbn_potential,
    hbn_layers,
    field_sign,
):
    """Negative-valley partner with the original LL ordering convention."""
    layer_potential = np.array([-1.0, 0.0, 1.0])
    w1 = kappa * w
    w2 = w
    phase = field_sign * 2.0 * np.pi / 3.0
    phase_factor = np.exp(1j * phase)
    tunneling0 = np.array([[w1, w2], [w2, w1]], dtype=np.complex128)
    tunneling1 = (
        np.array(
            [[w1, w2 * phase_factor], [w2 / phase_factor, w1]],
            dtype=np.complex128,
        )
        * np.exp(1j * 2.0 * np.pi * kphi1)
    )
    tunneling2 = (
        np.array(
            [[w1, w2 / phase_factor], [w2 * phase_factor, w1]],
            dtype=np.complex128,
        )
        * np.exp(1j * 2.0 * np.pi * kphi2)
    )
    tunneling = (tunneling1, tunneling2)

    dimension = 6 * (nll + 1) * p
    indices = np.arange(dimension).reshape((3, 2, nll + 1, p))
    hamiltonian = np.zeros((dimension, dimension), np.complex128)
    cone_shift = np.array([1.0, 0.0, 0.0])

    for magnetic_index in range(p):
        for layer in range(3):
            for n in range(nll):
                hamiltonian[
                    indices[layer, 0, n + 1, magnetic_index],
                    indices[layer, 1, n, magnetic_index],
                ] = (
                    field_sign
                    * np.sqrt(2.0)
                    * eps
                    / k_theta
                    / magnetic_length
                    * np.sqrt(n + 1)
                )
                hamiltonian[
                    indices[layer, 0, n, magnetic_index],
                    indices[layer, 1, n, magnetic_index],
                ] = -field_sign * 1j * eps * cone_shift[layer]
                hamiltonian[
                    indices[layer, 0, n, magnetic_index],
                    indices[layer, 0, n, magnetic_index],
                ] = (
                    hbn_potential * hbn_layers[layer] / 4.0
                    + displacement * layer_potential[layer] / 4.0
                )
                hamiltonian[
                    indices[layer, 1, n, magnetic_index],
                    indices[layer, 1, n, magnetic_index],
                ] = (
                    -hbn_potential * hbn_layers[layer] / 4.0
                    + displacement * layer_potential[layer] / 4.0
                )
            hamiltonian[
                indices[layer, 1, nll, magnetic_index],
                indices[layer, 1, nll, magnetic_index],
            ] = 5.0e5

    for magnetic_index in range(p):
        for n in range(nll):
            hamiltonian[
                indices[2, 0, n, magnetic_index],
                indices[1, 0, n + 1, magnetic_index],
            ] = (
                -field_sign
                * np.sqrt(2.0)
                * hbarv3
                / magnetic_length
                * np.sqrt(n + 1)
            )
            hamiltonian[
                indices[2, 1, n, magnetic_index],
                indices[1, 1, n + 1, magnetic_index],
            ] = (
                -field_sign
                * np.sqrt(2.0)
                * hbarv3
                / magnetic_length
                * np.sqrt(n + 1)
            )
            hamiltonian[
                indices[2, 1, n + 1, magnetic_index],
                indices[1, 0, n, magnetic_index],
            ] = (
                -field_sign
                * np.sqrt(2.0)
                * hbarv4
                / magnetic_length
                * np.sqrt(n + 1)
            )
            hamiltonian[
                indices[2, 0, n, magnetic_index],
                indices[1, 1, n, magnetic_index],
            ] = t1

    kx = k1 * qx
    ky = k2 * qy / q
    reciprocal_x = np.array([1, -1])
    reciprocal_y = np.array([1, 1])
    form_factor_argument = magnetic_length / np.sqrt(2.0) * np.array(
        [
            field_sign * reciprocal_x[index] * qx
            + 1j * reciprocal_y[index] * qy
            for index in range(2)
        ]
    )
    for magnetic_index in range(p):
        translated_index = np.array(
            [
                (magnetic_index - reciprocal_x[index]) % p
                for index in range(2)
            ]
        )
        for sublattice1 in range(2):
            for sublattice2 in range(2):
                for m in range(nll + 1):
                    hamiltonian[
                        indices[0, sublattice1, m, magnetic_index],
                        indices[1, sublattice2, m, magnetic_index],
                    ] += tunneling0[sublattice1, sublattice2]
                    for n in range(nll + 1):
                        for index in range(2):
                            hamiltonian[
                                indices[
                                    0,
                                    sublattice1,
                                    m,
                                    translated_index[index],
                                ],
                                indices[1, sublattice2, n, magnetic_index],
                            ] += (
                                tunneling[index][sublattice1, sublattice2]
                                * Lmn(form_factor_argument[index], m, n)
                                * np.exp(
                                    field_sign
                                    * 1j
                                    * ky
                                    * reciprocal_x[index]
                                    * delta
                                )
                                * np.exp(
                                    -field_sign
                                    * 1j
                                    * reciprocal_y[index]
                                    * kx
                                    * qy
                                    * magnetic_length**2
                                )
                                * np.exp(
                                    -field_sign
                                    * 1j
                                    * (
                                        reciprocal_y[index] * magnetic_index
                                        - reciprocal_x[index]
                                        * reciprocal_y[index]
                                        / 2.0
                                    )
                                    * 2.0
                                    * np.pi
                                    * q
                                    / p
                                )
                            )
    return hamiltonian + hamiltonian.T.conj()


def landau_transfer_matrix(nll: int, p: int, z: complex) -> np.ndarray:
    """Return the full analytic momentum-transfer matrix in the LL basis."""
    dimension = 6 * (nll + 1) * p
    indices = np.arange(dimension).reshape((3, 2, nll + 1, p))
    transfer = np.zeros((dimension, dimension), np.complex128)
    # This reproduces the paper calculation's cutoff convention: the last
    # orbital is a regulator and is not included in geometric transfer.
    for n in range(nll):
        for m in range(nll):
            for magnetic_index in range(p):
                for sublattice in range(2):
                    for layer in range(3):
                        transfer[
                            indices[layer, sublattice, m, magnetic_index],
                            indices[layer, sublattice, n, magnetic_index],
                        ] = Lmn(z, n, m)
    return transfer


def _plaquette_flux(
    point,
    states,
    transfer_x,
    transfer_y,
    p,
    rank,
    nk1,
    nk2,
):
    i1, i2 = point
    u00 = states[i1, i2]
    u10 = states[(i1 + 1) % nk1, i2]
    u01 = states[i1, (i2 + 1) % nk2]
    u11 = states[(i1 + 1) % nk1, (i2 + 1) % nk2]
    link_x_bottom = u00.conj().T @ transfer_x.conj().T @ u10
    link_y_right = u10.conj().T @ transfer_y @ u11
    link_x_top = u01.conj().T @ transfer_x.conj().T @ u11
    link_y_left = u00.conj().T @ transfer_y @ u01
    background = np.exp(-1j * 2.0 * np.pi / p / nk1 / nk2)
    wilson = (
        link_x_bottom
        @ link_y_right
        @ np.linalg.inv(link_x_top)
        @ np.linalg.inv(link_y_left)
        * background
    )
    flux = float(np.imag(np.log(np.linalg.det(wilson))))
    eigenphases = np.imag(np.log(np.linalg.eigvals(wilson)))
    return flux, eigenphases


def compute_non_abelian_geometry(
    model,
    first_band: int,
    rank: int,
    nk1: int,
    nk2: int,
    p: int,
    q: int,
    nll: int,
):
    """Compute non-Abelian curvature and metric for a selected LL projector."""
    states = model.eigsvec2D[:, :, :, first_band : first_band + rank]
    magnetic_length, _, _ = model.get_magnetic_strength(p, q)
    dkx = model.Qx / nk1
    dky = model.Qy / q / nk2
    transfer_x = landau_transfer_matrix(
        nll, p, dkx * magnetic_length / np.sqrt(2.0)
    )
    transfer_y = landau_transfer_matrix(
        nll, p, 1j * dky * magnetic_length / np.sqrt(2.0)
    )
    transfer_mixed = landau_transfer_matrix(
        nll, p, (-dkx + 1j * dky) * magnetic_length / np.sqrt(2.0)
    )

    curvature = np.zeros((nk1, nk2), dtype=float)
    metric_trace = np.zeros_like(curvature)
    determinant_condition = np.zeros_like(curvature)
    total_flux = 0.0
    positive_flux = 0.0
    negative_flux = 0.0
    wilson_spread_terms = []

    for i1 in range(nk1):
        for i2 in range(nk2):
            flux, eigenphases = _plaquette_flux(
                (i1, i2),
                states,
                transfer_x,
                transfer_y,
                p,
                rank,
                nk1,
                nk2,
            )
            total_flux += flux
            positive_flux += max(flux, 0.0)
            negative_flux += min(flux, 0.0)
            wilson_spread_terms.append(
                np.sum(eigenphases**2) / (4.0 * np.pi**2)
            )
            curvature[i1, i2] = flux / (dkx * dky)

            u = states[i1, i2]
            u_x = states[(i1 + 1) % nk1, i2]
            u_y = states[i1, (i2 + 1) % nk2]
            phase = np.exp(1j * magnetic_length**2 * dky * i1 * dkx)
            overlap_x = u_x.conj().T @ transfer_x @ u
            overlap_y = u.conj().T @ transfer_y @ u_y * phase
            overlap_mixed = (
                u_x.conj().T @ transfer_mixed @ u_y * phase
            )
            norm_x = np.trace(overlap_x @ overlap_x.conj().T).real
            norm_y = np.trace(overlap_y @ overlap_y.conj().T).real
            norm_mixed = np.trace(
                overlap_mixed @ overlap_mixed.conj().T
            ).real
            metric_xx = (rank - norm_x) / dkx**2
            metric_yy = (rank - norm_y) / dky**2
            metric_xy = (
                rank - norm_x - norm_y + norm_mixed
            ) / (2.0 * dkx * dky)
            metric_trace[i1, i2] = metric_xx + metric_yy
            determinant_condition[i1, i2] = (
                metric_xx * metric_yy
                - metric_xy**2
                - abs(curvature[i1, i2]) ** 2 / 4.0
            )

    trace_excess = metric_trace - np.abs(curvature)
    plaquette_area = dkx * dky
    return (
        total_flux / (2.0 * np.pi),
        curvature,
        positive_flux / (2.0 * np.pi),
        negative_flux / (2.0 * np.pi),
        float(np.sqrt(np.sum(wilson_spread_terms))),
        metric_trace,
        float(np.sum(trace_excess) * plaquette_area),
        determinant_condition,
        float(np.sum(determinant_condition)),
    )


class TMBGLandauLevelModel:
    """tMBG continuum Hamiltonian in a finite Landau-orbital basis."""

    def __init__(
        self,
        valley=1,
        signB=1,
        theta_d=1.4,
        kappa=0.6,
        hBN_pot=30.0,
        V_pot=50.0,
        LayerhBN=0,
    ):
        self.valley = valley
        self.signB = signB
        self.V_pot = V_pot
        self.a = 1.43
        self.hBN_pot = hBN_pot
        self.hBN_pot_layer = np.zeros(3)
        self.hBN_pot_layer[LayerhBN] = 1.0
        self.w = 110.0
        self.theta = np.deg2rad(theta_d)
        self.kappa = kappa
        self.k_theta = (
            8.0
            * np.pi
            / 3.0
            / np.sqrt(3.0)
            / self.a
            * np.sin(self.theta / 2.0)
        )
        self.A1 = (
            4.0
            * np.pi
            / 3.0
            / self.k_theta
            * np.array([np.sqrt(3.0) / 2.0, 1.0 / 2.0])
        )
        self.A2 = (
            4.0
            * np.pi
            / 3.0
            / self.k_theta
            * np.array([-np.sqrt(3.0) / 2.0, 1.0 / 2.0])
        )
        self.Qx = self.k_theta * np.sqrt(3.0) / 2.0
        self.Qy = self.k_theta * 3.0 / 2.0
        self.b1 = np.array([self.Qx, self.Qy])
        self.b2 = np.array([-self.Qx, self.Qy])
        gamma0 = 2610.0
        self.t1 = 361.0
        gamma3 = 140.0
        gamma4 = 283.0
        lattice_constant = self.a * np.sqrt(3.0)
        hbarv = np.sqrt(3.0) / 2.0 * lattice_constant * gamma0
        self.hbarv3 = np.sqrt(3.0) / 2.0 * lattice_constant * gamma3
        self.hbarv4 = np.sqrt(3.0) / 2.0 * lattice_constant * gamma4
        self.eps = hbarv * self.k_theta

    def get_magnetic_strength(self, p: int, q: int):
        magnetic_length = np.sqrt(
            2.0 * np.pi * q / p / self.Qx / self.Qy
        )
        delta = self.Qx * magnetic_length**2
        field_tesla = (
            self.signB
            * FLUX_QUANTUM
            / (magnetic_length * 1.0e-10) ** 2
            / (2.0 * np.pi)
        )
        return magnetic_length, delta, field_tesla

    def get_hamiltonian(
        self,
        k1,
        k2,
        kphi1,
        kphi2,
        p,
        q,
        nll,
    ):
        magnetic_length, delta, _ = self.get_magnetic_strength(p, q)
        builder = (
            _positive_valley_hamiltonian
            if self.valley == 1
            else _negative_valley_hamiltonian
        )
        return builder(
            k1,
            k2,
            kphi1,
            kphi2,
            p,
            q,
            nll,
            self.t1,
            self.hbarv3,
            self.hbarv4,
            self.kappa,
            self.w,
            self.eps,
            magnetic_length,
            self.k_theta,
            delta,
            self.Qx,
            self.Qy,
            self.V_pot,
            self.hBN_pot,
            self.hBN_pot_layer,
            self.signB,
        )


# Compatibility aliases for private notebooks that used the historical names.
LLexpModeltMBG = TMBGLandauLevelModel


def get_transfer_matrix(nll, p, z, sizeH=None):
    """Compatibility wrapper for the historical four-argument function."""
    transfer = landau_transfer_matrix(nll, p, z)
    if sizeH is not None and transfer.shape != (sizeH, sizeH):
        raise ValueError(
            f"Requested sizeH={sizeH}, but LL basis has shape {transfer.shape}"
        )
    return transfer


get_chern_number_non_ab = compute_non_abelian_geometry
