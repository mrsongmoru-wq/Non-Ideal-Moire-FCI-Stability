# Non-Ideal FCI Stability Check

This repository provides an independent hybrid-Wannier workflow for testing
the single-particle ingredients of fractional Chern-insulator stability in a
non-ideal magnetic flat band.  Starting from a continuum Hamiltonian, the code
constructs magnetic subbands, exports their spectrum, evaluates the Berry
curvature and Cartesian quantum metric of an isolated target subband, records
the trace-condition result and overlap diagnostics, and resolves the target
state into its momentum-dependent ideal component.

The implemented and validated case is twisted monolayer-bilayer graphene
(tMBG) at flux `p/q=1/20`.  No Landau-level expansion, fitted reference
wavefunction, or external geometry data are used.

## Associated article

This code accompanies and was developed for the numerical analysis in the
arXiv manuscript:

> Moru Song and Kai Chang, *Fractional Chern Insulators Transition in
> Non-ideal Flat Bands of Twisted Mono-bilayer Graphene* (2026).

In particular, the repository provides the independent hybrid-Wannier
calculation of the weak-field magnetic spectrum, the multiband quantum
geometry and trace condition of the isolated magnetic subband, and the
momentum-resolved ideal-component distribution discussed in that work.  An
arXiv identifier and direct article link will be added here once assigned.

## Origin, license, and scope

This project is based on Xiaoyu Wang's
[`Hofstadter-TBG`](https://github.com/xywang2017/Hofstadter-TBG) repository and
the Wang--Vafek hybrid-Wannier magnetic-translation construction.  It adapts
that construction to the tMBG continuum Hamiltonian and adds the generalized
nonorthogonal magnetic solve, full field-coupling matrix elements, common-basis
lifting, multiband quantum geometry, overlap diagnostics, and ideal-component
projection.

The code is distributed under the GNU General Public License v3.0 inherited
from the upstream project; see `LICENSE`.

The architecture can be extended to other continuum-model moire systems.  A
new backend must supply the continuum Hamiltonian, direct and reciprocal
lattices, momentum-linear field coupling, local orbital dimension, and
reciprocal-space sewing maps.  The same pipeline can then generate magnetic
spectra, topology, quantum geometry, trace-condition diagnostics, and the
single-particle inputs for subsequent FCI stability calculations.  At present,
only the tMBG backend and the case below have been numerically validated.

## Physical method

The calculation is one linear sequence.

1. The tMBG continuum Hamiltonian is diagonalized on a uniform zero-field
   Brillouin-zone mesh.  Non-Abelian parallel transport along the first
   reciprocal direction constructs hybrid-Wannier states localized along the
   corresponding direct-lattice direction.
2. Magnetic translations of the retained zero-field bands form a
   nonorthogonal basis at flux `p/q`.  The overlap matrix `O(k)` and magnetic
   Hamiltonian `H(k)` are evaluated in that same basis.  Restriction to the
   physical rank followed by `O(k)^(-1/2)` gives the orthonormal magnetic
   Hamiltonian and its complete subband spectrum.
3. The isolated target magnetic subband is put first in a local multiband
   correction frame `Phi(k)`.  Nearby subbands stabilize the frame but are not
   merged into the target projector.  For neighboring momenta,

       M(k,k+dk) = Phi(k)' T(dk) Phi(k+dk)

   is lifted to a common plane-wave Hilbert space.  If `M = W S V'`, its polar
   part `U = W V'` defines the covariant transport and `U[1,1]` is the target
   link.  This separates the intrinsic target geometry from nonunitary
   distortion caused by a finite hybrid-Wannier frame.
4. Wilson plaquettes of the target links give the Berry curvature.  Projector
   distances are evaluated along both reciprocal steps and their difference.
   The three distances determine the full metric tensor after transforming the
   nonorthogonal reciprocal coordinates to Cartesian `x,y` coordinates.  The
   reported trace is therefore the physical

       Tr g(k) = g_xx(k) + g_yy(k),

   not the sum of diagonal components in the oblique lattice coordinates.  For
   the C3-symmetric tMBG case, the locally saturating Kahler component is

       g_ideal(k) = |Omega(k)| I / 2.

   The intrinsic metric is written as
   `g_intrinsic = g_ideal + g_residual`.  The saved integrated diagnostics are

       eta = integral Tr(g_residual),
       ideal_fraction = integral |Omega| / integral Tr(g_intrinsic).

5. To obtain the momentum-resolved ideal component, the isolated magnetic
   state is lifted into the same six-local-orbital plane-wave Hilbert space as
   the zero-field eigenstates.  A unitary real-space phase converts the
   hybrid-Wannier Landau-gauge convention to symmetric gauge, and the result is
   projected onto the true zero-field energy eigenbasis.  The scalar target
   weight is finally projected onto its C3-invariant representation by
   averaging each three-point momentum orbit.  This operation preserves the
   integrated projection weight exactly.  Both the raw and C3-restored weights
   are exported for auditing.

The six local orbitals are the physical three-layer times two-sublattice
degrees of freedom of tMBG.  They are distinct from the six retained
zero-field bands.  Plane-wave convergence is controlled by the reciprocal
grid `lg`; the validated calculation uses `lg=11`, or `6*11^2=726` local
plane-wave orbitals at each momentum.  Magnetic-basis convergence is
controlled separately by `smax=5`.

## Code flow used by the calculation

There is only one numerical entry point:
`examples/run_tmbg_stability_check.jl`.  It calls the source files in the
following order; these are stages of one calculation, not alternative
methods.

| Stage | Source code and main call | Result |
| --- | --- | --- |
| Continuum model | `TMBGZeroField.jl`, `WangContinuum.jl`; `build_tmbg_wang` | tMBG Hamiltonian, lattice, plane-wave basis, and field coupling |
| Hybrid-Wannier basis | `TMBGHybridWannier.jl`; called inside `build_magnetic_hw` | Six-band zero-field hybrid-Wannier basis |
| Magnetic spectrum | `TMBGMagneticHW.jl`; `build_magnetic_hw` | Nonorthogonal magnetic problem, orthonormalized subbands, and `magnetic_spectrum.csv` |
| Multiband quantum geometry | `TMBGCommonBasis.jl`, `TMBGMagneticGeometry.jl`, `TMBGIntrinsicIdealGeometry.jl`, and `TMBGIdealComponent.jl`; `compute_ideal_component` | Berry curvature, Cartesian quantum metric, trace condition, and overlap diagnostics |
| Momentum-resolved ideal component | `TMBGProjection.jl`, `TMBGSymmetricGaugeProjection.jl`, and `TMBGIdealComponentProjection.jl`; `write_ideal_component_projection` | Common-basis zero-field projection and C3-restored ideal-component distribution |

`src/TMBGMagneticHybridWannier.jl` only assembles these modules into the Julia
package.  The two Python files do not perform any physics calculation:
`plot_magnetic_spectrum.py` reads `magnetic_spectrum.csv`, while
`plot_ideal_component_projection.py` reads the exported projection tables.
Thus a complete numerical rerun requires only the Julia entry point; plotting
is a separate final step.

## Reproduce the validated tMBG case

Julia 1.10 or newer is recommended.

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. examples/run_tmbg_stability_check.jl
```

The single entry point fixes the complete validated parameter set:

- twist angle `theta=1.04 deg`;
- tunneling ratio `kappa=0.70`;
- displacement potential `U=50 meV`;
- monolayer hBN sublattice potential `-30 meV`;
- physical magnetic flux `p/q=1/20`;
- six retained zero-field bands;
- reciprocal plane-wave grid `lg=11`;
- magnetic-translation cutoff `smax=5`.

The output directory is `results/tmbg_pq1_20_lg11_smax5/`.  A different
directory may be supplied as the only command-line argument.  The program
writes, in order:

- `case_parameters.txt`: all physical and numerical parameters;
- `magnetic_hw_cache.jls`: reusable full magnetic model cache;
- `magnetic_spectrum.csv`: complete magnetic-subband energies;
- `trace_condition_summary.txt`: Chern numbers, integrated trace condition,
  and ideal fraction;
- `quantum_geometry.csv`: local Berry curvature and ideal, residual, and
  intrinsic Cartesian quantum metrics;
- `overlap_diagnostics.csv` and `overlap_diagnostics_summary.txt`: raw and
  covariant target links and multiband-frame singular-value distributions;
- `ideal_component_projection/momentum_resolved_target_weight.csv`: raw and
  C3-restored momentum-resolved target weights;
- `ideal_component_projection/zero_field_band_weights.csv`: integrated
  projection onto each retained zero-field band;
- `ideal_component_projection/zero_field_path_*.csv` and metadata: zero-field
  path spectrum used by the projection visualization;
- `ideal_component_projection/projection_summary.txt`: integrated target
  fraction and C3 diagnostics.

The cache is large and is deliberately excluded from Git.  The lightweight
numerical tables and plots of the validated run are versioned under
`reference/tmbg_pq1_20_lg11_smax5/`.

### Saved trace-condition result

The validated example always saves the full local quantum-geometry table and
the integrated trace-condition summary.  The final `lg=11, smax=5` values are
listed here after the reference calculation so that a rerun can be checked
without interpreting a plot.

<!-- TRACE_RESULT_BEGIN -->
For `p/q=1/20`, `lg=11`, and `smax=5`, the saved result is:

- isolated target magnetic band `79`;
- multiband correction frame `(79, 78)`;
- target Chern number `C=1` and frame Chern number `C_frame=2`;
- ideal integrated trace `6.283185307179585`;
- residual integrated trace and trace excess
  `eta=0.00022470507632712814`;
- total integrated trace `6.283410012255914`;
- integrated ideal fraction `0.9999642383553055`;
- minimum correction-frame link singular value `0.8379962538093386`.

The momentum-resolved projection independently places
`0.919390352703438` of the central zero-field pair in the target band.  Its
raw C3 relative error is `0.12447892706843225`; the explicit C3 group
projection reduces it to `7.209450195207299e-17` while preserving all
integrated weights.

The exact values are stored in
`reference/tmbg_pq1_20_lg11_smax5/trace_condition_summary.txt` and
`reference/tmbg_pq1_20_lg11_smax5/ideal_component_projection/projection_summary.txt`.
<!-- TRACE_RESULT_END -->

## Plot the two physical observables

Install the lightweight plotting dependencies once:

```bash
python -m pip install -r requirements-plot.txt
```

Plot the magnetic hybrid-Wannier spectrum in the target-band energy window:

```bash
python examples/plot_magnetic_spectrum.py
```

The default `15--40 meV` window reproduces the saved reference plot.  The
flags `--energy-min` and `--energy-max` select another window, while
`--full-range` shows every energy in the complete exported spectrum.

Plot the momentum-resolved ideal-component projection:

```bash
python examples/plot_ideal_component_projection.py
```

Each program writes both PDF and PNG.  Berry curvature, quantum metric, trace
condition, and overlap distributions remain available as numerical CSV and
summary outputs; no additional plotting entry points are included in the
minimal workflow.

## Tests

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

The tests cover lattice-coordinate conventions, continuum-basis equivalence,
magnetic overlap and Hamiltonian hermiticity, Cartesian metric reconstruction,
multiband polar-frame covariance, and C3 group projection.

## Method references

- Xiaoyu Wang,
  [`Hofstadter-TBG`](https://github.com/xywang2017/Hofstadter-TBG), the
  upstream GPL-3.0 hybrid-Wannier implementation.
- Xiaoyu Wang and Oskar Vafek, *Narrow bands in magnetic field and
  strong-coupling Hofstadter spectra*, Phys. Rev. B **106**, L121111 (2022).
- Xiaoyu Wang and Oskar Vafek, *Revisiting Bloch electrons in a magnetic
  field: Hofstadter physics via hybrid Wannier states*, Phys. Rev. B **108**,
  245109 (2023).
- Patrick J. Ledwith, Ashvin Vishwanath, and Daniel E. Parker,
  *Vortexability: A Unifying Criterion for Ideal Fractional Chern
  Insulators*, Phys. Rev. B **108**, 205144 (2023).
