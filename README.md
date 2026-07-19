# Non-Ideal Moiré FCI Stability

## Purpose and physical background

Establishing a fractional Chern insulator (FCI) normally requires a
strongly-correlated many-body calculation, such as exact diagonalization or
DMRG, because the FCI is an interacting many-electron phase.  This repository
provides a complementary diagnostic at the single-electron continuum-model
level.  It asks whether a non-ideal zero-field flat band contains a latent
near-ideal Chern component that supplies favorable single-particle geometry
for an FCI.  The many-body calculation establishes the phase itself; the code
here detects, characterizes, and locates its candidate ideal component.

### Why the non-ideal problem is essential

The ideal trace condition is a constructive reference point, not a realistic
material requirement.  It makes the Landau-level mapping and vortex
attachment exact, but continuum models of moiré materials also contain
lattice relaxation, remote hopping, displacement fields, and substrate
potentials.  These ingredients generically produce Berry-curvature
inhomogeneity and finite trace-condition violation.  A diagnostic restricted
to globally ideal bands would therefore omit much of the experimentally
accessible parameter space before asking whether correlations can stabilize
an FCI.

The associated tMBG article gives a concrete example.  As the relaxation
parameter crosses $\kappa_c\simeq0.55$, the non-interacting parent band remains
$C=2$, while its global quantum geometry degrades and its hybrid-Wannier flow
develops a local geometric instability.  Exact diagonalization nevertheless
finds a transition from a Halperin-(112)-like phase to a Laughlin-like
$1/3$ FCI.  The latter persists on the strongly non-ideal side, so its origin
is not captured by assigning one global ideality score to the entire parent
band.

The relevant question is instead component resolved.  A non-ideal Bloch band
can contain sectors with sharply different quantum geometry.  In the tMBG
regime studied here, interactions favor a hidden near-ideal $C=1$ component
embedded in the non-ideal $C=2$ zero-field states, while its strongly
non-ideal partner contributes much less to the Laughlin-like state.  The weak
magnetic field does not create this favorable component: it resolves the
pre-existing components into separate magnetic subbands, after which the
near-ideal branch can be characterized and projected back onto the original
zero-field band.  This extends ideal-band design principles into a practical
diagnostic for realistic non-ideal moiré bands.

The central procedure proposed in the associated article is:

1. Construct the zero-field non-ideal band directly from the continuum
   Hamiltonian.
2. Add a weak perpendicular magnetic field.  The field mixes nearby momenta
   and resolves components embedded in the same zero-field Bloch band into
   distinct magnetic subbands.  This weak-field resolution is the
   single-particle **color-separation probe**.
3. Select the separated subband with near-ideal quantum geometry and compute
   its Chern number, Berry curvature, Quantum Metric, trace condition, and
   multiband overlap diagnostics.
4. Lift that magnetic subband into a common plane-wave Hilbert space, convert
   it to symmetric gauge, and project it onto the true zero-field energy
   eigenbasis.  This produces the momentum-resolved distribution of the ideal
   component inside the original zero-field non-ideal band.

Thus the magnetic field is a controlled resolution tool: the object being
tested is the original non-ideal band, and the final output identifies where
its potential ideal component is distributed at zero field.  The implemented
and validated example is twisted monolayer-bilayer graphene (tMBG) at
`p/q=1/20`.  No Landau-level expansion, fitted reference wavefunction, or
external geometry data are used.

## What "ideal" means

For a rank-one Chern band, the quantum geometric tensor is

~~~math
Q_{\mu\nu}(\mathbf{k})=
\langle \partial_{k_\mu}u_{\mathbf{k}}|
(1-|u_{\mathbf{k}}\rangle\langle u_{\mathbf{k}}|)
|\partial_{k_\nu}u_{\mathbf{k}}\rangle .
~~~

Its real part is the Quantum Metric
$g_{\mu\nu}=\mathrm{Re}Q_{\mu\nu}$, and its antisymmetric imaginary
part gives the Berry curvature $\Omega$.  They obey the pointwise geometric
bound

~~~math
\mathrm{Tr}g(\mathbf{k})\geq |\Omega(\mathbf{k})|.
~~~

An ideal Chern band saturates this trace condition throughout the Brillouin
zone.  Under C3 symmetry, the locally saturating component used here is
$g_{\mathrm{ideal}}=|\Omega|I/2$.  The code measures the departure from this
limit through the residual trace

~~~math
\eta=\int_{\mathrm{BZ}}
[\mathrm{Tr}g(\mathbf{k})-|\Omega(\mathbf{k})|]\,d^2k,
~~~

and reports an integrated ideal fraction, which approaches one as the
residual contribution vanishes.  Berry-curvature uniformity is a useful
additional diagnostic, but pointwise trace saturation is the ideality
criterion used by this workflow.

This geometry has a constructive meaning.  The exact Landau-level mapping of
ideal flat bands and the vortexability criterion show when vortex attachment
can be carried out within the band subspace, providing the single-particle
structure needed to build FCI wavefunctions.  This is the physical reason the
repository uses ideal-component geometry as its FCI-stability diagnostic.

## Associated article

This code accompanies and was developed for the numerical analysis in the
arXiv manuscript:

> Moru Song and Kai Chang,
> *[Fractional Chern Insulators Transition in Non-ideal Flat Bands of Twisted
> Mono-bilayer Graphene](https://arxiv.org/abs/2511.12231)*,
> arXiv:2511.12231 [cond-mat.mes-hall] (2025).

In particular, the repository provides the independent hybrid-Wannier
calculation of the weak-field magnetic spectrum, the multiband quantum
geometry and trace condition of the isolated magnetic subband, and the
momentum-resolved ideal-component distribution discussed in that work.

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
   correction frame $\Phi(\mathbf{k})$.  Nearby subbands stabilize the frame
   but are not merged into the target projector.  For neighboring momenta,
   the link matrix $M$ is constructed from
   $\Phi_{\mathbf{k}}^\dagger T_{\delta\mathbf{k}}\Phi_{\mathbf{k}+\delta\mathbf{k}}$
   and lifted to a common plane-wave Hilbert space.  Here
   $\Phi_{\mathbf{k}}\equiv\Phi(\mathbf{k})$ and
   $T_{\delta\mathbf{k}}\equiv T(\delta\mathbf{k})$.  If
   $M=WSV^\dagger$, its polar part $U=WV^\dagger$ defines the covariant
   transport, and $U_{11}$ is the target link.  This separates the intrinsic
   target geometry from nonunitary distortion caused by a finite
   hybrid-Wannier frame.
4. Wilson plaquettes of the target links give the Berry curvature.  Projector
   distances are evaluated along both reciprocal steps and their difference.
   The three distances determine the full quantum metric after accounting for
   the angle between the two reciprocal-lattice directions.  The reported
   trace is therefore the physical quantity
   $\mathrm{Tr}g(\mathbf{k})=g_{xx}(\mathbf{k})+g_{yy}(\mathbf{k})$,
   not the sum of diagonal components in the oblique lattice coordinates.
   For the C3-symmetric tMBG case, the locally saturating ideal component is
   $g_{\mathrm{ideal}}(\mathbf{k})=|\Omega(\mathbf{k})|I/2$.  The intrinsic
   metric is written as
   $g_{\mathrm{intrinsic}}=g_{\mathrm{ideal}}+g_{\mathrm{residual}}$.
   The saved integrated diagnostics are the residual trace
   $\eta=\int_{\mathrm{BZ}}\mathrm{Tr}g_{\mathrm{residual}}d^2k$
   and the ideal fraction $f_{\mathrm{ideal}}=I_{\Omega}/I_g$, with
   $I_{\Omega}=\int_{\mathrm{BZ}}|\Omega|d^2k$ and
   $I_g=\int_{\mathrm{BZ}}\mathrm{Tr}g_{\mathrm{intrinsic}}d^2k$.

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
`examples/run_tmbg_stability_check.jl`.  It executes the following five stages
in order.  They are parts of one calculation, not alternative methods.

1. **Continuum model.** `TMBGZeroField.jl` and `WangContinuum.jl` call
   `build_tmbg_wang` to construct the tMBG Hamiltonian, lattice, plane-wave
   basis, and field coupling.

2. **Hybrid-Wannier basis.** `TMBGHybridWannier.jl` is called inside
   `build_magnetic_hw` to construct the six-band zero-field hybrid-Wannier
   basis.

3. **Magnetic spectrum.** `TMBGMagneticHW.jl` calls `build_magnetic_hw` to
   solve the nonorthogonal magnetic problem and export the orthonormalized
   subbands to `magnetic_spectrum.csv`.

4. **Multiband quantum geometry.** `TMBGCommonBasis.jl`,
   `TMBGMagneticGeometry.jl`, `TMBGIntrinsicIdealGeometry.jl`, and
   `TMBGIdealComponent.jl` call `compute_ideal_component` to obtain the Berry
   curvature, Quantum Metric, trace condition, and overlap diagnostics.

5. **Momentum-resolved ideal component.** `TMBGProjection.jl`,
   `TMBGSymmetricGaugeProjection.jl`, and
   `TMBGIdealComponentProjection.jl` call
   `write_ideal_component_projection` to produce the common-basis zero-field
   projection and the C3-restored ideal-component distribution.

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
  intrinsic quantum metrics;
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
magnetic overlap and Hamiltonian hermiticity, quantum-metric reconstruction
from nonorthogonal reciprocal coordinates, multiband polar-frame covariance,
and C3 group projection.

## Method references

- Moru Song and Kai Chang,
  *[Fractional Chern Insulators Transition in Non-ideal Flat Bands of Twisted
  Mono-bilayer Graphene](https://arxiv.org/abs/2511.12231)*,
  arXiv:2511.12231 [cond-mat.mes-hall] (2025), the article accompanied by this
  repository.
- Jie Wang, Jennifer Cano, Andrew J. Millis, Zhao Liu, and Bo Yang,
  *[Exact Landau Level Description of Geometry and Interaction in a
  Flatband](https://doi.org/10.1103/PhysRevLett.127.246403)*, Phys. Rev. Lett.
  **127**, 246403 (2021), the exact ideal-flat-band/Landau-level mapping.
- Patrick J. Ledwith, Ashvin Vishwanath, and Daniel E. Parker,
  *[Vortexability: A Unifying Criterion for Ideal Fractional Chern
  Insulators](https://doi.org/10.1103/PhysRevB.108.205144)*, Phys. Rev. B
  **108**, 205144 (2023), the constructive ideal-band criterion underlying
  the stability diagnostic.
- Zhen Ma, Shuai Li, Ya-Wen Zheng, Meng-Meng Xiao, Hua Jiang, Jin-Hua Gao,
  and X. C. Xie, *[Topological flat bands in twisted trilayer
  graphene](https://doi.org/10.1016/j.scib.2020.10.004)*, Science Bulletin
  **66**, 18-22 (2021), the six-orbital tMBG continuum-Hamiltonian and
  realistic bilayer-hopping reference used in `src/TMBGZeroField.jl`.
- Xiaoyu Wang,
  [`Hofstadter-TBG`](https://github.com/xywang2017/Hofstadter-TBG), the
  upstream GPL-3.0 hybrid-Wannier implementation.
- Xiaoyu Wang and Oskar Vafek, *Narrow bands in magnetic field and
  strong-coupling Hofstadter spectra*, Phys. Rev. B **106**, L121111 (2022).
- Xiaoyu Wang and Oskar Vafek, *Revisiting Bloch electrons in a magnetic
  field: Hofstadter physics via hybrid Wannier states*, Phys. Rev. B **108**,
  245109 (2023).
