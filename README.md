# Hyperlinear Check for FCI Non-Ideal Flat-Band Stability

**Hyperlinear FCI Stability Check** is a Julia framework for diagnosing FCI
stability in non-ideal flat bands through magnetic-subband spectroscopy,
multiband quantum geometry, trace-condition analysis, overlap diagnostics, and
ideal-component projection.  The currently validated implementation targets
twisted monolayer-bilayer graphene (tMBG).  Its paper example constructs the
magnetic Hilbert space from zero-field continuum-model bands, solves the
generalized magnetic eigenproblem, and extracts the ideal component of an
isolated magnetic subband through a multiband covariant frame.

## Origin and scope

This project is based on and extends Xiaoyu Wang's
[`Hofstadter-TBG`](https://github.com/xywang2017/Hofstadter-TBG) project, which
implements the Wang--Vafek hybrid-Wannier magnetic-translation construction
for twisted bilayer graphene.  The present package adapts that construction to
the tMBG continuum Hamiltonian and adds the nonorthogonal magnetic solve,
complete field-coupling matrix elements, common-basis symmetric-gauge
projection, multiband quantum geometry, overlap diagnostics, and the Fig. 4
ideal-component observable used in the accompanying paper.  Please cite the
original Wang--Vafek works listed below when using this implementation.

The code is organized so that the same workflow can be extended to other
continuum-model moire systems.  A new model backend must provide its continuum
Hamiltonian, reciprocal and direct lattices, momentum-linear field coupling,
local orbital dimension, and reciprocal-space sewing maps.  The resulting
magnetic subbands can then be tested through their spectrum, Chern number,
Berry curvature, quantum metric and trace condition, correction-frame
overlaps, and ideal-component distribution.  These are the single-particle
inputs used to assess candidate magnetic-subband FCI stability and can be
passed to a separate projected-interaction or exact-diagonalization workflow.
At present, the fully implemented and numerically validated paper example is
restricted to tMBG; other moire continuum models require a model-specific
backend and independent validation.

## Physical construction

The calculation has five steps.

1. The tMBG continuum Hamiltonian is diagonalized on a uniform zero-field
   Brillouin-zone mesh.  Non-Abelian parallel transport along `k1` produces
   hybrid-Wannier states, which are localized along the first direct-lattice
   direction and Bloch extended along the second.
2. Magnetic translations of these states form a nonorthogonal basis at flux
   `p/q` per moire unit cell.  The overlap matrix `O(k)` and magnetic
   Hamiltonian `H(k)` are evaluated in the same basis.  Retaining the physical
   rank and applying `O^(-1/2)` gives the orthonormal magnetic Hamiltonian.
3. The target magnetic subband is placed first in a local multiband correction
   frame `Phi(k)`.  For neighboring momenta, the common-basis overlap is

       M(k,k+dk) = Phi(k)' T(dk) Phi(k+dk).

   Its polar part `U = W V'`, obtained from `M = W S V'`, is a covariant frame
   transport.  The target link is `U[1,1]`; nearby bands define the connection
   but are not merged into the target projector.
4. Wilson plaquettes of the target links give the Berry curvature.  Projector
   distances along both reciprocal directions and their difference determine
   the metric in that oblique basis.  The code transforms it to Cartesian
   coordinates before evaluating `Tr g`.  For the C3-symmetric paper model,
   the locally saturating Kahler component is

       g_ideal(k) = |Omega(k)| I / 2.

   The residual metric is reconstructed from the nonnegative covariant
   projector distances left by the polar frame.  The intrinsic result is
   `g_intrinsic = g_ideal + g_residual`, with
   `eta = integral Tr(g_residual)` and
   `ideal_fraction = integral |Omega| / integral Tr(g_intrinsic)`.
5. For the Fig. 4 observable, the isolated weak-field branch is lifted into
   the same six-orbital plane-wave Hilbert space as the zero-field states.  A
   single exact gauge transformation converts the hWF Landau-gauge convention
   to symmetric gauge, after which the state is projected onto the true
   zero-field energy eigenbasis.  Since the continuum model and magnetic field
   preserve C3, the final scalar weight is projected onto its C3-invariant
   representation by averaging each three-point momentum orbit.  This keeps
   the total projection weight exactly unchanged while removing truncation-
   induced non-C3 components.  The exported `target_relative` field is this
   symmetry-restored target-band weight normalized by its maximum; the raw
   weight is retained as `target_weight_raw`, and the integrated target share
   within the central zero-field pair is reported separately.

The implementation uses the six-orbital plane-wave continuum basis as its
common Hilbert space.  No external basis-expansion program or reference
wavefunction is read by the spectrum or geometry calculation.

## Reproduce the paper example

Julia 1.10 or newer is recommended.

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. examples/paper_example.jl
```

The default parameters are the paper values
`theta=1.04 deg`, `kappa=0.70`, `U=50 meV`, hBN potential `-30 meV`, and
physical flux `p/q=1/20`.  The reference run uses a six-band hWF subspace,
reciprocal cutoff `lg=9`, and magnetic-translation cutoff `smax=5`; the latter
two replace the less converged `lg=7, smax=3` plotting cache.  The script
writes:

- `results/paper_q20_converged/magnetic_hw.jls`: reusable magnetic hWF cache;
- `results/paper_q20_converged/magnetic_spectrum.csv`: the complete magnetic
  spectrum;
- `results/paper_q20_converged/ideal_component_geometry.csv`: Berry curvature,
  ideal, residual, and total Cartesian metric;
- `results/paper_q20_converged/ideal_component_summary.txt`: Chern number,
  trace excess, and ideal fraction;
- `results/paper_q20_converged/overlap_distribution.csv`: all raw/covariant
  target-link moduli and multiband-frame singular values;
- `results/paper_q20_converged/overlap_distribution_summary.txt`: quantiles
  and moments of each overlap diagnostic;
- `results/paper_q20_converged/fig4/upper_component_grid.csv`: the Fig. 4
  momentum-resolved projection of the separated magnetic component onto the
  zero-field target band;
- `results/paper_q20_converged/fig4/ideal_component_projection.pdf`: the
  complete Fig. 4 path-and-hexagon visualization.

### Saved paper-case results

The paper example intentionally saves both the Fig. 4 distribution and the
quantum-geometric trace-condition result.  For the converged `p/q=1/20`
calculation included under `reference/paper_q20_converged/`, the isolated
magnetic target band has

- target magnetic band index `79` and correction frame `(79, 78)`;
- target-band Chern number `C=1` and two-band frame Chern number `C_frame=2`;
- ideal integrated trace `6.283185307179586 = 2*pi`;
- residual integrated trace, equivalently the trace-condition excess,
  `eta=2.2469106188858858e-4`;
- total integrated trace `6.283409998241475`;
- integrated ideal fraction `0.9999642405856133`;
- minimum correction-frame link singular value `0.8380051901128523`.

The corresponding Fig. 4 projection retains `0.919354982036541` of the
weight in the target band within the central zero-field pair.  Its
C3-averaged K-region ideal component is `0.12827514648901675` relative to the
symmetry-restored map maximum.  Before C3 group projection, the raw relative
C3 error is `0.12450748406769918`; after group projection it is
`6.446235956802648e-17`.  The relative norm of the removed non-C3 component is
`0.07188442944260906`, while the integrated target weight and the `91.9355%`
central-pair fraction are unchanged.  The trace-condition values are stored in
`reference/paper_q20_converged/ideal_component_summary.txt`; the local Berry
curvature and Cartesian metric data are stored in
`reference/paper_q20_converged/ideal_component_geometry.csv`.  These numerical
outputs are part of the reproducible example even though the repository only
ships plotting code for Fig. 4.

The common-plane-wave convergence was checked separately at fixed magnetic
flux and fixed momentum mesh.  Increasing the reciprocal grid from `lg=9`
(`6*81=486` plane-wave orbitals) to `lg=11` (`6*121=726`) at `smax=3`
changes the raw C3 error only from `0.12345514116635437` to
`0.12342601054800617`, and changes the central-pair weight only from
`0.9193434379517239` to `0.9193788367118987`.  The common-basis norm remains
unity.  Thus `lg=9` is already converged for this observable; the residual
non-C3 contribution is removed by the explicit physical C3 group projection,
not by momentum interpolation or plot smoothing.  The six local orbitals are
the fixed three-layer times two-sublattice degrees of freedom of tMBG; basis
convergence is controlled by the reciprocal `G` grid and magnetic-translation
cutoffs, not by changing this physical local-orbital count.

For a short installation check, use:

```bash
julia --project=. examples/paper_example.jl --quick
```

The repository provides one plotting entry point, reproducing the Fig. 4
path-and-hexagon ideal-component distribution:

```bash
python -m pip install -r requirements-plot.txt
python examples/plot_fig4_ideal_component_projection.py
```

The quantum-geometry and overlap calculations remain available as numerical
CSV and summary outputs, but no separate plotting example is shipped for
them.  The lightweight tables and the final Fig. 4 plot from the converged
reference run are versioned under `reference/paper_q20_converged/`.  The
approximately 500 MB serialized magnetic cache is regenerated locally and is
not committed.

## Tests

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

The tests cover reciprocal/direct lattice conventions, global central-band
indices, velocity matrix elements, magnetic overlap/Hamiltonian hermiticity,
the nonorthogonal-to-Cartesian metric transformation, polar-frame covariance,
and the C3 ideal trace condition.

## Main source files

- `src/TMBGZeroField.jl`: tMBG continuum Hamiltonian.
- `src/WangContinuum.jl`: reciprocal-frame continuum model and magnetic
  coupling conventions.
- `src/TMBGHybridWannier.jl`: Wilson-loop hybrid-Wannier construction.
- `src/TMBGMagneticHW.jl`: nonorthogonal magnetic basis and spectrum.
- `src/TMBGMagneticGeometry.jl`: common-basis lifting and Cartesian metric.
- `src/TMBGIdealComponent.jl`: multiband polar frame and ideal-component
  extraction.
- `src/TMBGSymmetricGaugeProjection.jl`: common-basis gauge transform and
  zero-field band projection.
- `src/TMBGFigure4Projection.jl`: complete Fig. 4 data export, including the
  path spectrum and hexagonal momentum distribution.

## Method references

- Xiaoyu Wang, [`Hofstadter-TBG`](https://github.com/xywang2017/Hofstadter-TBG),
  the upstream GPL-3.0 hybrid-Wannier implementation on which this project is
  based.
- Xiaoyu Wang and Oskar Vafek, *Narrow bands in magnetic field and
  strong-coupling Hofstadter spectra*, Phys. Rev. B **106**, L121111 (2022).
- Xiaoyu Wang and Oskar Vafek, *Revisiting Bloch electrons in a magnetic
  field: Hofstadter physics via hybrid Wannier states*, Phys. Rev. B **108**,
  245109 (2023).
- Patrick J. Ledwith, Ashvin Vishwanath, and Daniel E. Parker,
  *Vortexability: A Unifying Criterion for Ideal Fractional Chern
  Insulators*, Phys. Rev. B **108**, 205144 (2023).

## License

Because this project is derived from `Hofstadter-TBG`, it is distributed under
the GNU General Public License v3.0.  See `LICENSE` for the full terms.
