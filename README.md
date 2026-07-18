# tMBG magnetic hybrid-Wannier spectrum and geometry

This Julia package implements a self-contained hybrid-Wannier calculation of
the magnetic spectrum and quantum geometry of twisted monolayer-bilayer
graphene (tMBG).  The paper example constructs the magnetic Hilbert space from
zero-field continuum-model bands, solves the generalized magnetic eigenproblem,
and extracts the ideal component of an isolated magnetic subband through a
multiband covariant frame.

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
   zero-field energy eigenbasis.  The exported `target_relative` field is the
   target-band weight normalized by its maximum over momentum; the integrated
   target share within the central zero-field pair is reported separately.

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

For a short installation check, use:

```bash
julia --project=. examples/paper_example.jl --quick
```

To plot the paper geometry after the full calculation:

```bash
python -m pip install -r requirements-plot.txt
python examples/plot_paper_geometry.py
python examples/plot_overlap_distribution.py
python examples/plot_fig4_ideal_component_projection.py
```

The lightweight tables and final plots from the converged reference run are
versioned under `reference/paper_q20_converged/`.  The approximately 500 MB
serialized magnetic cache is regenerated locally and is not committed.

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

- Xiaoyu Wang and Oskar Vafek, *Narrow bands in magnetic field and
  strong-coupling Hofstadter spectra*, Phys. Rev. B **106**, L121111 (2022).
- Xiaoyu Wang and Oskar Vafek, *Revisiting Bloch electrons in a magnetic
  field: Hofstadter physics via hybrid Wannier states*, Phys. Rev. B **108**,
  245109 (2023).
- Patrick J. Ledwith, Ashvin Vishwanath, and Daniel E. Parker,
  *Vortexability: A Unifying Criterion for Ideal Fractional Chern
  Insulators*, Phys. Rev. B **108**, 205144 (2023).
