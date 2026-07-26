# Non-Ideal Moiré FCI Stability

This repository provides a single-particle screening workflow for a question
that normally requires a much more expensive many-body calculation: does a
realistic, non-ideal moiré band contain a geometrically favorable component
that can support a fractional Chern insulator (FCI)?

The calculation uses a weak perpendicular magnetic field as a resolution
probe.  The field separates components hidden inside the same zero-field band;
the code identifies the near-ideal magnetic subband and projects it back onto
the original zero-field eigenstates.  The main result is therefore not merely
a magnetic spectrum, but the **momentum-resolved distribution of the candidate
ideal component in the non-ideal zero-field band**.

This is an FCI-stability diagnostic, not a replacement for exact
diagonalization, DMRG, or another interacting calculation.  It identifies
favorable single-particle structure and supplies a controlled target for the
many-body analysis.

## Why non-ideal bands need a component-resolved test

Exact ideal geometry is a useful reference, but realistic continuum models
include relaxation, remote hopping, displacement fields, and substrate
potentials.  These effects make the Berry curvature and Quantum Metric
nonuniform, so testing only globally ideal bands would exclude much of the
physically relevant parameter space.

A single global geometry score can also hide internal structure.  The tMBG
case studied here contains a near-ideal component embedded in a non-ideal
parent band.  The weak-field calculation resolves that component without
assuming that the full zero-field band is ideal, then locates it again at zero
field.

## Choose the calculation you need

| Goal | Entry point | Role |
| --- | --- | --- |
| Identify and map the ideal component of a non-ideal band | `examples/run_tmbg_stability_check.jl` | Primary, independent hybrid-Wannier workflow |
| Export the conventional LL magnetic spectrum | `reference_methods/landau_level/compute_ll_magnetic_spectrum.py` | Standard-method comparison reference |
| Plot an exported LL magnetic spectrum | `reference_methods/landau_level/plot_ll_magnetic_spectrum.py` | Standard-method comparison reference |
| Compute LL Berry curvature, Quantum Metric, and trace condition | `reference_methods/landau_level/compute_ll_quantum_geometry.py` | Standard-method comparison reference |

The hybrid-Wannier calculation does not read LL eigenstates, LL geometry, or a
fitted LL wavefunction.  The LL implementation is included as a conventional
Landau-level-expansion benchmark for methodological comparison.

## Quick start: primary hybrid-Wannier workflow

Julia 1.10 or newer is recommended.

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. examples/run_tmbg_stability_check.jl
```

The second command executes the complete validated tMBG workflow in one
sequence:

1. construct the zero-field tMBG continuum bands;
2. build the magnetic hybrid-Wannier basis and export the full magnetic
   spectrum;
3. select the isolated target subband and compute its multiband-covariant
   geometry;
4. lift the target state to the common plane-wave basis;
5. project it onto the zero-field bands and save the ideal-component
   distribution.

Results are written to `results/tmbg_pq1_20_lg11_smax5/` by default.  The
complete magnetic model cache is large and is intentionally regenerated
locally.

To plot the two main observables:

```bash
python -m pip install -r requirements-plot.txt
python examples/plot_magnetic_spectrum.py
python examples/plot_ideal_component_projection.py
```

## Quick start: conventional LL comparison

Create a separate Python environment and install:

```bash
python -m pip install -r requirements-ll.txt
```

First run the small installation checks:

```bash
python reference_methods/landau_level/compute_ll_magnetic_spectrum.py --quick
python reference_methods/landau_level/compute_ll_quantum_geometry.py --quick
```

These quick calculations are deliberately unconverged.  The paper-scale
geometry reference is generated with:

```bash
python reference_methods/landau_level/compute_ll_quantum_geometry.py \
  --nll 400 --nk1 40 --nk2 16 --band-start 82 \
  --geometry-backend block --jobs 8 \
  --output results/ll_reference/quantum_geometry.npz

python reference_methods/landau_level/plot_ll_quantum_geometry.py \
  --input results/ll_reference/quantum_geometry.npz
```

When `--jobs` is larger than one, set the BLAS thread count to one to avoid
oversubscription.  The versioned converged output is in
`reference/tmbg_pq1_20_ll_nll400_40x16/`.

## What is saved

The primary workflow saves:

- the complete magnetic-subband spectrum;
- local Berry curvature, Quantum Metric, and trace-condition data;
- multiband overlap diagnostics;
- raw and C3-restored momentum-resolved ideal-component weights;
- integrated weights on each retained zero-field band;
- a zero-field path spectrum for the projection plot.

The LL comparison saves its magnetic spectrum separately and stores Berry
curvature, metric trace, trace excess, determinant condition, Chern number,
and integrated trace-condition measures in one portable `.npz` file.

## Validated tMBG example

The repository fixes the paper example at twist angle `1.04 deg`,
`kappa=0.70`, displacement potential `50 meV`, hBN sublattice potential
`-30 meV`, and physical flux `p/q=1/20`.

For the primary hybrid-Wannier calculation:

- six zero-field bands are retained;
- the plane-wave cutoff is `lg=11` and the magnetic-translation cutoff is
  `smax=5`;
- target magnetic band `79` has `C=1`;
- the integrated ideal fraction is `0.9999642383553055`;
- the target component occupies `0.919390352703438` of the central
  zero-field pair.

For the conventional LL reference, target band `82` at `NLL=400` on a
`40 x 16` mesh has `|C|=1` and integrated trace excess
`eta=0.098808121157`.  The Chern sign differs only by the orientation
convention used in the LL implementation.

## Repository guide

Readable entry points:

- `examples/run_tmbg_stability_check.jl` — complete primary calculation;
- `examples/plot_magnetic_spectrum.py` — hybrid-Wannier spectrum plot;
- `examples/plot_ideal_component_projection.py` — zero-field ideal-component
  distribution;
- `reference_methods/landau_level/README.md` — LL comparison instructions.

Technical material:

- `src/` — hybrid-Wannier continuum, magnetic, geometry, and projection
  modules;
- `reference_methods/landau_level/` — conventional LL Hamiltonian, spectrum,
  geometry, and plotting code;
- `docs/technical_method.md` — basis conventions, geometry definitions, and
  numerical details;
- `reference/` — lightweight outputs of the validated calculations.

## Scope, origin, and license

The implemented and validated material backend is twisted monolayer-bilayer
graphene (tMBG).  The architecture can be extended to other continuum-model
moiré systems by supplying their Hamiltonian, lattice conventions, local
orbital basis, field coupling, and reciprocal-space sewing maps.

The hybrid-Wannier implementation is based on Xiaoyu Wang's
[`Hofstadter-TBG`](https://github.com/xywang2017/Hofstadter-TBG) project and
the Wang--Vafek magnetic-translation construction.  This repository adds the
tMBG continuum model, generalized nonorthogonal solve, common-basis lifting,
multiband geometry, overlap diagnostics, and zero-field component projection.
The conventional LL implementation is supplied as an independent comparison
reference.

The code is distributed under GNU GPL v3.0; see `LICENSE`.

## Associated article and method references

- M. Song and K. Chang,
  *[Fractional Chern Insulators Transition in Non-ideal Flat Bands of Twisted
  Mono-bilayer Graphene](https://arxiv.org/abs/2511.12231)*,
  arXiv:2511.12231 (2025).
- X. Wang and O. Vafek, *Narrow bands in magnetic field and strong-coupling
  Hofstadter spectra*, Phys. Rev. B **106**, L121111 (2022).
- X. Wang and O. Vafek, *Revisiting Bloch electrons in a magnetic field:
  Hofstadter physics via hybrid Wannier states*, Phys. Rev. B **108**, 245109
  (2023).
- J. Wang, J. Cano, A. J. Millis, Z. Liu, and B. Yang, *Exact Landau Level
  Description of Geometry and Interaction in a Flatband*, Phys. Rev. Lett.
  **127**, 246403 (2021).
- P. J. Ledwith, A. Vishwanath, and D. E. Parker, *Vortexability: A Unifying
  Criterion for Ideal Fractional Chern Insulators*, Phys. Rev. B **108**,
  205144 (2023).
- Z. Ma *et al.*, *Topological flat bands in twisted trilayer graphene*,
  Science Bulletin **66**, 18–22 (2021).

BibTeX entries are collected in `CITATION.bib`.

## Tests

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

The LL installation checks shown above additionally exercise both new Python
entry points.
