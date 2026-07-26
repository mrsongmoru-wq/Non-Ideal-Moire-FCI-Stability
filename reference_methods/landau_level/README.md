# Conventional Landau-level comparison

This directory contains an independent tMBG Landau-level expansion used as a
standard-method reference for the hybrid-Wannier results.

It is intentionally separate from the primary workflow:

- no LL state is passed to the hybrid-Wannier calculation;
- the LL result is not used to fit or select a hybrid-Wannier wavefunction;
- the comparison is made only after both calculations have produced their
  physical spectrum and quantum-geometry observables.

## Files

- `tmbg_landau_level.py` — tMBG continuum Hamiltonian in the full
  layer/sublattice/Landau-orbital basis, with magnetic translations and
  Landau-orbital form factors;
- `compute_ll_magnetic_spectrum.py` — fixed-flux magnetic-spectrum export;
- `plot_ll_magnetic_spectrum.py` — magnetic-spectrum visualization;
- `compute_ll_quantum_geometry.py` — rank-one or multiband Berry curvature,
  Quantum Metric, Chern number, and trace-condition calculation;
- `plot_ll_quantum_geometry.py` — three-panel geometry visualization.

## Installation check

From the repository root:

```bash
python -m pip install -r requirements-ll.txt
python reference_methods/landau_level/compute_ll_magnetic_spectrum.py --quick
python reference_methods/landau_level/compute_ll_quantum_geometry.py --quick
python reference_methods/landau_level/plot_ll_magnetic_spectrum.py
```

The quick mode confirms that the Hamiltonian, diagonalization, export, and
geometry paths work.  It is not a convergence setting.

## Reproduce the versioned geometry reference

```bash
python reference_methods/landau_level/compute_ll_quantum_geometry.py \
  --nll 400 --nk1 40 --nk2 16 --band-start 82 \
  --geometry-backend block --jobs 8 \
  --output results/ll_reference/quantum_geometry.npz

python reference_methods/landau_level/plot_ll_quantum_geometry.py \
  --input results/ll_reference/quantum_geometry.npz
```

The selected target is band 82 in the central 120-level ordering.  A
multiband projector is selected by increasing `--nbands`; all chosen states
remain embedded in the full LL basis.

For this tMBG lattice convention, LL input flux `p/q=1/10` corresponds to a
physical flux of `1/20` per moiré cell.  See `docs/technical_method.md` for the
basis and geometry conventions.

The converged lightweight output and summary are stored in
`reference/tmbg_pq1_20_ll_nll400_40x16/`.
