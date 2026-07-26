# Technical method

This document collects the basis conventions and geometry formulas that are
intentionally kept out of the usage-oriented main README.

## 1. Primary hybrid-Wannier calculation

### Zero-field and magnetic bases

The tMBG continuum Hamiltonian has three layers and two sublattices, hence six
local orbitals for every plane-wave momentum.  These six local orbitals are
not the same object as the six retained zero-field energy bands.

The zero-field problem is solved on a uniform Brillouin-zone mesh.  Non-Abelian
parallel transport along the first reciprocal direction constructs
hybrid-Wannier states.  Magnetic translations of the retained bands form a
nonorthogonal basis with overlap matrix \(O(\mathbf{k})\) and Hamiltonian
\(H(\mathbf{k})\).  After restricting to the physical rank, the code solves
the orthonormalized problem

\[
\widetilde H(\mathbf{k})
=O(\mathbf{k})^{-1/2}H(\mathbf{k})O(\mathbf{k})^{-1/2}.
\]

The production calculation uses a reciprocal plane-wave grid `lg=11`, giving
\(6\times11^2=726\) local plane-wave orbitals at each momentum, and a
magnetic-translation cutoff `smax=5`.

### Covariant target projector

The final observable is associated with one isolated magnetic subband, so its
physical projector has rank one.  Nearby subbands enter only through a local
multiband correction frame \(\Phi(\mathbf{k})\).  For neighboring momenta the
frame link is evaluated in the common plane-wave basis,

\[
M=\Phi(\mathbf{k})^\dagger
T(\delta\mathbf{k})\Phi(\mathbf{k}+\delta\mathbf{k}).
\]

If \(M=WSV^\dagger\), the polar factor \(U=WV^\dagger\) defines covariant
transport.  The target link is taken from the target entry of \(U\).  This
preserves a rank-one target projector while removing nonunitary distortion
from the finite correction frame.

### Physical Cartesian trace

The finite differences are taken along two reciprocal-lattice directions
that are not Cartesian axes.  The code evaluates projector distances along
both steps and their difference, reconstructs the full metric tensor, and
then transforms it to physical Cartesian coordinates.  The reported trace is

\[
\mathrm{Tr}\,g=g_{xx}+g_{yy},
\]

not \(g_{11}+g_{22}\) in oblique reciprocal coordinates.

### Projection back to zero field

The selected magnetic state is lifted into the complete six-local-orbital
plane-wave Hilbert space.  A unitary real-space phase converts the
hybrid-Wannier Landau-gauge state to symmetric gauge.  The result is projected
onto the true zero-field energy eigenstates.  A final C3 group average acts on
the scalar momentum-resolved weight and preserves its integrated weight.

## 2. Geometry reported by both implementations

For a rank-one projector, the quantum geometric tensor is

\[
Q_{\mu\nu}(\mathbf{k})=
\langle \partial_{k_\mu}u_{\mathbf{k}}|
(1-|u_{\mathbf{k}}\rangle\langle u_{\mathbf{k}}|)
|\partial_{k_\nu}u_{\mathbf{k}}\rangle .
\]

Its real part is the Quantum Metric \(g_{\mu\nu}\), while its antisymmetric
imaginary part gives the Berry curvature \(\Omega\).  The pointwise trace
bound is

\[
\mathrm{Tr}\,g(\mathbf{k})\geq |\Omega(\mathbf{k})|.
\]

The integrated trace excess used in this repository is

\[
\eta=\int_{\mathrm{BZ}}
\left[\mathrm{Tr}\,g(\mathbf{k})-|\Omega(\mathbf{k})|\right]\,d^2k.
\]

Small \(\eta\), local trace saturation, the Chern number, curvature
uniformity, and overlap quality are reported together.  No single number is
treated as proof of an interacting FCI.

For a rank-\(r\) group of bands, the same construction is applied to the
projector \(P(\mathbf{k})=\sum_{a=1}^{r}|u_a\rangle\langle u_a|\).  Links,
Wilson plaquettes, and projector distances are then non-Abelian and invariant
under a momentum-dependent \(U(r)\) rotation within the selected subspace.

## 3. Conventional Landau-level comparison

The LL reference uses the basis

\[
|\mathrm{layer},\mathrm{sublattice},n,j\rangle,
\qquad
n=0,\ldots,N_{\mathrm{LL}},
\quad j=0,\ldots,p-1,
\]

with total dimension \(6(N_{\mathrm{LL}}+1)p\).  Thus the geometry is not
computed from a truncated list of energy bands alone: the wavefunctions retain
the complete layer, sublattice, and Landau-orbital basis at the chosen
\(N_{\mathrm{LL}}\).

Neighboring-momentum overlaps include analytic Landau-orbital form factors.
The Wilson plaquette also contains the Landau-gauge background factor

\[
\exp\!\left[-\frac{i2\pi}{pN_{k_1}N_{k_2}}\right].
\]

Metric elements are obtained from projector distances along \(x\), \(y\), and
the mixed displacement.  The LL coordinate system used here is Cartesian, so
the trace is directly \(g_{xx}+g_{yy}\).  The code supports a rank-one target
or a selected multiband projector through `--nbands`.

The LL implementation uses a rectangular magnetic-cell convention in which
the input `p/q=1/10` represents the same physical moiré-cell flux `1/20` used
by the hybrid-Wannier calculation.  This factor of two is purely a unit-cell
area convention.

The LL result is a conventional finite-field reference for methodological
comparison.  It does not calibrate the hybrid-Wannier target state and is not
read by the primary pipeline.

## 4. Convergence controls

The two implementations have independent convergence parameters:

| Implementation | Basis convergence | Momentum convergence |
| --- | --- | --- |
| Hybrid Wannier | plane-wave grid `lg`; magnetic translations `smax`; retained correction-frame bands | magnetic momentum mesh |
| Landau level | maximum orbital `NLL`; selected projector rank | `nk1 x nk2` magnetic momentum mesh |

The versioned production results use `lg=11`, `smax=5` for the hybrid-Wannier
calculation and `NLL=400`, `40 x 16` for the LL comparison.  The `--quick`
commands only verify installation and data flow.

## 5. Numerical output conventions

The hybrid-Wannier workflow writes human-readable CSV and text files.  The LL
geometry calculation writes a compressed NumPy archive containing:

- `energies`;
- `berry_curvature`;
- `metric_trace`;
- `trace_violation`;
- `determinant_condition`;
- `chern`, `eta`, and integrated curvature/trace values;
- mesh, flux, band, and cutoff metadata.

The corresponding text summary is generated beside the archive.
