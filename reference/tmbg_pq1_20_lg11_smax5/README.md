# Validated tMBG reference calculation

These files are the lightweight outputs of
`examples/run_tmbg_stability_check.jl` for the fixed tMBG case
`p/q=1/20`, `lg=11`, `smax=5`, and six retained zero-field bands.  The full
serialized magnetic hybrid-Wannier cache is about 567 MB and is regenerated
locally rather than committed.

The isolated target magnetic band is `79`, with multiband correction frame
`(79,78)`.  The saved trace-condition result is:

- target Chern number `1`;
- correction-frame Chern number `2`;
- ideal integrated trace `6.283185307179585`;
- residual integrated trace and trace excess
  `0.00022470507632712814`;
- total integrated trace `6.283410012255914`;
- ideal fraction `0.9999642383553055`;
- minimum frame-link singular value `0.8379962538093386`.

The independently resolved target contribution occupies
`0.919390352703438` of the central zero-field pair.  C3 group projection
reduces the scalar map's relative C3 error from `0.12447892706843225` to
`7.209450195207299e-17` without changing its integrated weight.

`magnetic_spectrum.csv` is the complete numerical spectrum; the accompanying
PDF and PNG show the `15--40 meV` window.  `quantum_geometry.csv` and the two
overlap files preserve all local numerical diagnostics.  The
`ideal_component_projection/` directory contains the momentum-resolved
projection tables, zero-field path spectrum, summary, and the corresponding
two-panel visualization.
