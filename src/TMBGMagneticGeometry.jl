module TMBGMagneticGeometry

using LinearAlgebra

using ..TMBGMagneticHW
using ..TMBGCommonBasis

export magnetic_shift_maps,
       cartesian_metric_from_projector_distances

"""
Construct the reciprocal-space sewing maps required to lift every translated
hybrid-Wannier basis state into one common plane-wave Hilbert space.
"""
function magnetic_shift_maps(
    result::MagneticHWModel,
    translations::UnitRange{Int};
    reciprocal_shift_periodic::Bool=result.model.periodic_G,
)
    l1, q, nq = result.l1, result.flux.q, result.nq
    p_signed = result.flux.signB * result.flux.p
    shifts = Set{Tuple{Int,Int}}()
    for ir_row in 1:q, i1_row in 1:l1, translation in translations
        i1_unwrapped = i1_row - translation * p_signed * nq
        ir_unwrapped = ir_row - translation * p_signed
        push!(shifts, (
            div(i1_unwrapped - mod1(i1_unwrapped, l1), l1),
            div(ir_unwrapped - mod1(ir_unwrapped, q), q),
        ))
    end
    return Dict(
        shift => TMBGCommonBasis.reciprocal_shift_map(
            result.model,
            shift...;
            periodic=reciprocal_shift_periodic,
        )
        for shift in shifts
    )
end

"""
Recover the Cartesian Fubini--Study metric from three gauge-invariant
projector distances measured along two possibly nonorthogonal momentum steps.

`distance_difference` is the distance between the projectors at
`k + delta1` and `k + delta2`, so it probes `delta2 - delta1`.  The returned
tensor is expressed in Cartesian `x,y` coordinates; its ordinary matrix trace
is therefore the physical trace condition.
"""
function cartesian_metric_from_projector_distances(
    delta1::AbstractVector{<:Real},
    delta2::AbstractVector{<:Real},
    distance1::Real,
    distance2::Real,
    distance_difference::Real,
)
    basis_step = hcat(delta1, delta2)
    abs(det(basis_step)) > 1e-15 || error("Degenerate momentum steps")
    cross = (distance1 + distance2 - distance_difference) / 2
    directional_metric = [distance1 cross; cross distance2]
    inverse_step = inv(basis_step)
    return transpose(inverse_step) * directional_metric * inverse_step
end

end
