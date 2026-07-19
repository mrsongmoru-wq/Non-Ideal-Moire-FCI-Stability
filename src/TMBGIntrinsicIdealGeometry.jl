module TMBGIntrinsicIdealGeometry

using LinearAlgebra

export polar_frame_transport,
       c3_kahler_metric

"""
Return the unitary polar transport and singular values of an inter-frame
overlap.  The polar transport is covariant under independent unitary gauge
rotations of the frames at the two endpoints.
"""
function polar_frame_transport(overlap::AbstractMatrix{<:Complex})
    factorization = svd(overlap)
    return factorization.U * factorization.Vt, factorization.S
end

"""
The unique Cartesian ideal Kähler metric compatible with C3 symmetry and a
rank-one Berry curvature `curvature`.  It saturates
`Tr(g) = abs(curvature)` pointwise.
"""
function c3_kahler_metric(curvature::Real)
    scale = abs(curvature) / 2
    return [scale 0.0; 0.0 scale]
end

end
