module TMBGProjection

using LinearAlgebra
using Statistics

export select_target_upper_group, c3_relative_error

function select_target_upper_group(
    spectrum::Array{Float64,2},
    q::Int;
    target_chern::Int=-2,
    upper_reference_energy::Float64=25.047,
)
    means = vec(mean(spectrum; dims=2))
    target_count = q + target_chern
    target_count > 1 || error("Target magnetic manifold must contain at least two branches")
    candidates = Int[]
    for last in target_count:length(means)
        23.0 <= means[last] <= 27.0 && push!(candidates, last)
    end
    isempty(candidates) &&
        error("No target-manifold endpoint found near the reference energy")
    last = candidates[argmin(abs.(means[candidates] .- upper_reference_energy))]
    target_block = (last - target_count + 1):last
    split = argmax(diff(means[target_block]))
    upper_group = collect((first(target_block) + split):last)
    return target_block, upper_group
end

function c3_relative_error(weights::AbstractMatrix{<:Real})
    n1, n2 = size(weights)
    n1 == n2 || error("C3 check requires a square full-Brillouin-zone grid")
    rotated = similar(weights)
    for i2 in 0:(n2 - 1), i1 in 0:(n1 - 1)
        j1 = mod(-i2, n1)
        j2 = mod(i1 - i2, n2)
        rotated[j1 + 1, j2 + 1] = weights[i1 + 1, i2 + 1]
    end
    return norm(weights - rotated) / max(norm(weights), eps(Float64))
end

end
