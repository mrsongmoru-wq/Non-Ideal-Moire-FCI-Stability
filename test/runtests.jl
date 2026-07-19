using Test
using LinearAlgebra

using TMBGMagneticHybridWannier
using TMBGMagneticHybridWannier.TMBGMagneticHW
using TMBGMagneticHybridWannier.TMBGProjection
using TMBGMagneticHybridWannier.TMBGCommonBasis
using TMBGMagneticHybridWannier.TMBGMagneticGeometry
using TMBGMagneticHybridWannier.TMBGIntrinsicIdealGeometry
using TMBGMagneticHybridWannier.TMBGSymmetricGaugeProjection
using TMBGMagneticHybridWannier.TMBGMagneticHW.TMBGHybridWannier
using TMBGMagneticHybridWannier.TMBGMagneticHW.TMBGHybridWannier.TMBGZeroField
using TMBGMagneticHybridWannier.TMBGMagneticHW.TMBGHybridWannier.WangContinuum

const PARAMS = TMBGParams(1.04, 3, 0.7, 1, -30.0, (1, 0, 0), 50.0)

@testset "Wang reciprocal frame" begin
    model = build_tmbg_wang(PARAMS; lg=7, periodic_G=false)
    reciprocal_angle = acosd(dot(model.g1, model.g2) / norm(model.g1) / norm(model.g2))
    direct_angle = acosd(dot(model.a1, model.a2) / norm(model.a1) / norm(model.a2))
    @test reciprocal_angle ≈ 120.0 atol=1e-12
    @test direct_angle ≈ 60.0 atol=1e-12
    @test norm(model.g1 - (build_model(PARAMS).g1 - build_model(PARAMS).g2)) < 1e-14
end

@testset "Wang high-symmetry coordinates" begin
    source = build_model(PARAMS)
    model = build_tmbg_wang(PARAMS; lg=7, periodic_G=false)
    source_points = (
        K=(1 / 3, 1 / 3),
        Gamma=(2 / 3, -1 / 3),
        M=(1 / 6, 1 / 6),
        Kprime=(0.0, 0.0),
    )
    expected_wang = (
        K=(1 / 3, 2 / 3),
        Gamma=(2 / 3, 1 / 3),
        M=(1 / 6, 1 / 3),
        Kprime=(0.0, 0.0),
    )
    for label in keys(source_points)
        source_reduced = collect(source_points[label])
        wang_reduced = collect(expected_wang[label])
        source_cartesian = source_reduced[1] .* source.g1 .+
                           source_reduced[2] .* source.g2
        wang_cartesian = wang_reduced[1] .* model.g1 .+
                         wang_reduced[2] .* model.g2
        @test norm(source_cartesian - wang_cartesian) < 1e-14
    end
    gamma = collect(expected_wang.Gamma)
    rotated_gamma = [-gamma[2], gamma[1] - gamma[2]]
    @test norm(rotated_gamma - gamma - [-1.0, 0.0]) < 1e-14
end

@testset "Zero-field basis equivalence" begin
    source = build_model(PARAMS)
    model = build_tmbg_wang(PARAMS; lg=7, periodic_G=false)
    for kred in ([0.0, 0.0], [0.2, 0.3], [1 / 3, 1 / 3])
        k = kred[1] .* model.g1 .+ kred[2] .* model.g2
        e_source = eigvals(Hermitian(TMBGZeroField.hamiltonian(source, k)))[110:113]
        e_wang_all = eigvals(Hermitian(WangContinuum.hamiltonian(model, k)))
        mid = length(e_wang_all) ÷ 2
        e_wang = e_wang_all[(mid - 1):(mid + 2)]
        @test maximum(abs.(e_source - e_wang)) < 3e-3
    end
end

@testset "Central bands preserve global indices" begin
    model = build_tmbg_wang(PARAMS; lg=7, periodic_G=false)
    k = [0.0, 0.0]
    hall = WangContinuum.hamiltonian(model, k)
    full = eigen(Hermitian(hall))
    mid = size(hall, 1) ÷ 2
    vals, vecs = central_eigensystem(model, k, 2)
    @test vals ≈ full.values[mid:(mid + 1)] atol=1e-10
    @test abs.(diag(vecs' * full.vectors[:, mid:(mid + 1)])) ≈ ones(2) atol=1e-8
end

@testset "All momentum-linear matrix elements" begin
    model = build_tmbg_wang(PARAMS; lg=5, periodic_G=false)
    k = [0.0123, -0.0071]
    delta = 1e-7
    finite_difference = (WangContinuum.hamiltonian(model, k .+ [0, delta]) -
                         WangContinuum.hamiltonian(model, k .- [0, delta])) / (2delta)
    analytic = kron(Matrix{ComplexF64}(I, length(model.Glist), length(model.Glist)),
                    velocity_y(model))
    @test norm(finite_difference - analytic) / norm(analytic) < 1e-9
end

@testset "Finite-box x matrix element" begin
    L = 17.0
    for n in (-3, -1, 1, 4)
        delta = 2pi * n / L
        expected = -im * L * (-1)^n / (2pi * n)
        @test finite_box_x_integral(delta, L) ≈ expected atol=1e-13
    end
    @test finite_box_x_integral(0.0, L) == 0
end

@testset "Flux convention" begin
    flux = FluxSpec(1, 20, 1)
    @test physical_flux(flux) == 1 // 20
    nq, l1, l2, _, _, _ = magnetic_grid(
        flux;
        mesh_factor1=2,
        mesh_factor2=4,
    )
    @test (nq, l1, l2) == (2, 80, 8)
    model = build_tmbg_wang(PARAMS; lg=5, periodic_G=false)
    @test field_coupling_sign(model, 1) == -1
    @test field_coupling_sign(model, -1) == 1
end

@testset "Nonorthogonal projector-distance metric" begin
    delta1 = [0.07, 0.0]
    delta2 = [-0.02, 0.04]
    expected = [2.0 0.3; 0.3 1.4]
    distance1 = dot(delta1, expected * delta1)
    distance2 = dot(delta2, expected * delta2)
    difference = delta2 - delta1
    distance_difference = dot(difference, expected * difference)
    recovered = cartesian_metric_from_projector_distances(
        delta1,
        delta2,
        distance1,
        distance2,
        distance_difference,
    )
    @test recovered ≈ expected atol=1e-13
    @test tr(recovered) ≈ expected[1, 1] + expected[2, 2] atol=1e-13
end

@testset "Intrinsic ideal-component geometry" begin
    left = qr(randn(ComplexF64, 4, 4)).Q |> Matrix
    right = qr(randn(ComplexF64, 4, 4)).Q |> Matrix
    overlap = randn(ComplexF64, 4, 4)
    transport, singular_values = polar_frame_transport(overlap)
    rotated_transport, rotated_singular_values =
        polar_frame_transport(left * overlap * right)
    @test rotated_transport ≈ left * transport * right atol=1e-12
    @test sort(rotated_singular_values) ≈ sort(singular_values) atol=1e-12

    curvature = 2.75
    metric = c3_kahler_metric(curvature)
    @test tr(metric) ≈ abs(curvature)
    @test det(metric) ≈ curvature^2 / 4

end

@testset "Chern-counted upper-group selection" begin
    spectrum = repeat([-5.0, 10.0, 18.0, 22.0, 25.05, 35.0], 1, 2)
    target, upper = select_target_upper_group(
        spectrum, 4; target_chern=-2, upper_reference_energy=25.047,
    )
    @test target == 4:5
    @test upper == [5]
end

@testset "C3 projection-grid diagnostic" begin
    n = 12
    symmetric = zeros(Float64, n, n)
    for i2 in 0:(n - 1), i1 in 0:(n - 1)
        orbit = Tuple{Int,Int}[]
        point = (i1, i2)
        for _ in 1:3
            push!(orbit, point)
            point = (mod(-point[2], n), mod(point[1] - point[2], n))
        end
        value = sum(first(x) + 2 * last(x) for x in orbit)
        symmetric[i1 + 1, i2 + 1] = value
    end
    @test c3_relative_error(symmetric) < 1e-14
    symmetric[2, 1] += 1
    @test c3_relative_error(symmetric) > 0
    restored = c3_symmetrize(symmetric)
    @test c3_relative_error(restored) < 1e-14
    @test sum(restored) ≈ sum(symmetric) atol=1e-12
    @test norm(restored - symmetric) > 0
end

@testset "Unitary real-space gauge-transform convention" begin
    fourier = TMBGSymmetricGaugeProjection.unitary_dft(9)
    @test norm(fourier' * fourier - I) < 1e-13
end

@testset "Strict magnetic trial matrices" begin
    model = build_tmbg_wang(PARAMS; lg=5, periodic_G=true)
    flux = FluxSpec(1, 6, 1)
    nq, l1, l2, k1, _, r = magnetic_grid(flux; mesh_factor=1)
    hw = build_hybrid_wannier(model, l1, l1; nactive=4)
    @test hw.chern == -2
    recovered_energy_error = 0.0
    recovered_orthonormality_error = 0.0
    transform_unitarity_error = 0.0
    frame_transform_difference = 0.0
    for i2 in 1:l1, i1 in 1:l1
        values, vectors = energy_eigenbasis(hw, i1, i2)
        transform = hybrid_to_energy_transform(hw, i1, i2)
        recovered_energy_error = max(
            recovered_energy_error,
            maximum(abs.(values .- hw.eigvals[:, i1, i2])),
        )
        recovered_orthonormality_error = max(
            recovered_orthonormality_error,
            norm(vectors' * vectors - I),
        )
        transform_unitarity_error = max(
            transform_unitarity_error,
            norm(transform' * transform - I),
        )
        frame_transform_difference = max(
            frame_transform_difference,
            norm(transform - hw.coeffs[:, :, i1, i2]),
        )
    end
    @test recovered_energy_error < 1e-10
    @test recovered_orthonormality_error < 1e-12
    @test transform_unitarity_error < 1e-12
    @test frame_transform_difference > 1e-3
    O, H0, HB = build_trial_matrices(
        model, hw, flux, nq, l1, l2, k1, r;
        nvec=[0, 1], svec=collect(-3:3), gc=8,
    )
    relherm(A) = maximum(norm(A[:, :, i] - A[:, :, i]') /
                         max(norm(A[:, :, i]), eps()) for i in axes(A, 3))
    @test relherm(O) < 1e-12
    @test relherm(H0 + HB) < 0.04
    @test norm(HB) > 0
end
