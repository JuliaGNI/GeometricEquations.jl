
using GeometricEquations
using Test

include("functions.jl")
include("initial_conditions.jl")


@testset "$(rpad("Stochastic Differential Equations (SDE)",80))" begin

    λ  = 2.
    μ  = 1.

    sde_v_params = (t, q, v) -> sde_v(λ, t, q, v)
    sde_B_params = (t, q, v) -> sde_B(μ, t, q, v)

    sde  = SDE(1, 1, sde_v_params, sde_B_params, t₀, x₀)
    sde1 = SDE(1, 1, sde_v_params, sde_B_params, x₀)
    sde2 = SDE(1, 3, sde_v_params, sde_B_params, x₀)
    sde3 = SDE(1, 1, sde_v_params, sde_B_params, x₁ₛ)

    @test ndims(sde) == 2
    @test periodicity(sde) == zero(x₀)
    @test functions(sde) == NamedTuple{(:v,:B)}((sde_v_params, sde_B_params))

    @test sde == sde1
    @test sde != sde2
    @test sde != sde3

    @test hash(sde) == hash(sde1)
    @test hash(sde) != hash(sde2)
    @test hash(sde) != hash(sde3)

    @test sde1.d == 2
    @test sde2.d == 2
    @test sde3.d == 2

    @test sde1.ns == 1
    @test sde2.ns == 3
    @test sde3.ns == 1

    @test nsamples(sde1) == 1
    @test nsamples(sde2) == 1
    @test nsamples(sde3) == 3

    @test sde == similar(sde, t₀, x₀, 1)
    @test sde == similar(sde, t₀, x₀)
    @test sde == similar(sde, x₀, 1)
    @test sde == similar(sde, x₀)

end


@testset "$(rpad("Partitioned Stochastic Differential Equations (PSDE)",80))" begin

    psde  = PSDE(1, 1, psde_v, psde_f, psde_B, psde_G, t₀, q₀, p₀)
    psde1 = PSDE(1, 1, psde_v, psde_f, psde_B, psde_G, q₀, p₀)
    psde2 = PSDE(1, 3, psde_v, psde_f, psde_B, psde_G, q₀, p₀)
    psde3 = PSDE(1, 1, psde_v, psde_f, psde_B, psde_G, q₁ₛ, p₁ₛ)

    @test ndims(psde) == 1
    @test periodicity(psde) == zero(q₀)
    @test functions(psde) == NamedTuple{(:v,:f,:B,:G)}((psde_v,psde_f,psde_B,psde_G))

    @test psde == psde1
    @test psde != psde2
    @test psde != psde3

    @test hash(psde) == hash(psde1)
    @test hash(psde) != hash(psde2)
    @test hash(psde) != hash(psde3)

    @test psde1.d == 1
    @test psde2.d == 1
    @test psde3.d == 1

    @test psde1.ns == 1
    @test psde2.ns == 3
    @test psde3.ns == 1

    @test nsamples(psde1) == 1
    @test nsamples(psde2) == 1
    @test nsamples(psde3) == 3

    @test psde == similar(psde, t₀, q₀, p₀, 1)
    @test psde == similar(psde, q₀, p₀, 1)

    @test psde != similar(psde, t₀, q₁ₛ, p₁ₛ)
    @test psde != similar(psde, q₁ₛ, p₁ₛ)

end


@testset "$(rpad("Split Partitioned Stochastic Differential Equations (SPSDE)",80))" begin

    spsde  = SPSDE(1, 1, spsde_v, spsde_f1, spsde_f2, spsde_B, spsde_G1, spsde_G2, t₀, q₀, p₀)
    spsde1 = SPSDE(1, 1, spsde_v, spsde_f1, spsde_f2, spsde_B, spsde_G1, spsde_G2, q₀, p₀)
    spsde2 = SPSDE(1, 3, spsde_v, spsde_f1, spsde_f2, spsde_B, spsde_G1, spsde_G2, q₀, p₀)
    spsde3 = SPSDE(1, 1, spsde_v, spsde_f1, spsde_f2, spsde_B, spsde_G1, spsde_G2, q₁ₛ, p₁ₛ)

    @test ndims(spsde) == 1
    @test periodicity(spsde) == zero(q₀)
    @test functions(spsde) == NamedTuple{(:v,:f1,:f2,:B,:G1,:G2)}((spsde_v,spsde_f1,spsde_f2,spsde_B,spsde_G1,spsde_G2))

    @test spsde == spsde1
    @test spsde != spsde2
    @test spsde != spsde3

    @test hash(spsde) == hash(spsde1)
    @test hash(spsde) != hash(spsde2)
    @test hash(spsde) != hash(spsde3)

    @test spsde1.d == 1
    @test spsde2.d == 1
    @test spsde3.d == 1

    @test spsde1.ns == 1
    @test spsde2.ns == 3
    @test spsde3.ns == 1

    @test nsamples(spsde1) == 1
    @test nsamples(spsde2) == 1
    @test nsamples(spsde3) == 3

    @test spsde == similar(spsde, t₀, q₀, p₀, 1)
    @test spsde == similar(spsde, q₀, p₀, 1)

    @test spsde != similar(spsde, t₀, q₁ₛ, p₁ₛ)
    @test spsde != similar(spsde, q₁ₛ, p₁ₛ)

end
