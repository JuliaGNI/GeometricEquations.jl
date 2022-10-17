
using GeometricEquations
using Test

include("functions.jl")
include("initial_conditions.jl")


@testset "$(rpad("Stochastic Differential Equations (SDE)",80))" begin

    sde  = SDE(sde_v, sde_B, NullInvariants(), NullParameters(), NullPeriodicity())
    sde1 = SDE(sde_v, sde_B)
    sde2 = SDE(sde_v, sde_B; invariants=NullInvariants())
    sde3 = SDE(sde_v, sde_B; parameters=NullParameters())
    sde4 = SDE(sde_v, sde_B; periodicity=NullPeriodicity())

    @test sde == sde1
    @test sde == sde2
    @test sde == sde3
    @test sde == sde4

    @test hash(sde) == hash(sde1)
    @test hash(sde) == hash(sde2)
    @test hash(sde) == hash(sde3)
    @test hash(sde) == hash(sde4)

    @test functions(sde) == NamedTuple{(:v,:B)}((sde_v,sde_B))
    @test solutions(sde) == NamedTuple()

    @test parameters(sde) == NullParameters()
    @test invariants(sde) == NullInvariants()
    @test periodicity(sde) == NullPeriodicity()

    @test hasvectorfield(sde) == true
    @test hassolution(sde) == false
    @test hasprimary(sde) == false
    @test hassecondary(sde) == false

    @test hasinvariants(sde) == false
    @test hasparameters(sde) == false
    @test hasperiodicity(sde) == false

    @test hashamiltonian(sde) == false
    @test haslagrangian(sde) == false


    sde  = SDE(sde_v, sde_B, NullInvariants(), sde_param_types, NullPeriodicity())

    funcs = functions(sde)

    @test_nowarn funcs.v(zero(x₀), t₀, x₀, sde_params)

    funcs = functions(sde, sde_params)

    @test_nowarn funcs.v(zero(x₀), t₀, x₀)

end


@testset "$(rpad("Partitioned Stochastic Differential Equations (PSDE)",80))" begin

    psde  = PSDE(psde_v, psde_f, psde_B, psde_G, NullInvariants(), NullParameters(), NullPeriodicity())
    psde1 = PSDE(psde_v, psde_f, psde_B, psde_G)
    psde2 = PSDE(psde_v, psde_f, psde_B, psde_G; invariants=NullInvariants())
    psde3 = PSDE(psde_v, psde_f, psde_B, psde_G; parameters=NullParameters())
    psde4 = PSDE(psde_v, psde_f, psde_B, psde_G; periodicity=NullPeriodicity())

    @test psde == psde1
    @test psde == psde2
    @test psde == psde3
    @test psde == psde4

    @test hash(psde) == hash(psde1)
    @test hash(psde) == hash(psde2)
    @test hash(psde) == hash(psde3)
    @test hash(psde) == hash(psde4)

    @test functions(psde) == NamedTuple{(:v,:f,:B,:G)}((psde_v,psde_f,psde_B,psde_G))
    @test solutions(psde) == NamedTuple()

    @test parameters(psde) == NullParameters()
    @test invariants(psde) == NullInvariants()
    @test periodicity(psde) == NullPeriodicity()

    @test hasvectorfield(psde) == true
    @test hassolution(psde) == false
    @test hasprimary(psde) == false
    @test hassecondary(psde) == false

    @test hasinvariants(psde) == false
    @test hasparameters(psde) == false
    @test hasperiodicity(psde) == false

    @test hashamiltonian(psde) == false
    @test haslagrangian(psde) == false


    psde = PSDE(psde_v, psde_f, psde_B, psde_G, NullInvariants(), sde_param_types, NullPeriodicity())

    funcs = functions(psde)

    @test_nowarn funcs.v(zero(q₀), t₀, q₀, p₀, sde_params)
    @test_nowarn funcs.f(zero(p₀), t₀, q₀, p₀, sde_params)

    funcs = functions(psde, sde_params)

    @test_nowarn funcs.v(zero(q₀), t₀, q₀, p₀)
    @test_nowarn funcs.f(zero(p₀), t₀, q₀, p₀)

end


@testset "$(rpad("Split Partitioned Stochastic Differential Equations (SPSDE)",80))" begin

    psde  = SPSDE(spsde_v, spsde_f1, spsde_f2, spsde_B, spsde_G1, spsde_G2, NullInvariants(), NullParameters(), NullPeriodicity())
    psde1 = SPSDE(spsde_v, spsde_f1, spsde_f2, spsde_B, spsde_G1, spsde_G2)
    psde2 = SPSDE(spsde_v, spsde_f1, spsde_f2, spsde_B, spsde_G1, spsde_G2; invariants=NullInvariants())
    psde3 = SPSDE(spsde_v, spsde_f1, spsde_f2, spsde_B, spsde_G1, spsde_G2; parameters=NullParameters())
    psde4 = SPSDE(spsde_v, spsde_f1, spsde_f2, spsde_B, spsde_G1, spsde_G2; periodicity=NullPeriodicity())

    @test psde == psde1
    @test psde == psde2
    @test psde == psde3
    @test psde == psde4

    @test hash(psde) == hash(psde1)
    @test hash(psde) == hash(psde2)
    @test hash(psde) == hash(psde3)
    @test hash(psde) == hash(psde4)

    @test functions(psde) == NamedTuple{(:v,:f1,:f2,:B,:G1,:G2)}((spsde_v,spsde_f1,spsde_f2,spsde_B,spsde_G1,spsde_G2))
    @test solutions(psde) == NamedTuple()

    @test parameters(psde) == NullParameters()
    @test invariants(psde) == NullInvariants()
    @test periodicity(psde) == NullPeriodicity()

    @test hasvectorfield(psde) == true
    @test hassolution(psde) == false
    @test hasprimary(psde) == false
    @test hassecondary(psde) == false

    @test hasinvariants(psde) == false
    @test hasparameters(psde) == false
    @test hasperiodicity(psde) == false

    @test hashamiltonian(psde) == false
    @test haslagrangian(psde) == false


    psde = SPSDE(spsde_v, spsde_f1, spsde_f2, spsde_B, spsde_G1, spsde_G2, NullInvariants(), sde_param_types, NullPeriodicity())

    funcs = functions(psde)

    @test_nowarn funcs.v(zero(q₀), t₀, q₀, p₀, sde_params)
    @test_nowarn funcs.f1(zero(p₀), t₀, q₀, p₀, sde_params)
    @test_nowarn funcs.f2(zero(p₀), t₀, q₀, p₀, sde_params)

    funcs = functions(psde, sde_params)

    @test_nowarn funcs.v(zero(q₀), t₀, q₀, p₀)
    @test_nowarn funcs.f1(zero(p₀), t₀, q₀, p₀, )
    @test_nowarn funcs.f2(zero(p₀), t₀, q₀, p₀, )

end
