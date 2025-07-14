
using GeometricEquations
using Test

include("functions.jl")
include("initial_conditions.jl")


@testset "$(rpad("Stochastic Differential Equations (SDE)",80))" begin

    sde  = SDE(sde_v, sde_B, TestNoise(), NullInvariants(), NullParameters(), NullPeriodicity())
    sde1 = SDE(sde_v, sde_B, TestNoise())
    sde2 = SDE(sde_v, sde_B, TestNoise(); invariants=NullInvariants())
    sde3 = SDE(sde_v, sde_B, TestNoise(); parameters=NullParameters())
    sde4 = SDE(sde_v, sde_B, TestNoise(); periodicity=NullPeriodicity())

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


    sde  = SDE(sde_v, sde_B, TestNoise(), NullInvariants(), sde_param_types, NullPeriodicity())

    funcs = functions(sde)

    @test_nowarn funcs.v(zero(x₀), t₀, x₀, sde_params)

    funcs = functions(sde, sde_params)

    @test_nowarn funcs.v(zero(x₀), t₀, x₀)


    # Test for periodicity
    sdep = SDE(sde_v, sde_B, TestNoise(); periodicity=([-π,0],[+π,2π]))

    @test periodicity(sdep) == (Float64[-π,0],Float64[+π,2π])
    @test getperiodicity(sdep) == BitArray([true,true])
    @test hasperiodicity(sdep) == true

    sdep = SDE(sde_v, sde_B, TestNoise(); periodicity=([-Inf,0],[+Inf,2π]))

    @test periodicity(sdep) == (Float64[-Inf,0],Float64[+Inf,2π])
    @test getperiodicity(sdep) == BitArray([false,true])
    @test hasperiodicity(sdep) == true

    sdep = SDE(sde_v, sde_B, TestNoise(); periodicity=([-Inf,-Inf],[+Inf,+Inf]))

    @test periodicity(sdep) == NullPeriodicity()
    @test ismissing(getperiodicity(sdep))
    @test hasperiodicity(sdep) == false

end


@testset "$(rpad("Partitioned Stochastic Differential Equations (PSDE)",80))" begin

    psde  = PSDE(psde_v, psde_f, psde_B, psde_G, TestNoise(), NullInvariants(), NullParameters(), NullPeriodicity())
    psde1 = PSDE(psde_v, psde_f, psde_B, psde_G, TestNoise())
    psde2 = PSDE(psde_v, psde_f, psde_B, psde_G, TestNoise(); invariants=NullInvariants())
    psde3 = PSDE(psde_v, psde_f, psde_B, psde_G, TestNoise(); parameters=NullParameters())
    psde4 = PSDE(psde_v, psde_f, psde_B, psde_G, TestNoise(); periodicity=NullPeriodicity())

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


    psde = PSDE(psde_v, psde_f, psde_B, psde_G, TestNoise(), NullInvariants(), sde_param_types, NullPeriodicity())

    funcs = functions(psde)

    @test_nowarn funcs.v(zero(q₀), t₀, q₀, p₀, sde_params)
    @test_nowarn funcs.f(zero(p₀), t₀, q₀, p₀, sde_params)

    funcs = functions(psde, sde_params)

    @test_nowarn funcs.v(zero(q₀), t₀, q₀, p₀)
    @test_nowarn funcs.f(zero(p₀), t₀, q₀, p₀)


    # Test for periodicity
    psdep = PSDE(psde_v, psde_f, psde_B, psde_G, TestNoise(); periodicity=([0.0,],[2π]))

    @test periodicity(psdep) == (Float64[0],Float64[2π])
    @test getperiodicity(psdep) == BitArray([true])
    @test hasperiodicity(psdep) == true

    psdep = PSDE(psde_v, psde_f, psde_B, psde_G, TestNoise(); periodicity=([-Inf],[+Inf]))

    @test periodicity(psdep) == NullPeriodicity()
    @test ismissing(getperiodicity(psdep))
    @test hasperiodicity(psdep) == false

end


@testset "$(rpad("Split Partitioned Stochastic Differential Equations (SPSDE)",80))" begin

    psde  = SPSDE(spsde_v, spsde_f1, spsde_f2, spsde_B, spsde_G1, spsde_G2, TestNoise(), NullInvariants(), NullParameters(), NullPeriodicity())
    psde1 = SPSDE(spsde_v, spsde_f1, spsde_f2, spsde_B, spsde_G1, spsde_G2, TestNoise())
    psde2 = SPSDE(spsde_v, spsde_f1, spsde_f2, spsde_B, spsde_G1, spsde_G2, TestNoise(); invariants=NullInvariants())
    psde3 = SPSDE(spsde_v, spsde_f1, spsde_f2, spsde_B, spsde_G1, spsde_G2, TestNoise(); parameters=NullParameters())
    psde4 = SPSDE(spsde_v, spsde_f1, spsde_f2, spsde_B, spsde_G1, spsde_G2, TestNoise(); periodicity=NullPeriodicity())

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


    psde = SPSDE(spsde_v, spsde_f1, spsde_f2, spsde_B, spsde_G1, spsde_G2, TestNoise(), NullInvariants(), sde_param_types, NullPeriodicity())

    funcs = functions(psde)

    @test_nowarn funcs.v(zero(q₀), t₀, q₀, p₀, sde_params)
    @test_nowarn funcs.f1(zero(p₀), t₀, q₀, p₀, sde_params)
    @test_nowarn funcs.f2(zero(p₀), t₀, q₀, p₀, sde_params)

    funcs = functions(psde, sde_params)

    @test_nowarn funcs.v(zero(q₀), t₀, q₀, p₀)
    @test_nowarn funcs.f1(zero(p₀), t₀, q₀, p₀, )
    @test_nowarn funcs.f2(zero(p₀), t₀, q₀, p₀, )


    # Test for periodicity
    psdep = SPSDE(spsde_v, spsde_f1, spsde_f2, spsde_B, spsde_G1, spsde_G2, TestNoise(); periodicity=([0.0,],[2π]))

    @test periodicity(psdep) == (Float64[0],Float64[2π])
    @test getperiodicity(psdep) == BitArray([true])
    @test hasperiodicity(psdep) == true

    psdep = SPSDE(spsde_v, spsde_f1, spsde_f2, spsde_B, spsde_G1, spsde_G2, TestNoise(); periodicity=([-Inf],[+Inf]))

    @test periodicity(psdep) == NullPeriodicity()
    @test ismissing(getperiodicity(psdep))
    @test hasperiodicity(psdep) == false

end
