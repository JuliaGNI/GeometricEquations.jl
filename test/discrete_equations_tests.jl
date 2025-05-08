
using GeometricEquations
using Test

include("functions.jl")
include("initial_conditions.jl")


@testset "$(rpad("Discrete Euler-Lagrange Equations (DELE)",80))" begin

    dele  = DELE(dele_eqs..., NullInvariants(), NullParameters(), NullPeriodicity())
    # dele  = DELE(dele_eqs..., dele_igs..., NullInvariants(), NullParameters(), NullPeriodicity())
    dele1 = DELE(dele_eqs...)
    dele2 = DELE(dele_eqs...; invariants=NullInvariants())
    dele3 = DELE(dele_eqs...; parameters=NullParameters())
    dele4 = DELE(dele_eqs...; periodicity=NullPeriodicity())

    @test dele == dele1
    @test dele == dele2
    @test dele == dele3
    @test dele == dele4

    @test hash(dele) == hash(dele1)
    @test hash(dele) == hash(dele2)
    @test hash(dele) == hash(dele3)
    @test hash(dele) == hash(dele4)

    @test functions(dele) == NamedTuple{(:Ld,:D1Ld,:D2Ld)}(dele_eqs)
    @test solutions(dele) == NamedTuple()
    @test initialguess(dele) == NamedTuple()

    @test parameters(dele) == NullParameters()
    @test invariants(dele) == NullInvariants()
    @test periodicity(dele) == NullPeriodicity()

    @test hasvectorfield(dele) == true
    @test hassolution(dele) == false
    @test hasinitialguess(dele) == false
    @test hasprimary(dele) == false
    @test hassecondary(dele) == false

    @test hasinvariants(dele) == false
    @test hasparameters(dele) == false
    @test hasperiodicity(dele) == false

    @test hashamiltonian(dele) == false
    @test haslagrangian(dele) == false

    funcs = functions(dele)

    @test_nowarn funcs.Ld(t₀, t₁, x₀, x₁, NullParameters())
    @test_nowarn funcs.D1Ld(zero(x₀), t₀, t₁, x₀, x₁, NullParameters())
    @test_nowarn funcs.D2Ld(zero(x₀), t₀, t₁, x₀, x₁, NullParameters())

    # funcs = functions(dele, NullParameters())

    # @test_nowarn funcs.Ld(t₀, t₁, x₀, x₁)
    # @test_nowarn funcs.D1Ld(zero(x₀), t₀, t₁, x₀, x₁)
    # @test_nowarn funcs.D2Ld(zero(x₀), t₀, t₁, x₀, x₁)

end
