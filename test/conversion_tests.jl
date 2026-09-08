using GeometricEquations
using Test

include("functions.jl")
include("initial_conditions.jl")

@testset "$(rpad("Conversion between Problem Types",80))" begin
    pode = PODEProblem(pode_eqs..., (t₀, t₁), Δt, pode_ics)
    hode = HODEProblem(hode_eqs..., (t₀, t₁), Δt, hode_ics)
    lode = LODEProblem(lode_eqs..., (t₀, t₁), Δt, lode_ics)

    # PODEProblem -> ODEProblem

    ode = convert(ODEProblem, pode)

    @test typeof(ode) <: ODEProblem
    @test timespan(ode) == timespan(pode)
    @test timestep(ode) == timestep(pode)
    @test initial_conditions(ode) == (t = t₀, q = vcat(q₀, p₀))

    ẋ₁ = zero(x₀)
    ẋ₂ = zero(x₀)
    ode_v(ẋ₁, t₀, x₀, NullParameters())
    equation(ode).v(ẋ₂, t₀, x₀, NullParameters())
    @test ẋ₁ == ẋ₂

    # PODEProblem -> SODEProblem

    sode = convert(SODEProblem, pode)

    @test typeof(sode) <: SODEProblem
    @test timespan(sode) == timespan(pode)
    @test timestep(sode) == timestep(pode)
    @test initial_conditions(sode) == (t = t₀, q = vcat(q₀, p₀))

    ẋ₁ = zero(x₀)
    ẋ₂ = zero(x₀)
    sode_v1(ẋ₁, t₀, x₀, NullParameters())
    equation(sode).v[1](ẋ₂, t₀, x₀, NullParameters())
    @test ẋ₁ == ẋ₂

    ẋ₁ = zero(x₀)
    ẋ₂ = zero(x₀)
    sode_v2(ẋ₁, t₀, x₀, NullParameters())
    equation(sode).v[2](ẋ₂, t₀, x₀, NullParameters())
    @test ẋ₁ == ẋ₂

    # HODEProblem -> ODEProblem and SODEProblem

    ode = convert(ODEProblem, hode)

    @test typeof(ode) <: ODEProblem
    @test timespan(ode) == timespan(hode)
    @test timestep(ode) == timestep(hode)
    @test initial_conditions(ode) == (t = t₀, q = vcat(q₀, p₀))

    sode = convert(SODEProblem, hode)

    @test typeof(sode) <: SODEProblem
    @test timespan(sode) == timespan(hode)
    @test timestep(sode) == timestep(hode)
    @test initial_conditions(sode) == (t = t₀, q = vcat(q₀, p₀))

    # HODEProblem -> PODEProblem

    pode = convert(PODEProblem, hode)

    @test typeof(pode) <: PODEProblem
    @test timespan(pode) == timespan(hode)
    @test timestep(pode) == timestep(hode)
    @test initial_conditions(pode) == (t = t₀, q = q₀, p = p₀)

    @test equation(pode).v === equation(hode).v
    @test equation(pode).f === equation(hode).f

    # LODEProblem -> IODEProblem

    iode = convert(IODEProblem, lode)

    @test typeof(iode) <: IODEProblem
    @test timespan(iode) == timespan(lode)
    @test timestep(iode) == timestep(lode)
    @test initial_conditions(iode) == (t = t₀, q = q₀, p = p₀, v = v₀)

    @test equation(iode).ϑ === equation(lode).ϑ
    @test equation(iode).f === equation(lode).f
    @test equation(iode).g === equation(lode).g
    @test equation(iode).v̄ === equation(lode).v̄
    @test equation(iode).f̄ === equation(lode).f̄

    # Periodicity is extended when the state vectors are concatenated: the momentum block
    # is not periodic, which `getperiodicity` reads off as `(-Inf, +Inf)`.
    # `extend_periodicity` still assumes periodicity is a single vector and calls
    # `zero` on the `(lower, upper)` tuple, so both conversions currently throw.

    podep = PODEProblem(pode_eqs..., (t₀, t₁), Δt, pode_ics;
        periodicity = ([0.0], [2π]))

    @test_broken periodicity(equation(convert(ODEProblem, podep))) ==
                 ([0.0, -Inf], [2π, +Inf])
    @test_broken periodicity(equation(convert(SODEProblem, podep))) ==
                 ([0.0, -Inf], [2π, +Inf])
end
