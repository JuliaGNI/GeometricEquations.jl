using GeometricEquations
using Test

include("functions.jl")
include("initial_conditions.jl")

@testset "$(rpad("Conversion between Problem Types",80))" begin
    pode = PODEProblem(pode_eqs..., (t₀, t₁), Δt, pode_ics)
    hode = HODEProblem(hode_eqs..., (t₀, t₁), Δt, hode_ics)
    lode = LODEProblem(lode_eqs..., (t₀, t₁), Δt, lode_ics)

    # PODEProblem -> ODEProblem

    ode_from_pode = convert(ODEProblem, pode)

    @test typeof(ode_from_pode) <: ODEProblem
    @test timespan(ode_from_pode) == timespan(pode)
    @test timestep(ode_from_pode) == timestep(pode)
    @test initial_conditions(ode_from_pode) == (t = t₀, q = vcat(q₀, p₀))

    ẋ₁ = zero(x₀)
    ẋ₂ = zero(x₀)
    ode_v(ẋ₁, t₀, x₀, NullParameters())
    equation(ode_from_pode).v(ẋ₂, t₀, x₀, NullParameters())
    @test ẋ₁ == ẋ₂

    # PODEProblem -> SODEProblem

    sode_from_pode = convert(SODEProblem, pode)

    @test typeof(sode_from_pode) <: SODEProblem
    @test timespan(sode_from_pode) == timespan(pode)
    @test timestep(sode_from_pode) == timestep(pode)
    @test initial_conditions(sode_from_pode) == (t = t₀, q = vcat(q₀, p₀))

    ẋ₁ = zero(x₀)
    ẋ₂ = zero(x₀)
    sode_v1(ẋ₁, t₀, x₀, NullParameters())
    equation(sode_from_pode).v[1](ẋ₂, t₀, x₀, NullParameters())
    @test ẋ₁ == ẋ₂

    ẋ₁ = zero(x₀)
    ẋ₂ = zero(x₀)
    sode_v2(ẋ₁, t₀, x₀, NullParameters())
    equation(sode_from_pode).v[2](ẋ₂, t₀, x₀, NullParameters())
    @test ẋ₁ == ẋ₂

    # HODEProblem -> ODEProblem and SODEProblem

    ode_from_hode = convert(ODEProblem, hode)

    @test typeof(ode_from_hode) <: ODEProblem
    @test timespan(ode_from_hode) == timespan(hode)
    @test timestep(ode_from_hode) == timestep(hode)
    @test initial_conditions(ode_from_hode) == (t = t₀, q = vcat(q₀, p₀))

    ẋ₁ = zero(x₀)
    ẋ₂ = zero(x₀)
    ode_v(ẋ₁, t₀, x₀, NullParameters())
    equation(ode_from_hode).v(ẋ₂, t₀, x₀, NullParameters())
    @test ẋ₁ == ẋ₂

    sode_from_hode = convert(SODEProblem, hode)

    @test typeof(sode_from_hode) <: SODEProblem
    @test timespan(sode_from_hode) == timespan(hode)
    @test timestep(sode_from_hode) == timestep(hode)
    @test initial_conditions(sode_from_hode) == (t = t₀, q = vcat(q₀, p₀))

    ẋ₁ = zero(x₀)
    ẋ₂ = zero(x₀)
    sode_v1(ẋ₁, t₀, x₀, NullParameters())
    equation(sode_from_hode).v[1](ẋ₂, t₀, x₀, NullParameters())
    @test ẋ₁ == ẋ₂

    ẋ₁ = zero(x₀)
    ẋ₂ = zero(x₀)
    sode_v2(ẋ₁, t₀, x₀, NullParameters())
    equation(sode_from_hode).v[2](ẋ₂, t₀, x₀, NullParameters())
    @test ẋ₁ == ẋ₂

    # HODEProblem -> PODEProblem

    pode_from_hode = convert(PODEProblem, hode)

    @test typeof(pode_from_hode) <: PODEProblem
    @test timespan(pode_from_hode) == timespan(hode)
    @test timestep(pode_from_hode) == timestep(hode)
    @test initial_conditions(pode_from_hode) == (t = t₀, q = q₀, p = p₀)

    @test equation(pode_from_hode).v === equation(hode).v
    @test equation(pode_from_hode).f === equation(hode).f

    # LODEProblem -> IODEProblem

    iode_from_lode = convert(IODEProblem, lode)

    @test typeof(iode_from_lode) <: IODEProblem
    @test timespan(iode_from_lode) == timespan(lode)
    @test timestep(iode_from_lode) == timestep(lode)
    @test initial_conditions(iode_from_lode) == (t = t₀, q = q₀, p = p₀, v = v₀)

    @test equation(iode_from_lode).ϑ === equation(lode).ϑ
    @test equation(iode_from_lode).f === equation(lode).f
    @test equation(iode_from_lode).g === equation(lode).g
    @test equation(iode_from_lode).v̄ === equation(lode).v̄
    @test equation(iode_from_lode).f̄ === equation(lode).f̄

    # Periodicity is extended when the state vectors are concatenated: the momentum block
    # is not periodic, which the package encodes as the bounds `(-Inf, +Inf)`.

    podep = PODEProblem(pode_eqs..., (t₀, t₁), Δt, pode_ics;
        periodicity = ([0.0], [2π]))

    @test periodicity(equation(convert(ODEProblem, podep))) ==
          ([0.0, -Inf], [2π, +Inf])
    @test periodicity(equation(convert(SODEProblem, podep))) ==
          ([0.0, -Inf], [2π, +Inf])
end
