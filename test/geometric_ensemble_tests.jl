using GeometricEquations
using Test

include("functions.jl")
include("initial_conditions.jl")


@testset "$(rpad("Geometric Ensemble",80))" begin

    ics = [(q=x₀,), (q=rand(length(x₀)),)]
    
    @test_nowarn GeometricEnsemble(ODE(ode_v), (t₀,t₁), Δt, ics)
    @test_nowarn GeometricEnsemble(ODE(ode_v), (t₀,t₁), Δt, ics, nothing)
    @test_nowarn GeometricEnsemble(ODE(ode_v), (t₀,t₁), Δt, ics, NullParameters())
    @test_nowarn GeometricEnsemble(ODE(ode_v), (t₀,t₁), Δt, ics; parameters=NullParameters())


    ode = ODE(ode_v, parameters = ode_param_types)
    ics = [(q=x₀,), (q=rand(length(x₀)),)]
    params = ode_params
    
    @test_nowarn GeometricEnsemble(ode, (t₀,t₁), Δt, ics, params)
    @test_nowarn GeometricEnsemble(ode, (t₀,t₁), Δt, ics; parameters=params)

    ens = GeometricEnsemble(ode, (t₀,t₁), Δt, ics, params)

    @test typeof(ens) <: GeometricEnsemble
    @test typeof(ens) <: ODEEnsemble
    @test typeof(ens).parameters[1] == ODE
    @test typeof(ens).parameters[5] == typeof(ode)

    @test datatype(ens) == eltype(x₀)
    @test timetype(ens) == typeof(t₀)
    @test arrtype(ens) == typeof(x₀)
    @test equtype(ens) == ODE

    @test tspan(ens) == (t₀,t₁)
    @test tbegin(ens) == t₀
    @test tend(ens) == t₁
    @test tstep(ens) == Δt

    @test equation(ens) == ode
    @test functions(ens) == functions(ode)
    @test solutions(ens) == solutions(ode)
    @test parameters(ens) == [params, params]
    @test initial_conditions(ens) == ics
    @test nsamples(ens) == 2

    probs = (
        GeometricProblem(ode, (t₀,t₁), Δt, ics[1], params),
        GeometricProblem(ode, (t₀,t₁), Δt, ics[2], params),
    )

    for prob in ens
        @test prob ∈ probs
    end


    ode = ODE(ode_v, parameters = ode_param_types)
    ics = (q=x₀,)
    params = [(α=1,), (α=2,)]

    @test_nowarn GeometricEnsemble(ode, (t₀,t₁), Δt, ics, params)
    @test_nowarn GeometricEnsemble(ode, (t₀,t₁), Δt, ics; parameters=params)

    ens = GeometricEnsemble(ode, (t₀,t₁), Δt, ics, params)

    @test equation(ens) == ode
    @test functions(ens) == functions(ode)
    @test solutions(ens) == solutions(ode)
    @test parameters(ens) == params

    @test initial_conditions(ens) == [ics, ics]
    @test nsamples(ens) == 2

    probs = (
        GeometricProblem(ode, (t₀,t₁), Δt, ics, params[1]),
        GeometricProblem(ode, (t₀,t₁), Δt, ics, params[2]),
    )

    for prob in ens
        @test prob ∈ probs
    end


    ode = ODE(ode_v, parameters = ode_param_types)
    ics = [(q=x₀,), (q=rand(length(x₀)),)]
    params = [(α=1,), (α=2,)]

    @test_nowarn GeometricEnsemble(ode, (t₀,t₁), Δt, ics, params)
    @test_nowarn GeometricEnsemble(ode, (t₀,t₁), Δt, ics; parameters=params)

    ens = GeometricEnsemble(ode, (t₀,t₁), Δt, ics, params)

    @test equation(ens) == ode
    @test functions(ens) == functions(ode)
    @test solutions(ens) == solutions(ode)
    @test parameters(ens) == params

    @test initial_conditions(ens) == ics
    @test nsamples(ens) == 2

    probs = (
        GeometricProblem(ode, (t₀,t₁), Δt, ics[1], params[1]),
        GeometricProblem(ode, (t₀,t₁), Δt, ics[2], params[2]),
    )

    for prob in ens
        @test prob ∈ probs
    end

end
