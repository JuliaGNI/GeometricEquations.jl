using GeometricEquations
using Test

include("functions.jl")
include("initial_conditions.jl")

_copy(x, n) = [x for _ in 1:n]


@testset "$(rpad("Geometric Ensemble",80))" begin

    ics_tpl = [(q=StateVariable(x₀),), (q=StateVariable(rand(length(x₀))),)]
    ics_sva = [StateVariable(x₀), StateVariable(rand(length(x₀)))]
    ics_arr = [x₀, rand(length(x₀))]
    
    @test_nowarn EnsembleProblem(ODE(ode_v), (t₀,t₁), Δt, ics_tpl)
    @test_nowarn EnsembleProblem(ODE(ode_v), (t₀,t₁), Δt, ics_tpl, nothing)
    @test_nowarn EnsembleProblem(ODE(ode_v), (t₀,t₁), Δt, ics_tpl, NullParameters())
    @test_nowarn EnsembleProblem(ODE(ode_v), (t₀,t₁), Δt, ics_tpl; parameters=NullParameters())


    ode = ODE(ode_v, parameters = ode_param_types)
    
    @test_nowarn EnsembleProblem(ode, (t₀,t₁), Δt, ics_tpl, ode_params)
    @test_nowarn EnsembleProblem(ode, (t₀,t₁), Δt, ics_tpl; parameters = ode_params)

    ens = EnsembleProblem(ode, (t₀,t₁), Δt, ics_tpl, ode_params)

    @test typeof(ens) <: EnsembleProblem
    @test typeof(ens) <: ODEEnsemble
    @test typeof(ens).parameters[1] == ODE
    @test typeof(ens).parameters[5] == typeof(ode)

    @test datatype(ens) == eltype(x₀)
    @test timetype(ens) == typeof(t₀)
    @test arrtype(ens) == typeof(StateVariable(x₀))
    @test equtype(ens) == ODE

    @test tspan(ens) == (t₀,t₁)
    @test tbegin(ens) == t₀
    @test tend(ens) == t₁
    @test tstep(ens) == Δt

    @test equation(ens) == ode
    @test functions(ens) == functions(ode)
    @test solutions(ens) == solutions(ode)
    @test parameters(ens) == [ode_params, ode_params]
    @test initial_conditions(ens) == ics_tpl
    @test nsamples(ens) == length(ens) == 2

    @test problem(ens, 2) == ens[2] == ens[CartesianIndex(2)] == ens[CartesianIndex(2,1)]

    probs = (
        EquationProblem(ode, (t₀,t₁), Δt, ics_tpl[1], ode_params),
        EquationProblem(ode, (t₀,t₁), Δt, ics_tpl[2], ode_params),
    )

    for prob in ens
        @test prob ∈ probs
    end


    ode = ODE(ode_v, parameters = ode_param_types)
    params = [(α=1,), (α=2,)]

    @test_nowarn EnsembleProblem(ode, (t₀,t₁), Δt, ode_ics, params)
    @test_nowarn EnsembleProblem(ode, (t₀,t₁), Δt, ode_ics; parameters = params)

    ens = EnsembleProblem(ode, (t₀,t₁), Δt, ode_ics, params)

    @test equation(ens) == ode
    @test functions(ens) == functions(ode)
    @test solutions(ens) == solutions(ode)
    @test parameters(ens) == params

    @test initial_conditions(ens) == [ode_ics, ode_ics]
    @test nsamples(ens) == length(ens) == 2

    @test problem(ens, 2) == ens[2] == ens[CartesianIndex(2)] == ens[CartesianIndex(2,1)]

    probs = (
        EquationProblem(ode, (t₀,t₁), Δt, ode_ics, params[1]),
        EquationProblem(ode, (t₀,t₁), Δt, ode_ics, params[2]),
    )

    for prob in ens
        @test prob ∈ probs
    end


    ode = ODE(ode_v, parameters = ode_param_types)
    params = [(α=1,), (α=2,)]

    @test_nowarn EnsembleProblem(ode, (t₀,t₁), Δt, ics_tpl, params)
    @test_nowarn EnsembleProblem(ode, (t₀,t₁), Δt, ics_tpl; parameters = params)

    ens = EnsembleProblem(ode, (t₀,t₁), Δt, ics_tpl, params)

    @test equation(ens) == ode
    @test functions(ens) == functions(ode)
    @test solutions(ens) == solutions(ode)
    @test parameters(ens) == params

    @test initial_conditions(ens) == ics_tpl
    @test nsamples(ens) == length(ens) == 2

    probs = (
        EquationProblem(ode, (t₀,t₁), Δt, ics_tpl[1], params[1]),
        EquationProblem(ode, (t₀,t₁), Δt, ics_tpl[2], params[2]),
    )

    for prob in ens
        @test prob ∈ probs
    end

end


@testset "$(rpad("ODE Ensemble",80))" begin

    @test_nowarn ODEEnsemble(ode_eqs..., (t₀,t₁), Δt, _copy(ode_ics, 3))
    @test_nowarn ODEEnsemble(ode_eqs..., (t₀,t₁), Δt, _copy(ode_ics, 3); parameters = _copy((α=1,), 3))



    ics_arr = [x₀, rand(length(x₀))]
    ics_sva = [StateVariable(ics_arr[1]), StateVariable(ics_arr[2])]
    ics_tpl = [(q=StateVariable(ics_arr[1]),), (q=StateVariable(ics_arr[2]),)]
    
    ens  = EnsembleProblem(ODE(ode_v), (t₀,t₁), Δt, ics_tpl)

    ens1 = ODEEnsemble(ode_eqs..., (t₀,t₁), Δt, ics_tpl)
    ens2 = ODEEnsemble(ode_eqs..., (t₀,t₁), Δt, ics_tpl; parameters = NullParameters())
    ens3 = ODEEnsemble(ode_eqs..., (t₀,t₁), Δt, ics_sva)
    ens4 = ODEEnsemble(ode_eqs..., (t₀,t₁), Δt, ics_sva; parameters = NullParameters())
    ens5 = ODEEnsemble(ode_eqs..., (t₀,t₁), Δt, ics_arr)
    ens6 = ODEEnsemble(ode_eqs..., (t₀,t₁), Δt, ics_arr; parameters = NullParameters())

    @test ens1 == ens
    @test ens2 == ens
    @test ens3 == ens
    @test ens4 == ens
    @test ens5 == ens
    @test ens6 == ens

end


@testset "$(rpad("SODE Ensemble",80))" begin

    @test_nowarn SODEEnsemble(sode_eqs, (t₀,t₁), Δt, _copy(ode_ics, 3))
    @test_nowarn SODEEnsemble(sode_eqs, (t₀,t₁), Δt, _copy(ode_ics, 3); parameters = _copy((α=1,), 3))
    @test_nowarn SODEEnsemble(sode_eqs, sode_sols, (t₀,t₁), Δt, _copy(ode_ics, 3))
    @test_nowarn SODEEnsemble(sode_eqs, sode_sols, (t₀,t₁), Δt, _copy(ode_ics, 3); parameters = _copy((α=1,), 3))

end


@testset "$(rpad("PODE Ensemble",80))" begin

    @test_nowarn PODEEnsemble(pode_eqs..., (t₀,t₁), Δt, _copy(pode_ics, 3))
    @test_nowarn PODEEnsemble(pode_eqs..., (t₀,t₁), Δt, _copy(pode_ics, 3); parameters = _copy((α=1,), 3))

end


@testset "$(rpad("IODE Ensemble",80))" begin

    @test_nowarn IODEEnsemble(iode_eqs..., (t₀,t₁), Δt, _copy(iode_ics, 3))
    @test_nowarn IODEEnsemble(iode_eqs..., (t₀,t₁), Δt, _copy(iode_ics, 3); parameters = _copy((α=1,), 3))

end


@testset "$(rpad("HODE Ensemble",80))" begin

    @test_nowarn HODEEnsemble(hode_eqs..., (t₀,t₁), Δt, _copy(hode_ics, 3))
    @test_nowarn HODEEnsemble(hode_eqs..., (t₀,t₁), Δt, _copy(hode_ics, 3); parameters = _copy((α=1,), 3))

end


@testset "$(rpad("LODE Ensemble",80))" begin

    @test_nowarn LODEEnsemble(lode_eqs..., (t₀,t₁), Δt, _copy(lode_ics, 3))
    @test_nowarn LODEEnsemble(lode_eqs..., (t₀,t₁), Δt, _copy(lode_ics, 3); parameters = _copy((α=1,), 3))

end


@testset "$(rpad("DAE Ensemble",80))" begin

    @test_nowarn DAEEnsemble(dae_eqs..., (t₀,t₁), Δt, _copy(dae_ics, 3))
    @test_nowarn DAEEnsemble(dae_eqs..., (t₀,t₁), Δt, _copy(dae_ics, 3); parameters = _copy((α=1,), 3))
    @test_nowarn DAEEnsemble(dae_eqs_full..., (t₀,t₁), Δt, _copy(dae_ics_full, 3))
    @test_nowarn DAEEnsemble(dae_eqs_full..., (t₀,t₁), Δt, _copy(dae_ics_full, 3); parameters = _copy((α=1,), 3))

end


@testset "$(rpad("PDAE Ensemble",80))" begin

    @test_nowarn PDAEEnsemble(pdae_eqs..., (t₀,t₁), Δt, _copy(pdae_ics, 3))
    @test_nowarn PDAEEnsemble(pdae_eqs..., (t₀,t₁), Δt, _copy(pdae_ics, 3); parameters = _copy((α=1,), 3))
    @test_nowarn PDAEEnsemble(pdae_eqs_full..., (t₀,t₁), Δt, _copy(pdae_ics_full, 3))
    @test_nowarn PDAEEnsemble(pdae_eqs_full..., (t₀,t₁), Δt, _copy(pdae_ics_full, 3); parameters = _copy((α=1,), 3))

end


@testset "$(rpad("IDAE Ensemble",80))" begin

    @test_nowarn IDAEEnsemble(idae_eqs..., (t₀,t₁), Δt, _copy(idae_ics, 3))
    @test_nowarn IDAEEnsemble(idae_eqs..., (t₀,t₁), Δt, _copy(idae_ics, 3); parameters = _copy((α=1,), 3))
    @test_nowarn IDAEEnsemble(idae_eqs_full..., (t₀,t₁), Δt, _copy(idae_ics_full, 3))
    @test_nowarn IDAEEnsemble(idae_eqs_full..., (t₀,t₁), Δt, _copy(idae_ics_full, 3); parameters = _copy((α=1,), 3))

end


@testset "$(rpad("HDAE Ensemble",80))" begin

    @test_nowarn HDAEEnsemble(hdae_eqs..., (t₀,t₁), Δt, _copy(hdae_ics, 3))
    @test_nowarn HDAEEnsemble(hdae_eqs..., (t₀,t₁), Δt, _copy(hdae_ics, 3); parameters = _copy((α=1,), 3))
    @test_nowarn HDAEEnsemble(hdae_eqs_full..., (t₀,t₁), Δt, _copy(hdae_ics_full, 3))
    @test_nowarn HDAEEnsemble(hdae_eqs_full..., (t₀,t₁), Δt, _copy(hdae_ics_full, 3); parameters = _copy((α=1,), 3))

end


@testset "$(rpad("LDAE Ensemble",80))" begin

    @test_nowarn LDAEEnsemble(ldae_eqs..., (t₀,t₁), Δt, _copy(ldae_ics, 3))
    @test_nowarn LDAEEnsemble(ldae_eqs..., (t₀,t₁), Δt, _copy(ldae_ics, 3); parameters = _copy((α=1,), 3))
    @test_nowarn LDAEEnsemble(ldae_eqs_full..., (t₀,t₁), Δt, _copy(ldae_ics_full, 3))
    @test_nowarn LDAEEnsemble(ldae_eqs_full..., (t₀,t₁), Δt, _copy(ldae_ics_full, 3); parameters = _copy((α=1,), 3))

end
