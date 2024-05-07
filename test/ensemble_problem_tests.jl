using GeometricEquations
using Test

include("functions.jl")
include("initial_conditions.jl")

_copy(x, n) = [x for _ in 1:n]

const size_one_warning = (:warn, "You created an EnsembleProblem with a single initial condition and a single set of parameters. You probably want to create a GeometricProblem instead.")


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

    @test_logs size_one_warning ODEEnsemble(ode_eqs..., (t₀,t₁), Δt, ode_ics)
    @test_logs size_one_warning ODEEnsemble(ode_eqs..., (t₀,t₁), Δt, ode_ics...)
    @test_logs size_one_warning ODEEnsemble(ode_eqs..., (t₀,t₁), Δt, ode_ics_raw...)
    @test_nowarn ODEEnsemble(ode_eqs..., (t₀,t₁), Δt, ode_ics_tpl)
    @test_nowarn ODEEnsemble(ode_eqs..., (t₀,t₁), Δt, ode_ics_sva...)

    @test_logs size_one_warning ODEEnsemble(ode_eqs..., (t₀,t₁), Δt, ode_ics; parameters = ode_params)
    @test_logs size_one_warning ODEEnsemble(ode_eqs..., (t₀,t₁), Δt, ode_ics...; parameters = ode_params)
    @test_logs size_one_warning ODEEnsemble(ode_eqs..., (t₀,t₁), Δt, ode_ics_raw...; parameters = ode_params)
    @test_nowarn ODEEnsemble(ode_eqs..., (t₀,t₁), Δt, ode_ics_tpl; parameters = ode_params)
    @test_nowarn ODEEnsemble(ode_eqs..., (t₀,t₁), Δt, ode_ics_sva...; parameters = ode_params)

    @test_nowarn ODEEnsemble(ode_eqs..., (t₀,t₁), Δt, ode_ics; parameters = _copy(ode_params, 3))
    @test_nowarn ODEEnsemble(ode_eqs..., (t₀,t₁), Δt, ode_ics...; parameters = _copy(ode_params, 3))
    @test_nowarn ODEEnsemble(ode_eqs..., (t₀,t₁), Δt, ode_ics_raw...; parameters = _copy(ode_params, 3))
    @test_nowarn ODEEnsemble(ode_eqs..., (t₀,t₁), Δt, ode_ics_tpl; parameters = _copy(ode_params, 2))
    @test_nowarn ODEEnsemble(ode_eqs..., (t₀,t₁), Δt, ode_ics_sva...; parameters = _copy(ode_params, 2))
    @test_throws AssertionError ODEEnsemble(ode_eqs..., (t₀,t₁), Δt, ode_ics_tpl; parameters = _copy(ode_params, 3))
    @test_throws AssertionError ODEEnsemble(ode_eqs..., (t₀,t₁), Δt, ode_ics_sva...; parameters = _copy(ode_params, 3))

    ens = EnsembleProblem(ODE(ode_eqs...), (t₀,t₁), Δt, ode_ics_tpl)

    @test ens == ODEEnsemble(ode_eqs..., (t₀,t₁), Δt, ode_ics_tpl)
    @test ens == ODEEnsemble(ode_eqs..., (t₀,t₁), Δt, ode_ics_tpl; parameters = NullParameters())
    @test ens == ODEEnsemble(ode_eqs..., (t₀,t₁), Δt, ode_ics_sva...)
    @test ens == ODEEnsemble(ode_eqs..., (t₀,t₁), Δt, ode_ics_sva...; parameters = NullParameters())
    @test ens == ODEEnsemble(ode_eqs..., (t₀,t₁), Δt, ode_ics_arr...)
    @test ens == ODEEnsemble(ode_eqs..., (t₀,t₁), Δt, ode_ics_arr...; parameters = NullParameters())

end


@testset "$(rpad("SODE Ensemble",80))" begin

    @test_logs size_one_warning SODEEnsemble(sode_eqs, (t₀,t₁), Δt, ode_ics)
    @test_logs size_one_warning SODEEnsemble(sode_eqs, (t₀,t₁), Δt, ode_ics...)
    @test_logs size_one_warning SODEEnsemble(sode_eqs, (t₀,t₁), Δt, ode_ics_raw...)
    @test_nowarn SODEEnsemble(sode_eqs, (t₀,t₁), Δt, ode_ics_tpl)
    @test_nowarn SODEEnsemble(sode_eqs, (t₀,t₁), Δt, ode_ics_sva...)
    @test_nowarn SODEEnsemble(sode_eqs, (t₀,t₁), Δt, ode_ics_tpl; parameters = _copy(ode_params, 2))
    @test_nowarn SODEEnsemble(sode_eqs, (t₀,t₁), Δt, ode_ics_sva...; parameters = _copy(ode_params, 2))
    @test_throws AssertionError SODEEnsemble(sode_eqs, (t₀,t₁), Δt, ode_ics_sva...; parameters = _copy(ode_params, 3))

    @test_logs size_one_warning SODEEnsemble(sode_eqs, sode_sols, (t₀,t₁), Δt, ode_ics)
    @test_logs size_one_warning SODEEnsemble(sode_eqs, sode_sols, (t₀,t₁), Δt, ode_ics...)
    @test_logs size_one_warning SODEEnsemble(sode_eqs, sode_sols, (t₀,t₁), Δt, ode_ics_raw...)
    @test_nowarn SODEEnsemble(sode_eqs, sode_sols, (t₀,t₁), Δt, ode_ics_tpl)
    @test_nowarn SODEEnsemble(sode_eqs, sode_sols, (t₀,t₁), Δt, ode_ics_sva...)
    @test_nowarn SODEEnsemble(sode_eqs, sode_sols, (t₀,t₁), Δt, ode_ics_tpl; parameters = _copy(ode_params, 2))
    @test_nowarn SODEEnsemble(sode_eqs, sode_sols, (t₀,t₁), Δt, ode_ics_sva...; parameters = _copy(ode_params, 2))
    @test_throws AssertionError SODEEnsemble(sode_eqs, sode_sols, (t₀,t₁), Δt, ode_ics_sva...; parameters = _copy(ode_params, 3))
    
    ens = EnsembleProblem(SODE(sode_eqs), (t₀,t₁), Δt, ode_ics_tpl)

    @test ens == SODEEnsemble(sode_eqs, (t₀,t₁), Δt, ode_ics_tpl)
    @test ens == SODEEnsemble(sode_eqs, (t₀,t₁), Δt, ode_ics_tpl; parameters = NullParameters())
    @test ens == SODEEnsemble(sode_eqs, (t₀,t₁), Δt, ode_ics_sva...)
    @test ens == SODEEnsemble(sode_eqs, (t₀,t₁), Δt, ode_ics_sva...; parameters = NullParameters())
    @test ens == SODEEnsemble(sode_eqs, (t₀,t₁), Δt, ode_ics_arr...)
    @test ens == SODEEnsemble(sode_eqs, (t₀,t₁), Δt, ode_ics_arr...; parameters = NullParameters())

    ens = EnsembleProblem(SODE(sode_eqs, sode_sols), (t₀,t₁), Δt, ode_ics_tpl)

    @test ens == SODEEnsemble(sode_eqs, sode_sols, (t₀,t₁), Δt, ode_ics_tpl)
    @test ens == SODEEnsemble(sode_eqs, sode_sols, (t₀,t₁), Δt, ode_ics_tpl; parameters = NullParameters())
    @test ens == SODEEnsemble(sode_eqs, sode_sols, (t₀,t₁), Δt, ode_ics_sva...)
    @test ens == SODEEnsemble(sode_eqs, sode_sols, (t₀,t₁), Δt, ode_ics_sva...; parameters = NullParameters())
    @test ens == SODEEnsemble(sode_eqs, sode_sols, (t₀,t₁), Δt, ode_ics_arr...)
    @test ens == SODEEnsemble(sode_eqs, sode_sols, (t₀,t₁), Δt, ode_ics_arr...; parameters = NullParameters())

end


@testset "$(rpad("PODE Ensemble",80))" begin

    @test_logs size_one_warning PODEEnsemble(pode_eqs..., (t₀,t₁), Δt, pode_ics)
    @test_logs size_one_warning PODEEnsemble(pode_eqs..., (t₀,t₁), Δt, pode_ics...)
    @test_logs size_one_warning PODEEnsemble(pode_eqs..., (t₀,t₁), Δt, pode_ics_raw...)
    @test_nowarn PODEEnsemble(pode_eqs..., (t₀,t₁), Δt, pode_ics_tpl)
    @test_nowarn PODEEnsemble(pode_eqs..., (t₀,t₁), Δt, pode_ics_sva...)

    @test_logs size_one_warning PODEEnsemble(pode_eqs..., (t₀,t₁), Δt, pode_ics; parameters = ode_params)
    @test_logs size_one_warning PODEEnsemble(pode_eqs..., (t₀,t₁), Δt, pode_ics...; parameters = ode_params)
    @test_logs size_one_warning PODEEnsemble(pode_eqs..., (t₀,t₁), Δt, pode_ics_raw...; parameters = ode_params)
    @test_nowarn PODEEnsemble(pode_eqs..., (t₀,t₁), Δt, pode_ics_tpl; parameters = ode_params)
    @test_nowarn PODEEnsemble(pode_eqs..., (t₀,t₁), Δt, pode_ics_sva...; parameters = ode_params)

    @test_nowarn PODEEnsemble(pode_eqs..., (t₀,t₁), Δt, pode_ics; parameters = _copy(ode_params, 3))
    @test_nowarn PODEEnsemble(pode_eqs..., (t₀,t₁), Δt, pode_ics...; parameters = _copy(ode_params, 3))
    @test_nowarn PODEEnsemble(pode_eqs..., (t₀,t₁), Δt, pode_ics_raw...; parameters = _copy(ode_params, 3))
    @test_nowarn PODEEnsemble(pode_eqs..., (t₀,t₁), Δt, pode_ics_tpl; parameters = _copy(ode_params, 2))
    @test_nowarn PODEEnsemble(pode_eqs..., (t₀,t₁), Δt, pode_ics_sva...; parameters = _copy(ode_params, 2))
    @test_throws AssertionError PODEEnsemble(pode_eqs..., (t₀,t₁), Δt, pode_ics_tpl; parameters = _copy(ode_params, 3))
    @test_throws AssertionError PODEEnsemble(pode_eqs..., (t₀,t₁), Δt, pode_ics_sva...; parameters = _copy(ode_params, 3))

    ens = EnsembleProblem(PODE(pode_eqs...), (t₀,t₁), Δt, pode_ics_tpl)

    @test ens == PODEEnsemble(pode_eqs..., (t₀,t₁), Δt, pode_ics_tpl)
    @test ens == PODEEnsemble(pode_eqs..., (t₀,t₁), Δt, pode_ics_tpl; parameters = NullParameters())
    @test ens == PODEEnsemble(pode_eqs..., (t₀,t₁), Δt, pode_ics_sva...)
    @test ens == PODEEnsemble(pode_eqs..., (t₀,t₁), Δt, pode_ics_sva...; parameters = NullParameters())
    @test ens == PODEEnsemble(pode_eqs..., (t₀,t₁), Δt, pode_ics_arr...)
    @test ens == PODEEnsemble(pode_eqs..., (t₀,t₁), Δt, pode_ics_arr...; parameters = NullParameters())

end


@testset "$(rpad("IODE Ensemble",80))" begin

    @test_logs size_one_warning IODEEnsemble(iode_eqs..., (t₀,t₁), Δt, iode_ics)
    # @test_logs size_one_warning IODEEnsemble(iode_eqs..., (t₀,t₁), Δt, iode_ics...)
    # @test_logs size_one_warning IODEEnsemble(iode_eqs..., (t₀,t₁), Δt, iode_ics_raw...)
    @test_nowarn IODEEnsemble(iode_eqs..., (t₀,t₁), Δt, iode_ics_tpl)
    @test_nowarn IODEEnsemble(iode_eqs..., (t₀,t₁), Δt, iode_ics_sva...)

    @test_logs size_one_warning IODEEnsemble(iode_eqs..., (t₀,t₁), Δt, iode_ics; parameters = ode_params)
    # @test_logs size_one_warning IODEEnsemble(iode_eqs..., (t₀,t₁), Δt, iode_ics...; parameters = ode_params)
    # @test_logs size_one_warning IODEEnsemble(iode_eqs..., (t₀,t₁), Δt, iode_ics_raw...; parameters = ode_params)
    @test_nowarn IODEEnsemble(iode_eqs..., (t₀,t₁), Δt, iode_ics_tpl; parameters = ode_params)
    @test_nowarn IODEEnsemble(iode_eqs..., (t₀,t₁), Δt, iode_ics_sva...; parameters = ode_params)

    @test_nowarn IODEEnsemble(iode_eqs..., (t₀,t₁), Δt, iode_ics; parameters = _copy(ode_params, 3))
    # @test_nowarn IODEEnsemble(iode_eqs..., (t₀,t₁), Δt, iode_ics...; parameters = _copy(ode_params, 3))
    # @test_nowarn IODEEnsemble(iode_eqs..., (t₀,t₁), Δt, iode_ics_raw...; parameters = _copy(ode_params, 3))
    @test_nowarn IODEEnsemble(iode_eqs..., (t₀,t₁), Δt, iode_ics_tpl; parameters = _copy(ode_params, 2))
    @test_nowarn IODEEnsemble(iode_eqs..., (t₀,t₁), Δt, iode_ics_sva...; parameters = _copy(ode_params, 2))
    @test_throws AssertionError IODEEnsemble(iode_eqs..., (t₀,t₁), Δt, iode_ics_tpl; parameters = _copy(ode_params, 3))
    @test_throws AssertionError IODEEnsemble(iode_eqs..., (t₀,t₁), Δt, iode_ics_sva...; parameters = _copy(ode_params, 3))

    ens = EnsembleProblem(IODE(iode_eqs...), (t₀,t₁), Δt, iode_ics_tpl)

    @test ens == IODEEnsemble(iode_eqs..., (t₀,t₁), Δt, iode_ics_tpl)
    @test ens == IODEEnsemble(iode_eqs..., (t₀,t₁), Δt, iode_ics_tpl; parameters = NullParameters())
    @test ens == IODEEnsemble(iode_eqs..., (t₀,t₁), Δt, q_sva, p_sva, λ_sva)
    @test ens == IODEEnsemble(iode_eqs..., (t₀,t₁), Δt, q_sva, p_sva, λ_sva; parameters = NullParameters())
    @test ens == IODEEnsemble(iode_eqs..., (t₀,t₁), Δt, q_arr, p_arr, λ_arr)
    @test ens == IODEEnsemble(iode_eqs..., (t₀,t₁), Δt, q_arr, p_arr, λ_arr; parameters = NullParameters())
    @test ens == IODEEnsemble(iode_eqs..., (t₀,t₁), Δt, q_sva, p_sva)
    @test ens == IODEEnsemble(iode_eqs..., (t₀,t₁), Δt, q_sva, p_sva; parameters = NullParameters())
    @test ens == IODEEnsemble(iode_eqs..., (t₀,t₁), Δt, q_arr, p_arr)
    @test ens == IODEEnsemble(iode_eqs..., (t₀,t₁), Δt, q_arr, p_arr; parameters = NullParameters())

    ens = EnsembleProblem(IODE(iode_eqs_default_g...), (t₀,t₁), Δt, iode_ics_tpl)

    @test ens == IODEEnsemble(iode_eqs_without_g..., (t₀,t₁), Δt, iode_ics_tpl)
    @test ens == IODEEnsemble(iode_eqs_without_g..., (t₀,t₁), Δt, iode_ics_tpl; parameters = NullParameters())
    @test ens == IODEEnsemble(iode_eqs_without_g..., (t₀,t₁), Δt, q_sva, p_sva, λ_sva)
    @test ens == IODEEnsemble(iode_eqs_without_g..., (t₀,t₁), Δt, q_sva, p_sva, λ_sva; parameters = NullParameters())
    @test ens == IODEEnsemble(iode_eqs_without_g..., (t₀,t₁), Δt, q_arr, p_arr, λ_arr)
    @test ens == IODEEnsemble(iode_eqs_without_g..., (t₀,t₁), Δt, q_arr, p_arr, λ_arr; parameters = NullParameters())
    @test ens == IODEEnsemble(iode_eqs_without_g..., (t₀,t₁), Δt, q_sva, p_sva)
    @test ens == IODEEnsemble(iode_eqs_without_g..., (t₀,t₁), Δt, q_sva, p_sva; parameters = NullParameters())
    @test ens == IODEEnsemble(iode_eqs_without_g..., (t₀,t₁), Δt, q_arr, p_arr)
    @test ens == IODEEnsemble(iode_eqs_without_g..., (t₀,t₁), Δt, q_arr, p_arr; parameters = NullParameters())

end


@testset "$(rpad("HODE Ensemble",80))" begin

    @test_logs size_one_warning HODEEnsemble(hode_eqs..., (t₀,t₁), Δt, hode_ics)
    @test_logs size_one_warning HODEEnsemble(hode_eqs..., (t₀,t₁), Δt, hode_ics...)
    @test_logs size_one_warning HODEEnsemble(hode_eqs..., (t₀,t₁), Δt, hode_ics_raw...)
    @test_nowarn HODEEnsemble(hode_eqs..., (t₀,t₁), Δt, hode_ics_tpl)
    @test_nowarn HODEEnsemble(hode_eqs..., (t₀,t₁), Δt, hode_ics_sva...)

    @test_logs size_one_warning HODEEnsemble(hode_eqs..., (t₀,t₁), Δt, hode_ics; parameters = ode_params)
    @test_logs size_one_warning HODEEnsemble(hode_eqs..., (t₀,t₁), Δt, hode_ics...; parameters = ode_params)
    @test_logs size_one_warning HODEEnsemble(hode_eqs..., (t₀,t₁), Δt, hode_ics_raw...; parameters = ode_params)
    @test_nowarn HODEEnsemble(hode_eqs..., (t₀,t₁), Δt, hode_ics_tpl; parameters = ode_params)
    @test_nowarn HODEEnsemble(hode_eqs..., (t₀,t₁), Δt, hode_ics_sva...; parameters = ode_params)

    @test_nowarn HODEEnsemble(hode_eqs..., (t₀,t₁), Δt, hode_ics; parameters = _copy(ode_params, 3))
    @test_nowarn HODEEnsemble(hode_eqs..., (t₀,t₁), Δt, hode_ics...; parameters = _copy(ode_params, 3))
    @test_nowarn HODEEnsemble(hode_eqs..., (t₀,t₁), Δt, hode_ics_raw...; parameters = _copy(ode_params, 3))
    @test_nowarn HODEEnsemble(hode_eqs..., (t₀,t₁), Δt, hode_ics_tpl; parameters = _copy(ode_params, 2))
    @test_nowarn HODEEnsemble(hode_eqs..., (t₀,t₁), Δt, hode_ics_sva...; parameters = _copy(ode_params, 2))
    @test_throws AssertionError HODEEnsemble(hode_eqs..., (t₀,t₁), Δt, hode_ics_tpl; parameters = _copy(ode_params, 3))
    @test_throws AssertionError HODEEnsemble(hode_eqs..., (t₀,t₁), Δt, hode_ics_sva...; parameters = _copy(ode_params, 3))

    ens = EnsembleProblem(HODE(hode_eqs...), (t₀,t₁), Δt, hode_ics_tpl)

    @test ens == HODEEnsemble(hode_eqs..., (t₀,t₁), Δt, hode_ics_tpl)
    @test ens == HODEEnsemble(hode_eqs..., (t₀,t₁), Δt, hode_ics_tpl; parameters = NullParameters())
    @test ens == HODEEnsemble(hode_eqs..., (t₀,t₁), Δt, q_sva, p_sva)
    @test ens == HODEEnsemble(hode_eqs..., (t₀,t₁), Δt, q_sva, p_sva; parameters = NullParameters())
    @test ens == HODEEnsemble(hode_eqs..., (t₀,t₁), Δt, q_arr, p_arr)
    @test ens == HODEEnsemble(hode_eqs..., (t₀,t₁), Δt, q_arr, p_arr; parameters = NullParameters())

end


@testset "$(rpad("LODE Ensemble",80))" begin

    @test_logs size_one_warning LODEEnsemble(lode_eqs..., (t₀,t₁), Δt, lode_ics)
    # @test_logs size_one_warning LODEEnsemble(lode_eqs..., (t₀,t₁), Δt, lode_ics...)
    # @test_logs size_one_warning LODEEnsemble(lode_eqs..., (t₀,t₁), Δt, lode_ics_raw...)
    @test_nowarn LODEEnsemble(lode_eqs..., (t₀,t₁), Δt, lode_ics_tpl)
    @test_nowarn LODEEnsemble(lode_eqs..., (t₀,t₁), Δt, lode_ics_sva...)

    @test_logs size_one_warning LODEEnsemble(lode_eqs..., (t₀,t₁), Δt, lode_ics; parameters = ode_params)
    # @test_logs size_one_warning LODEEnsemble(lode_eqs..., (t₀,t₁), Δt, lode_ics...; parameters = ode_params)
    # @test_logs size_one_warning LODEEnsemble(lode_eqs..., (t₀,t₁), Δt, lode_ics_raw...; parameters = ode_params)
    @test_nowarn LODEEnsemble(lode_eqs..., (t₀,t₁), Δt, lode_ics_tpl; parameters = ode_params)
    @test_nowarn LODEEnsemble(lode_eqs..., (t₀,t₁), Δt, lode_ics_sva...; parameters = ode_params)

    @test_nowarn LODEEnsemble(lode_eqs..., (t₀,t₁), Δt, lode_ics; parameters = _copy(ode_params, 3))
    # @test_nowarn LODEEnsemble(lode_eqs..., (t₀,t₁), Δt, lode_ics...; parameters = _copy(ode_params, 3))
    # @test_nowarn LODEEnsemble(lode_eqs..., (t₀,t₁), Δt, lode_ics_raw...; parameters = _copy(ode_params, 3))
    @test_nowarn LODEEnsemble(lode_eqs..., (t₀,t₁), Δt, lode_ics_tpl; parameters = _copy(ode_params, 2))
    @test_nowarn LODEEnsemble(lode_eqs..., (t₀,t₁), Δt, lode_ics_sva...; parameters = _copy(ode_params, 2))
    @test_throws AssertionError LODEEnsemble(lode_eqs..., (t₀,t₁), Δt, lode_ics_tpl; parameters = _copy(ode_params, 3))
    @test_throws AssertionError LODEEnsemble(lode_eqs..., (t₀,t₁), Δt, lode_ics_sva...; parameters = _copy(ode_params, 3))

    ens = EnsembleProblem(LODE(lode_eqs...), (t₀,t₁), Δt, lode_ics_tpl)

    @test ens == LODEEnsemble(lode_eqs..., (t₀,t₁), Δt, lode_ics_tpl)
    @test ens == LODEEnsemble(lode_eqs..., (t₀,t₁), Δt, lode_ics_tpl; parameters = NullParameters())
    @test ens == LODEEnsemble(lode_eqs..., (t₀,t₁), Δt, q_sva, p_sva, λ_sva)
    @test ens == LODEEnsemble(lode_eqs..., (t₀,t₁), Δt, q_sva, p_sva, λ_sva; parameters = NullParameters())
    @test ens == LODEEnsemble(lode_eqs..., (t₀,t₁), Δt, q_arr, p_arr, λ_arr)
    @test ens == LODEEnsemble(lode_eqs..., (t₀,t₁), Δt, q_arr, p_arr, λ_arr; parameters = NullParameters())
    @test ens == LODEEnsemble(lode_eqs..., (t₀,t₁), Δt, q_sva, p_sva)
    @test ens == LODEEnsemble(lode_eqs..., (t₀,t₁), Δt, q_sva, p_sva; parameters = NullParameters())
    @test ens == LODEEnsemble(lode_eqs..., (t₀,t₁), Δt, q_arr, p_arr)
    @test ens == LODEEnsemble(lode_eqs..., (t₀,t₁), Δt, q_arr, p_arr; parameters = NullParameters())

    ens = EnsembleProblem(LODE(lode_eqs_default_g...), (t₀,t₁), Δt, lode_ics_tpl)

    @test ens == LODEEnsemble(lode_eqs_without_g..., (t₀,t₁), Δt, lode_ics_tpl)
    @test ens == LODEEnsemble(lode_eqs_without_g..., (t₀,t₁), Δt, lode_ics_tpl; parameters = NullParameters())
    @test ens == LODEEnsemble(lode_eqs_without_g..., (t₀,t₁), Δt, q_sva, p_sva, λ_sva)
    @test ens == LODEEnsemble(lode_eqs_without_g..., (t₀,t₁), Δt, q_sva, p_sva, λ_sva; parameters = NullParameters())
    @test ens == LODEEnsemble(lode_eqs_without_g..., (t₀,t₁), Δt, q_arr, p_arr, λ_arr)
    @test ens == LODEEnsemble(lode_eqs_without_g..., (t₀,t₁), Δt, q_arr, p_arr, λ_arr; parameters = NullParameters())
    @test ens == LODEEnsemble(lode_eqs_without_g..., (t₀,t₁), Δt, q_sva, p_sva)
    @test ens == LODEEnsemble(lode_eqs_without_g..., (t₀,t₁), Δt, q_sva, p_sva; parameters = NullParameters())
    @test ens == LODEEnsemble(lode_eqs_without_g..., (t₀,t₁), Δt, q_arr, p_arr)
    @test ens == LODEEnsemble(lode_eqs_without_g..., (t₀,t₁), Δt, q_arr, p_arr; parameters = NullParameters())

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
