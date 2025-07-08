using GeometricEquations
using Test

include("functions.jl")
include("initial_conditions.jl")


@testset "$(rpad("Geometric Problem",80))" begin

    ode  = ODE(ode_v)
    prob = EquationProblem(ode, (t₀,t₁), Δt, ode_ics)

    @test typeof(prob) <: EquationProblem
    @test typeof(prob) <: ODEProblem
    @test typeof(prob).parameters[1] == ODE
    @test typeof(prob).parameters[5] == typeof(ode)

    @test datatype(prob) == eltype(x₀)
    @test timetype(prob) == typeof(t₀)
    @test arrtype(prob) == StateVariable{eltype(x₀), ndims(x₀), typeof(x₀), Tuple{eltype(x₀),eltype(x₀)}, Missing}
    @test equtype(prob) == ODE

    @test timespan(prob) == (t₀,t₁)
    @test initialtime(prob) == t₀
    @test finaltime(prob) == t₁
    @test timestep(prob) == Δt

    @test initial_conditions(prob) == (t = t₀, q = x₀)
    @test parameters(prob) == NullParameters()
    @test equation(prob) == ode
    @test nsamples(prob) == 1

    @test functions(prob) == functions(ode)

    @test prob == EquationProblem(ode, (t₀,t₁), Δt, ode_ics)
    @test prob == EquationProblem(ode, (t₀,t₁), Δt, ode_ics; parameters = NullParameters())
    @test prob == EquationProblem(ode, (t₀,t₁), Δt, ode_ics, NullParameters())
    @test prob == EquationProblem(ode, (t₀,t₁), Δt, ode_ics, nothing)

    @test prob == similar(prob)
    @test prob == similar(prob, (t₀,t₁))
    @test prob == similar(prob, (t₀,t₁), Δt)
    @test prob == similar(prob, (t₀,t₁), Δt, ode_ics)
    @test prob == similar(prob, (t₀,t₁), Δt, ode_ics, NullParameters())
    @test prob == similar(prob; timespan = (t₀,t₁))
    @test prob == similar(prob; timestep = Δt)
    @test prob == similar(prob; ics = ode_ics)
    @test prob == similar(prob; ics = ode_ics_raw)
    @test prob == similar(prob; ics = Tuple(ode_ics))
    @test prob == similar(prob; ics = Tuple(ode_ics_raw))
    @test prob == similar(prob; parameters = NullParameters())

end


@testset "$(rpad("ODE Problem",80))" begin

    ode  = ODE(ode_eqs...)
    prob = EquationProblem(ode, (t₀,t₁), Δt, ode_ics)

    @test typeof(prob) <: EquationProblem
    @test typeof(prob) <: ODEProblem
    @test typeof(prob) <: AbstractProblemODE
    @test equtype(prob) == ODE

    @test periodicity(prob).q == periodicity(equation(prob))

    @test prob == ODEProblem(ode_eqs..., (t₀,t₁), Δt, ode_ics)
    @test prob == ODEProblem(ode_eqs..., (t₀,t₁), Δt, ode_ics; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    @test prob == ODEProblem(ode_eqs..., (t₀,t₁), Δt, ode_ics...)
    @test prob == ODEProblem(ode_eqs..., (t₀,t₁), Δt, ode_ics...; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    @test prob == ODEProblem(ode_eqs..., (t₀,t₁), Δt, ode_ics_raw...)
    @test prob == ODEProblem(ode_eqs..., (t₀,t₁), Δt, ode_ics_raw...; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())

end


@testset "$(rpad("SODE Problem",80))" begin

    sode = SODE(sode_eqs)
    prob = EquationProblem(sode, (t₀,t₁), Δt, ode_ics)

    @test typeof(prob) <: EquationProblem
    @test typeof(prob) <: SODEProblem
    @test equtype(prob) == SODE

    @test periodicity(prob).q == periodicity(equation(prob))

    @test prob == SODEProblem(sode_eqs, (t₀,t₁), Δt, ode_ics...)
    @test prob == SODEProblem(sode_eqs, (t₀,t₁), Δt, ode_ics...; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    @test prob == SODEProblem(sode_eqs, (t₀,t₁), Δt, ode_ics_raw...)
    @test prob == SODEProblem(sode_eqs, (t₀,t₁), Δt, ode_ics_raw...; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())


    sode = SODE(sode_eqs, sode_sols)
    prob = EquationProblem(sode, (t₀,t₁), Δt, ode_ics)

    @test typeof(prob) <: EquationProblem
    @test typeof(prob) <: SODEProblem
    @test equtype(prob) == SODE

    @test prob == SODEProblem(sode_eqs, sode_sols, (t₀,t₁), Δt, ode_ics)
    @test prob == SODEProblem(sode_eqs, sode_sols, (t₀,t₁), Δt, ode_ics; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    @test prob == SODEProblem(sode_eqs, sode_sols, (t₀,t₁), Δt, ode_ics...)
    @test prob == SODEProblem(sode_eqs, sode_sols, (t₀,t₁), Δt, ode_ics...; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    @test prob == SODEProblem(sode_eqs, sode_sols, (t₀,t₁), Δt, ode_ics_raw...)
    @test prob == SODEProblem(sode_eqs, sode_sols, (t₀,t₁), Δt, ode_ics_raw...; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())


    prob1 = SubstepProblem(prob, one(timestep(prob))/2, 1)
    prob2 = SubstepProblem(prob, one(timestep(prob)), 2)
    prob3 = SubstepProblem(prob, one(timestep(prob))/2, 1)

    @test solutions(prob).q[1] == solutions(prob1).q
    @test solutions(prob).q[2] == solutions(prob2).q
    @test solutions(prob).q[1] == solutions(prob3).q

    @test problem(prob1) == prob
    @test problem(prob2) == prob
    @test problem(prob3) == prob

    @test timestep(prob1) == timestep(prob)/2
    @test timestep(prob2) == timestep(prob)
    @test timestep(prob3) == timestep(prob)/2

    @test datatype(prob1) == datatype(prob)
    @test timetype(prob1) == timetype(prob)
    @test arrtype(prob1) == arrtype(prob)
    @test equtype(prob1) == equtype(prob)

    @test ntime(prob1) == ntime(prob)
    @test timespan(prob1) == timespan(prob)
    @test timestep(prob1) == timestep(prob)/2

    @test equation(prob1) == equation(prob)
    @test invariants(prob1) == invariants(prob)
    @test parameters(prob1) == parameters(prob)
    @test nsamples(prob1) == nsamples(prob)

end


@testset "$(rpad("PODE Problem",80))" begin

    pode = PODE(pode_eqs...)
    prob = EquationProblem(pode, (t₀,t₁), Δt, pode_ics)

    @test typeof(prob) <: EquationProblem
    @test typeof(prob) <: PODEProblem
    @test typeof(prob) <: AbstractProblemPODE
    @test equtype(prob) == PODE

    @test periodicity(prob).q == periodicity(equation(prob))
    @test periodicity(prob).p == NullPeriodicity()

    @test prob == PODEProblem(pode_eqs..., (t₀,t₁), Δt, pode_ics)
    @test prob == PODEProblem(pode_eqs..., (t₀,t₁), Δt, pode_ics; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    @test prob == PODEProblem(pode_eqs..., (t₀,t₁), Δt, pode_ics...)
    @test prob == PODEProblem(pode_eqs..., (t₀,t₁), Δt, pode_ics...; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    @test prob == PODEProblem(pode_eqs..., (t₀,t₁), Δt, pode_ics_raw...)
    @test prob == PODEProblem(pode_eqs..., (t₀,t₁), Δt, pode_ics_raw...; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())

end


@testset "$(rpad("IODE Problem",80))" begin

    iode = IODE(iode_eqs...)
    prob = EquationProblem(iode, (t₀,t₁), Δt, iode_ics)

    @test typeof(prob) <: EquationProblem
    @test typeof(prob) <: IODEProblem
    @test typeof(prob) <: AbstractProblemIODE
    @test equtype(prob) == IODE

    @test periodicity(prob).q == periodicity(equation(prob))
    @test periodicity(prob).p == NullPeriodicity()
    @test periodicity(prob).v == NullPeriodicity()

    @test prob == IODEProblem(iode_eqs..., (t₀,t₁), Δt, iode_ics)
    @test prob == IODEProblem(iode_eqs..., (t₀,t₁), Δt, iode_ics; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    @test prob == IODEProblem(iode_eqs..., (t₀,t₁), Δt, iode_ics...)
    @test prob == IODEProblem(iode_eqs..., (t₀,t₁), Δt, iode_ics...; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    @test prob == IODEProblem(iode_eqs..., (t₀,t₁), Δt, iode_ics_raw...)
    @test prob == IODEProblem(iode_eqs..., (t₀,t₁), Δt, iode_ics_raw...; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())

end


@testset "$(rpad("HODE Problem",80))" begin

    hode = HODE(hode_eqs...)
    prob = EquationProblem(hode, (t₀,t₁), Δt, hode_ics)

    @test typeof(prob) <: EquationProblem
    @test typeof(prob) <: HODEProblem
    @test typeof(prob) <: AbstractProblemPODE
    @test typeof(prob) <: AbstractProblemHODE
    @test equtype(prob) == HODE

    @test periodicity(prob).q == periodicity(equation(prob))
    @test periodicity(prob).p == NullPeriodicity()

    @test prob == HODEProblem(hode_eqs..., (t₀,t₁), Δt, hode_ics)
    @test prob == HODEProblem(hode_eqs..., (t₀,t₁), Δt, hode_ics; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    @test prob == HODEProblem(hode_eqs..., (t₀,t₁), Δt, hode_ics...)
    @test prob == HODEProblem(hode_eqs..., (t₀,t₁), Δt, hode_ics...; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    @test prob == HODEProblem(hode_eqs..., (t₀,t₁), Δt, hode_ics_raw...)
    @test prob == HODEProblem(hode_eqs..., (t₀,t₁), Δt, hode_ics_raw...; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())

end


@testset "$(rpad("LODE Problem",80))" begin

    lode = LODE(lode_eqs...)
    prob = EquationProblem(lode, (t₀,t₁), Δt, lode_ics)

    @test typeof(prob) <: EquationProblem
    @test typeof(prob) <: LODEProblem
    @test typeof(prob) <: AbstractProblemIODE
    @test typeof(prob) <: AbstractProblemLODE
    @test equtype(prob) == LODE

    @test periodicity(prob).q == periodicity(equation(prob))
    @test periodicity(prob).p == NullPeriodicity()
    @test periodicity(prob).v == NullPeriodicity()

    @test prob == LODEProblem(lode_eqs..., (t₀,t₁), Δt, lode_ics)
    @test prob == LODEProblem(lode_eqs..., (t₀,t₁), Δt, lode_ics; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    @test prob == LODEProblem(lode_eqs..., (t₀,t₁), Δt, lode_ics...)
    @test prob == LODEProblem(lode_eqs..., (t₀,t₁), Δt, lode_ics...; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    @test prob == LODEProblem(lode_eqs..., (t₀,t₁), Δt, lode_ics_raw...)
    @test prob == LODEProblem(lode_eqs..., (t₀,t₁), Δt, lode_ics_raw...; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())

end


@testset "$(rpad("DAE Problem",80))" begin

    dae  = DAE(dae_eqs...)
    prob = EquationProblem(dae, (t₀,t₁), Δt, dae_ics_full)

    @test typeof(prob) <: EquationProblem
    @test typeof(prob) <: DAEProblem
    @test typeof(prob) <: AbstractProblemODE
    @test typeof(prob) <: AbstractProblemDAE
    @test equtype(prob) == DAE

    @test prob == DAEProblem(dae_eqs..., (t₀,t₁), Δt, dae_ics)
    @test prob == DAEProblem(dae_eqs..., (t₀,t₁), Δt, dae_ics; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    @test prob == DAEProblem(dae_eqs..., (t₀,t₁), Δt, dae_ics_full)
    @test prob == DAEProblem(dae_eqs..., (t₀,t₁), Δt, dae_ics_full; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    @test prob == DAEProblem(dae_eqs..., (t₀,t₁), Δt, dae_ics...)
    @test prob == DAEProblem(dae_eqs..., (t₀,t₁), Δt, dae_ics...; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    @test prob == DAEProblem(dae_eqs..., (t₀,t₁), Δt, dae_ics_full...)
    @test prob == DAEProblem(dae_eqs..., (t₀,t₁), Δt, dae_ics_full...; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    @test prob == DAEProblem(dae_eqs..., (t₀,t₁), Δt, dae_ics_raw...)
    @test prob == DAEProblem(dae_eqs..., (t₀,t₁), Δt, dae_ics_raw...; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())

    # @test initial_conditions(dae) == (t₀, [x₀], [λ₀], [λ₀])
    # @test eltype(dae)
    # @test arrtype(dae)
    # @test axes(dae) == axes(x₀)
    # @test ndims(dae) == 2
    # @test nsamples(dae) == 1
    # @test nconstraints(dae) == 1


    dae  = DAE(dae_eqs_full...)
    prob = EquationProblem(dae, (t₀,t₁), Δt, dae_ics_full)

    @test typeof(prob) <: EquationProblem
    @test typeof(prob) <: DAEProblem
    @test equtype(prob) == DAE

    @test prob == DAEProblem(dae_eqs_full..., (t₀,t₁), Δt, dae_ics)
    @test prob == DAEProblem(dae_eqs_full..., (t₀,t₁), Δt, dae_ics; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    @test prob == DAEProblem(dae_eqs_full..., (t₀,t₁), Δt, dae_ics_full)
    @test prob == DAEProblem(dae_eqs_full..., (t₀,t₁), Δt, dae_ics_full; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    @test prob == DAEProblem(dae_eqs_full..., (t₀,t₁), Δt, dae_ics...)
    @test prob == DAEProblem(dae_eqs_full..., (t₀,t₁), Δt, dae_ics...; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    @test prob == DAEProblem(dae_eqs_full..., (t₀,t₁), Δt, dae_ics_full...)
    @test prob == DAEProblem(dae_eqs_full..., (t₀,t₁), Δt, dae_ics_full...; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    @test prob == DAEProblem(dae_eqs_full..., (t₀,t₁), Δt, dae_ics_raw...)
    @test prob == DAEProblem(dae_eqs_full..., (t₀,t₁), Δt, dae_ics_raw...; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())

end


@testset "$(rpad("PDAE Problem",80))" begin

    pdae = PDAE(pdae_eqs...)
    prob = EquationProblem(pdae, (t₀,t₁), Δt, pdae_ics_full)

    @test typeof(prob) <: EquationProblem
    @test typeof(prob) <: PDAEProblem
    @test typeof(prob) <: AbstractProblemPODE
    @test typeof(prob) <: AbstractProblemPDAE
    @test equtype(prob) == PDAE

    @test prob == PDAEProblem(pdae_eqs..., (t₀,t₁), Δt, pdae_ics)
    @test prob == PDAEProblem(pdae_eqs..., (t₀,t₁), Δt, pdae_ics; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    @test prob == PDAEProblem(pdae_eqs..., (t₀,t₁), Δt, pdae_ics_full)
    @test prob == PDAEProblem(pdae_eqs..., (t₀,t₁), Δt, pdae_ics_full; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    @test prob == PDAEProblem(pdae_eqs..., (t₀,t₁), Δt, pdae_ics...)
    @test prob == PDAEProblem(pdae_eqs..., (t₀,t₁), Δt, pdae_ics...; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    @test prob == PDAEProblem(pdae_eqs..., (t₀,t₁), Δt, pdae_ics_full...)
    @test prob == PDAEProblem(pdae_eqs..., (t₀,t₁), Δt, pdae_ics_full...; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    @test prob == PDAEProblem(pdae_eqs..., (t₀,t₁), Δt, pdae_ics_raw...)
    @test prob == PDAEProblem(pdae_eqs..., (t₀,t₁), Δt, pdae_ics_raw...; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())

    # @test initial_conditions(dae) == (t₀, [x₀], [λ₀], [λ₀])
    # @test eltype(dae)
    # @test arrtype(dae)
    # @test axes(dae) == axes(x₀)
    # @test ndims(dae) == 2
    # @test nsamples(dae) == 1
    # @test nconstraints(dae) == 1


    pdae = PDAE(pdae_eqs_full...)
    prob = EquationProblem(pdae, (t₀,t₁), Δt, pdae_ics_full)

    @test typeof(prob) <: EquationProblem
    @test typeof(prob) <: PDAEProblem
    @test equtype(prob) == PDAE

    @test prob == PDAEProblem(pdae_eqs_full..., (t₀,t₁), Δt, pdae_ics)
    @test prob == PDAEProblem(pdae_eqs_full..., (t₀,t₁), Δt, pdae_ics; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    @test prob == PDAEProblem(pdae_eqs_full..., (t₀,t₁), Δt, pdae_ics_full)
    @test prob == PDAEProblem(pdae_eqs_full..., (t₀,t₁), Δt, pdae_ics_full; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    @test prob == PDAEProblem(pdae_eqs_full..., (t₀,t₁), Δt, pdae_ics...)
    @test prob == PDAEProblem(pdae_eqs_full..., (t₀,t₁), Δt, pdae_ics...; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    @test prob == PDAEProblem(pdae_eqs_full..., (t₀,t₁), Δt, pdae_ics_full...)
    @test prob == PDAEProblem(pdae_eqs_full..., (t₀,t₁), Δt, pdae_ics_full...; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    @test prob == PDAEProblem(pdae_eqs_full..., (t₀,t₁), Δt, pdae_ics_raw...)
    @test prob == PDAEProblem(pdae_eqs_full..., (t₀,t₁), Δt, pdae_ics_raw...; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())

end


@testset "$(rpad("HDAE Problem",80))" begin

    hdae = HDAE(hdae_eqs...)
    prob = EquationProblem(hdae, (t₀,t₁), Δt, hdae_ics_full)

    @test typeof(prob) <: EquationProblem
    @test typeof(prob) <: HDAEProblem
    @test typeof(prob) <: AbstractProblemPODE
    @test typeof(prob) <: AbstractProblemHODE
    @test typeof(prob) <: AbstractProblemPDAE
    @test equtype(prob) == HDAE

    @test prob == HDAEProblem(hdae_eqs..., (t₀,t₁), Δt, hdae_ics)
    @test prob == HDAEProblem(hdae_eqs..., (t₀,t₁), Δt, hdae_ics; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    @test prob == HDAEProblem(hdae_eqs..., (t₀,t₁), Δt, hdae_ics_full)
    @test prob == HDAEProblem(hdae_eqs..., (t₀,t₁), Δt, hdae_ics_full; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    @test prob == HDAEProblem(hdae_eqs..., (t₀,t₁), Δt, hdae_ics...)
    @test prob == HDAEProblem(hdae_eqs..., (t₀,t₁), Δt, hdae_ics...; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    @test prob == HDAEProblem(hdae_eqs..., (t₀,t₁), Δt, hdae_ics_full...)
    @test prob == HDAEProblem(hdae_eqs..., (t₀,t₁), Δt, hdae_ics_full...; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    @test prob == HDAEProblem(hdae_eqs..., (t₀,t₁), Δt, hdae_ics_raw...)
    @test prob == HDAEProblem(hdae_eqs..., (t₀,t₁), Δt, hdae_ics_raw...; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())

    # @test initial_conditions(dae) == (t₀, [x₀], [λ₀], [λ₀])
    # @test eltype(dae)
    # @test arrtype(dae)
    # @test axes(dae) == axes(x₀)
    # @test ndims(dae) == 2
    # @test nsamples(dae) == 1
    # @test nconstraints(dae) == 1


    hdae = HDAE(hdae_eqs_full...)
    prob = EquationProblem(hdae, (t₀,t₁), Δt, hdae_ics_full)

    @test typeof(prob) <: EquationProblem
    @test typeof(prob) <: HDAEProblem
    @test equtype(prob) == HDAE

    @test prob == HDAEProblem(hdae_eqs_full..., (t₀,t₁), Δt, hdae_ics)
    @test prob == HDAEProblem(hdae_eqs_full..., (t₀,t₁), Δt, hdae_ics; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    @test prob == HDAEProblem(hdae_eqs_full..., (t₀,t₁), Δt, hdae_ics_full)
    @test prob == HDAEProblem(hdae_eqs_full..., (t₀,t₁), Δt, hdae_ics_full; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    @test prob == HDAEProblem(hdae_eqs_full..., (t₀,t₁), Δt, hdae_ics...)
    @test prob == HDAEProblem(hdae_eqs_full..., (t₀,t₁), Δt, hdae_ics...; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    @test prob == HDAEProblem(hdae_eqs_full..., (t₀,t₁), Δt, hdae_ics_full...)
    @test prob == HDAEProblem(hdae_eqs_full..., (t₀,t₁), Δt, hdae_ics_full...; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    @test prob == HDAEProblem(hdae_eqs_full..., (t₀,t₁), Δt, hdae_ics_raw...)
    @test prob == HDAEProblem(hdae_eqs_full..., (t₀,t₁), Δt, hdae_ics_raw...; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())

end


@testset "$(rpad("IDAE Problem",80))" begin

    idae = IDAE(idae_eqs...; v̄=idae_v, f̄=idae_f)
    prob = EquationProblem(idae, (t₀,t₁), Δt, idae_ics_full)

    @test typeof(prob) <: EquationProblem
    @test typeof(prob) <: IDAEProblem
    @test typeof(prob) <: AbstractProblemIODE
    @test typeof(prob) <: AbstractProblemIDAE
    @test equtype(prob) == IDAE

    @test prob == IDAEProblem(idae_eqs..., (t₀,t₁), Δt, idae_ics)
    @test prob == IDAEProblem(idae_eqs..., (t₀,t₁), Δt, idae_ics; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    @test prob == IDAEProblem(idae_eqs..., (t₀,t₁), Δt, idae_ics_full)
    @test prob == IDAEProblem(idae_eqs..., (t₀,t₁), Δt, idae_ics_full; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    @test prob == IDAEProblem(idae_eqs..., (t₀,t₁), Δt, idae_ics...)
    @test prob == IDAEProblem(idae_eqs..., (t₀,t₁), Δt, idae_ics...; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    @test prob == IDAEProblem(idae_eqs..., (t₀,t₁), Δt, idae_ics_full...)
    @test prob == IDAEProblem(idae_eqs..., (t₀,t₁), Δt, idae_ics_full...; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    @test prob == IDAEProblem(idae_eqs..., (t₀,t₁), Δt, idae_ics_raw...)
    @test prob == IDAEProblem(idae_eqs..., (t₀,t₁), Δt, idae_ics_raw...; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())

    # @test initial_conditions(dae) == (t₀, [x₀], [λ₀], [λ₀])
    # @test eltype(dae)
    # @test arrtype(dae)
    # @test axes(dae) == axes(x₀)
    # @test ndims(dae) == 2
    # @test nsamples(dae) == 1
    # @test nconstraints(dae) == 1


    idae = IDAE(idae_eqs_full...)
    prob = EquationProblem(idae, (t₀,t₁), Δt, idae_ics_full)

    @test typeof(prob) <: EquationProblem
    @test typeof(prob) <: IDAEProblem
    @test equtype(prob) == IDAE

    @test prob == IDAEProblem(idae_eqs_full..., (t₀,t₁), Δt, idae_ics)
    @test prob == IDAEProblem(idae_eqs_full..., (t₀,t₁), Δt, idae_ics; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    @test prob == IDAEProblem(idae_eqs_full..., (t₀,t₁), Δt, idae_ics_full)
    @test prob == IDAEProblem(idae_eqs_full..., (t₀,t₁), Δt, idae_ics_full; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    @test prob == IDAEProblem(idae_eqs_full..., (t₀,t₁), Δt, idae_ics)
    @test prob == IDAEProblem(idae_eqs_full..., (t₀,t₁), Δt, idae_ics; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    @test prob == IDAEProblem(idae_eqs_full..., (t₀,t₁), Δt, idae_ics_full)
    @test prob == IDAEProblem(idae_eqs_full..., (t₀,t₁), Δt, idae_ics_full; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    @test prob == IDAEProblem(idae_eqs_full..., (t₀,t₁), Δt, idae_ics_raw...)
    @test prob == IDAEProblem(idae_eqs_full..., (t₀,t₁), Δt, idae_ics_raw...; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())

end


@testset "$(rpad("LDAE Problem",80))" begin

    ldae = LDAE(ldae_eqs...)
    prob = EquationProblem(ldae, (t₀,t₁), Δt, ldae_ics_full)

    @test typeof(prob) <: EquationProblem
    @test typeof(prob) <: LDAEProblem
    @test typeof(prob) <: AbstractProblemIODE
    @test typeof(prob) <: AbstractProblemLODE
    @test typeof(prob) <: AbstractProblemIDAE
    @test equtype(prob) == LDAE

    @test prob == LDAEProblem(ldae_eqs..., (t₀,t₁), Δt, ldae_ics)
    @test prob == LDAEProblem(ldae_eqs..., (t₀,t₁), Δt, ldae_ics; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    @test prob == LDAEProblem(ldae_eqs..., (t₀,t₁), Δt, ldae_ics_full)
    @test prob == LDAEProblem(ldae_eqs..., (t₀,t₁), Δt, ldae_ics_full; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    @test prob == LDAEProblem(ldae_eqs..., (t₀,t₁), Δt, ldae_ics...)
    @test prob == LDAEProblem(ldae_eqs..., (t₀,t₁), Δt, ldae_ics...; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    @test prob == LDAEProblem(ldae_eqs..., (t₀,t₁), Δt, ldae_ics_full...)
    @test prob == LDAEProblem(ldae_eqs..., (t₀,t₁), Δt, ldae_ics_full...; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    @test prob == LDAEProblem(ldae_eqs..., (t₀,t₁), Δt, ldae_ics_raw...)
    @test prob == LDAEProblem(ldae_eqs..., (t₀,t₁), Δt, ldae_ics_raw...; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())

    # @test initial_conditions(dae) == (t₀, [x₀], [λ₀], [λ₀])
    # @test eltype(dae)
    # @test arrtype(dae)
    # @test axes(dae) == axes(x₀)
    # @test ndims(dae) == 2
    # @test nsamples(dae) == 1
    # @test nconstraints(dae) == 1


    ldae = LDAE(ldae_eqs_full...)
    prob = EquationProblem(ldae, (t₀,t₁), Δt, ldae_ics_full)

    @test typeof(prob) <: EquationProblem
    @test typeof(prob) <: LDAEProblem
    @test equtype(prob) == LDAE

    @test prob == LDAEProblem(ldae_eqs_full..., (t₀,t₁), Δt, ldae_ics)
    @test prob == LDAEProblem(ldae_eqs_full..., (t₀,t₁), Δt, ldae_ics; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    @test prob == LDAEProblem(ldae_eqs_full..., (t₀,t₁), Δt, ldae_ics_full)
    @test prob == LDAEProblem(ldae_eqs_full..., (t₀,t₁), Δt, ldae_ics_full; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    @test prob == LDAEProblem(ldae_eqs_full..., (t₀,t₁), Δt, ldae_ics...)
    @test prob == LDAEProblem(ldae_eqs_full..., (t₀,t₁), Δt, ldae_ics...; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    @test prob == LDAEProblem(ldae_eqs_full..., (t₀,t₁), Δt, ldae_ics_full...)
    @test prob == LDAEProblem(ldae_eqs_full..., (t₀,t₁), Δt, ldae_ics_full...; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    @test prob == LDAEProblem(ldae_eqs_full..., (t₀,t₁), Δt, ldae_ics_raw...)
    @test prob == LDAEProblem(ldae_eqs_full..., (t₀,t₁), Δt, ldae_ics_raw...; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())

end


@testset "$(rpad("SDE Problem",80))" begin

    sde  = SDE(sde_v, sde_B, TestNoise())
    prob = EquationProblem(sde, (t₀,t₁), Δt, sde_ics)

    @test typeof(prob) <: EquationProblem
    @test typeof(prob) <: SDEProblem
    @test typeof(prob) <: AbstractProblemSDE
    @test equtype(prob) == SDE

    @test periodicity(prob).q == periodicity(equation(prob))

    @test prob == SDEProblem(sde_v, sde_B, TestNoise(), (t₀,t₁), Δt, sde_ics)
    @test prob == SDEProblem(sde_v, sde_B, TestNoise(), (t₀,t₁), Δt, sde_ics; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    @test prob == SDEProblem(sde_v, sde_B, TestNoise(), (t₀,t₁), Δt, sde_ics...)
    @test prob == SDEProblem(sde_v, sde_B, TestNoise(), (t₀,t₁), Δt, sde_ics...; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    @test prob == SDEProblem(sde_v, sde_B, TestNoise(), (t₀,t₁), Δt, sde_ics_raw...)
    @test prob == SDEProblem(sde_v, sde_B, TestNoise(), (t₀,t₁), Δt, sde_ics_raw...; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())

end


@testset "$(rpad("PSDE Problem",80))" begin

    psde = PSDE(psde_v, psde_f, psde_B, psde_G, TestNoise())
    prob = EquationProblem(psde, (t₀,t₁), Δt, psde_ics)

    @test typeof(prob) <: EquationProblem
    @test typeof(prob) <: PSDEProblem
    @test typeof(prob) <: AbstractProblemPSDE
    @test equtype(prob) == PSDE

    @test periodicity(prob).q == periodicity(equation(prob))
    @test periodicity(prob).p == NullPeriodicity()

    @test prob == PSDEProblem(psde_v, psde_f, psde_B, psde_G, TestNoise(), (t₀,t₁), Δt, psde_ics)
    @test prob == PSDEProblem(psde_v, psde_f, psde_B, psde_G, TestNoise(), (t₀,t₁), Δt, psde_ics; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    @test prob == PSDEProblem(psde_v, psde_f, psde_B, psde_G, TestNoise(), (t₀,t₁), Δt, psde_ics...)
    @test prob == PSDEProblem(psde_v, psde_f, psde_B, psde_G, TestNoise(), (t₀,t₁), Δt, psde_ics...; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    @test prob == PSDEProblem(psde_v, psde_f, psde_B, psde_G, TestNoise(), (t₀,t₁), Δt, psde_ics_raw...)
    @test prob == PSDEProblem(psde_v, psde_f, psde_B, psde_G, TestNoise(), (t₀,t₁), Δt, psde_ics_raw...; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())

end


@testset "$(rpad("SPSDE Problem",80))" begin

    psde = SPSDE(spsde_v, spsde_f1, spsde_f2, spsde_B, spsde_G1, spsde_G2, TestNoise())
    prob = EquationProblem(psde, (t₀,t₁), Δt, spsde_ics)

    @test typeof(prob) <: EquationProblem
    @test typeof(prob) <: SPSDEProblem
    @test typeof(prob) <: AbstractProblemSPSDE
    @test equtype(prob) == SPSDE

    @test periodicity(prob).q == periodicity(equation(prob))
    @test periodicity(prob).p == NullPeriodicity()

    @test prob == SPSDEProblem(spsde_v, spsde_f1, spsde_f2, spsde_B, spsde_G1, spsde_G2, TestNoise(), (t₀,t₁), Δt, spsde_ics)
    @test prob == SPSDEProblem(spsde_v, spsde_f1, spsde_f2, spsde_B, spsde_G1, spsde_G2, TestNoise(), (t₀,t₁), Δt, spsde_ics; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    @test prob == SPSDEProblem(spsde_v, spsde_f1, spsde_f2, spsde_B, spsde_G1, spsde_G2, TestNoise(), (t₀,t₁), Δt, spsde_ics...)
    @test prob == SPSDEProblem(spsde_v, spsde_f1, spsde_f2, spsde_B, spsde_G1, spsde_G2, TestNoise(), (t₀,t₁), Δt, spsde_ics...; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    @test prob == SPSDEProblem(spsde_v, spsde_f1, spsde_f2, spsde_B, spsde_G1, spsde_G2, TestNoise(), (t₀,t₁), Δt, spsde_ics_raw...)
    @test prob == SPSDEProblem(spsde_v, spsde_f1, spsde_f2, spsde_B, spsde_G1, spsde_G2, TestNoise(), (t₀,t₁), Δt, spsde_ics_raw...; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())

end
