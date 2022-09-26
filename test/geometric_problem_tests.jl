using GeometricBase
using GeometricEquations
using Test

include("functions.jl")
include("initial_conditions.jl")


@testset "$(rpad("Geometric Problem",80))" begin

    ics  = (q=x₀,)
    ode  = ODE(ode_v)
    prob = GeometricProblem(ode, (t₀,t₁), Δt, ics)

    @test typeof(prob) <: GeometricProblem
    @test typeof(prob) <: ODEProblem
    @test typeof(prob).parameters[1] == ODE
    @test typeof(prob).parameters[5] == typeof(ode)

    @test datatype(prob) == eltype(x₀)
    @test timetype(prob) == typeof(t₀)
    @test arrtype(prob) == typeof(x₀)
    @test equtype(prob) == ODE

    @test tspan(prob) == (t₀,t₁)
    @test tbegin(prob) == t₀
    @test tend(prob) == t₁
    @test tstep(prob) == Δt

    @test initial_conditions(prob) == (t = t₀, q = x₀)
    @test parameters(prob) == NullParameters()
    @test equation(prob) == ode
    @test nsamples(prob) == 1

    @test functions(prob) == functions(ode, parameters(prob))

    prob1 = GeometricProblem(ode, (t₀,t₁), Δt, ics)
    prob2 = GeometricProblem(ode, (t₀,t₁), Δt, ics; parameters=NullParameters())
    prob3 = GeometricProblem(ode, (t₀,t₁), Δt, ics, NullParameters())
    prob4 = GeometricProblem(ode, (t₀,t₁), Δt, ics, nothing)

    @test prob == prob1
    @test prob == prob2
    @test prob == prob3
    @test prob == prob4

    @test prob == similar(prob)
    @test prob == similar(prob, (t₀,t₁))
    @test prob == similar(prob, (t₀,t₁), Δt)
    @test prob == similar(prob, (t₀,t₁), Δt, ics)
    @test prob == similar(prob, (t₀,t₁), Δt, ics, NullParameters())
    @test prob == similar(prob; tspan=(t₀,t₁))
    @test prob == similar(prob; tstep=Δt)
    @test prob == similar(prob; ics=ics)
    @test prob == similar(prob; parameters=NullParameters())

end


@testset "$(rpad("ODE Problem",80))" begin

    ics  = (q=x₀,)
    ode  = ODE(ode_v)
    prob = GeometricProblem(ode, (t₀,t₁), Δt, ics)

    @test typeof(prob) <: GeometricProblem
    @test typeof(prob) <: ODEProblem
    @test equtype(prob) == ODE

    @test periodicity(prob).q == periodicity(equation(prob))

    prob1 = ODEProblem(ode_v, (t₀,t₁), Δt, ics)
    prob2 = ODEProblem(ode_v, (t₀,t₁), Δt, ics; invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity())
    prob3 = ODEProblem(ode_v, (t₀,t₁), Δt, ics...)
    prob4 = ODEProblem(ode_v, (t₀,t₁), Δt, ics...; invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity())

    @test prob1 == prob
    @test prob2 == prob
    @test prob3 == prob
    @test prob4 == prob
    
end


@testset "$(rpad("SODE Problem",80))" begin

    ics  = (q=x₀,)
    sode = SODE(sode_v)
    prob = GeometricProblem(sode, (t₀,t₁), Δt, ics)

    @test typeof(prob) <: GeometricProblem
    @test typeof(prob) <: SODEProblem
    @test equtype(prob) == SODE

    @test periodicity(prob).q == periodicity(equation(prob))

    prob1 = SODEProblem(sode_v, (t₀,t₁), Δt, ics)
    prob2 = SODEProblem(sode_v, (t₀,t₁), Δt, ics; invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity())
    prob3 = SODEProblem(sode_v, (t₀,t₁), Δt, ics...)
    prob4 = SODEProblem(sode_v, (t₀,t₁), Δt, ics...; invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity())

    @test prob1 == prob
    @test prob2 == prob
    @test prob3 == prob
    @test prob4 == prob
    

    sode = SODE(sode_v, sode_q)
    prob = GeometricProblem(sode, (t₀,t₁), Δt, ics)

    @test typeof(prob) <: GeometricProblem
    @test typeof(prob) <: SODEProblem
    @test equtype(prob) == SODE

    prob1 = SODEProblem(sode_v, sode_q, (t₀,t₁), Δt, ics)
    prob2 = SODEProblem(sode_v, sode_q, (t₀,t₁), Δt, ics; invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity())
    prob3 = SODEProblem(sode_v, sode_q, (t₀,t₁), Δt, ics...)
    prob4 = SODEProblem(sode_v, sode_q, (t₀,t₁), Δt, ics...; invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity())

    @test prob1 == prob
    @test prob2 == prob
    @test prob3 == prob
    @test prob4 == prob
    
end


@testset "$(rpad("PODE Problem",80))" begin

    ics  = (q=q₀, p=p₀)
    pode = PODE(pode_eqs...)
    prob = GeometricProblem(pode, (t₀,t₁), Δt, ics)

    @test typeof(prob) <: GeometricProblem
    @test typeof(prob) <: PODEProblem
    @test equtype(prob) == PODE

    @test periodicity(prob).q == periodicity(equation(prob))
    @test periodicity(prob).p == NullPeriodicity()

    prob1 = PODEProblem(pode_eqs..., (t₀,t₁), Δt, ics)
    prob2 = PODEProblem(pode_eqs..., (t₀,t₁), Δt, ics; invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity())
    prob3 = PODEProblem(pode_eqs..., (t₀,t₁), Δt, ics...)
    prob4 = PODEProblem(pode_eqs..., (t₀,t₁), Δt, ics...; invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity())

    @test prob1 == prob
    @test prob2 == prob
    @test prob3 == prob
    @test prob4 == prob
    
end


@testset "$(rpad("IODE Problem",80))" begin

    ics  = (q=q₀, p=p₀, λ=λ₀)
    iode = IODE(iode_eqs...)
    prob = GeometricProblem(iode, (t₀,t₁), Δt, ics)

    @test typeof(prob) <: GeometricProblem
    @test typeof(prob) <: IODEProblem
    @test equtype(prob) == IODE

    @test periodicity(prob).q == periodicity(equation(prob))
    @test periodicity(prob).p == NullPeriodicity()
    @test periodicity(prob).λ == NullPeriodicity()

    prob1 = IODEProblem(iode_eqs..., (t₀,t₁), Δt, ics)
    prob2 = IODEProblem(iode_eqs..., (t₀,t₁), Δt, ics; invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity())
    prob3 = IODEProblem(iode_eqs..., (t₀,t₁), Δt, ics...)
    prob4 = IODEProblem(iode_eqs..., (t₀,t₁), Δt, ics...; invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity())

    @test prob1 == prob
    @test prob2 == prob
    @test prob3 == prob
    @test prob4 == prob
    
end


@testset "$(rpad("HODE Problem",80))" begin

    ics  = (q=q₀, p=p₀)
    hode = HODE(hode_eqs...)
    prob = GeometricProblem(hode, (t₀,t₁), Δt, ics)

    @test typeof(prob) <: GeometricProblem
    @test typeof(prob) <: HODEProblem
    @test equtype(prob) == HODE

    @test periodicity(prob).q == periodicity(equation(prob))
    @test periodicity(prob).p == NullPeriodicity()

    prob1 = HODEProblem(hode_eqs..., (t₀,t₁), Δt, ics)
    prob2 = HODEProblem(hode_eqs..., (t₀,t₁), Δt, ics; invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity())
    prob3 = HODEProblem(hode_eqs..., (t₀,t₁), Δt, ics...)
    prob4 = HODEProblem(hode_eqs..., (t₀,t₁), Δt, ics...; invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity())

    @test prob1 == prob
    @test prob2 == prob
    @test prob3 == prob
    @test prob4 == prob
    
end


@testset "$(rpad("LODE Problem",80))" begin

    ics  = (q=q₀, p=p₀, λ=λ₀)
    lode = LODE(lode_eqs...)
    prob = GeometricProblem(lode, (t₀,t₁), Δt, ics)

    @test typeof(prob) <: GeometricProblem
    @test typeof(prob) <: LODEProblem
    @test equtype(prob) == LODE

    @test periodicity(prob).q == periodicity(equation(prob))
    @test periodicity(prob).p == NullPeriodicity()
    @test periodicity(prob).λ == NullPeriodicity()

    prob1 = LODEProblem(lode_eqs..., (t₀,t₁), Δt, ics)
    prob2 = LODEProblem(lode_eqs..., (t₀,t₁), Δt, ics; invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity())
    prob3 = LODEProblem(lode_eqs..., (t₀,t₁), Δt, ics...)
    prob4 = LODEProblem(lode_eqs..., (t₀,t₁), Δt, ics...; invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity())

    @test prob1 == prob
    @test prob2 == prob
    @test prob3 == prob
    @test prob4 == prob
    
end


@testset "$(rpad("DAE Problem",80))" begin

    eqs  = (dae_v, dae_u, dae_ϕ)
    ics  = (q=x₀, λ=λ₀)
    dae  = DAE(eqs...)
    prob = GeometricProblem(dae, (t₀,t₁), Δt, ics)

    @test typeof(prob) <: GeometricProblem
    @test typeof(prob) <: DAEProblem
    @test equtype(prob) == DAE

    prob1 = DAEProblem(eqs..., (t₀,t₁), Δt, ics)
    prob2 = DAEProblem(eqs..., (t₀,t₁), Δt, ics; invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity())
    prob3 = DAEProblem(eqs..., (t₀,t₁), Δt, ics...)
    prob4 = DAEProblem(eqs..., (t₀,t₁), Δt, ics...; invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity())

    @test prob1 == prob
    @test prob2 == prob
    @test prob3 == prob
    @test prob4 == prob

    # @test initial_conditions(dae) == (t₀, [x₀], [λ₀], [λ₀])
    # @test eltype(dae)
    # @test arrtype(dae)
    # @test axes(dae) == axes(x₀)
    # @test ndims(dae) == 2
    # @test nsamples(dae) == 1
    # @test nconstraints(dae) == 1


    eqs = (dae_v, dae_u, dae_ϕ, dae_u, dae_ϕ)
    ics  = (q=x₀, λ=λ₀, μ=λ₀)
    dae  = DAE(eqs...)
    prob = GeometricProblem(dae, (t₀,t₁), Δt, ics)

    @test typeof(prob) <: GeometricProblem
    @test typeof(prob) <: DAEProblem
    @test equtype(prob) == DAE

    prob1 = DAEProblem(eqs..., (t₀,t₁), Δt, ics)
    prob2 = DAEProblem(eqs..., (t₀,t₁), Δt, ics; invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity())
    prob3 = DAEProblem(eqs..., (t₀,t₁), Δt, ics...)
    prob4 = DAEProblem(eqs..., (t₀,t₁), Δt, ics...; invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity())

    @test prob1 == prob
    @test prob2 == prob
    @test prob3 == prob
    @test prob4 == prob
    
end


@testset "$(rpad("PDAE Problem",80))" begin

    eqs  = (pdae_v, pdae_f, pdae_u, pdae_g, pdae_ϕ)
    ics  = (q=q₀, p=p₀, λ=λ₀)
    pdae = PDAE(eqs...)
    prob = GeometricProblem(pdae, (t₀,t₁), Δt, ics)

    @test typeof(prob) <: GeometricProblem
    @test typeof(prob) <: PDAEProblem
    @test equtype(prob) == PDAE

    prob1 = PDAEProblem(eqs..., (t₀,t₁), Δt, ics)
    prob2 = PDAEProblem(eqs..., (t₀,t₁), Δt, ics; invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity())
    prob3 = PDAEProblem(eqs..., (t₀,t₁), Δt, ics...)
    prob4 = PDAEProblem(eqs..., (t₀,t₁), Δt, ics...; invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity())

    @test prob1 == prob
    @test prob2 == prob
    @test prob3 == prob
    @test prob4 == prob

    # @test initial_conditions(dae) == (t₀, [x₀], [λ₀], [λ₀])
    # @test eltype(dae)
    # @test arrtype(dae)
    # @test axes(dae) == axes(x₀)
    # @test ndims(dae) == 2
    # @test nsamples(dae) == 1
    # @test nconstraints(dae) == 1


    eqs  = (pdae_v, pdae_f, pdae_u, pdae_g, pdae_ϕ, pdae_u, pdae_g, pdae_ϕ)
    ics  = (q=q₀, p=p₀, λ=λ₀, μ=λ₀)
    pdae = PDAE(eqs...)
    prob = GeometricProblem(pdae, (t₀,t₁), Δt, ics)

    @test typeof(prob) <: GeometricProblem
    @test typeof(prob) <: PDAEProblem
    @test equtype(prob) == PDAE

    prob1 = PDAEProblem(eqs..., (t₀,t₁), Δt, ics)
    prob2 = PDAEProblem(eqs..., (t₀,t₁), Δt, ics; invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity())
    prob3 = PDAEProblem(eqs..., (t₀,t₁), Δt, ics...)
    prob4 = PDAEProblem(eqs..., (t₀,t₁), Δt, ics...; invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity())

    @test prob1 == prob
    @test prob2 == prob
    @test prob3 == prob
    @test prob4 == prob
    
end


@testset "$(rpad("HDAE Problem",80))" begin

    eqs  = (pdae_v, pdae_f, pdae_u, pdae_g, pdae_ϕ, pdae_h)
    ics  = (q=q₀, p=p₀, λ=λ₀)
    hdae = HDAE(eqs...)
    prob = GeometricProblem(hdae, (t₀,t₁), Δt, ics)

    @test typeof(prob) <: GeometricProblem
    @test typeof(prob) <: HDAEProblem
    @test equtype(prob) == HDAE

    prob1 = HDAEProblem(eqs..., (t₀,t₁), Δt, ics)
    prob2 = HDAEProblem(eqs..., (t₀,t₁), Δt, ics; invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity())
    prob3 = HDAEProblem(eqs..., (t₀,t₁), Δt, ics...)
    prob4 = HDAEProblem(eqs..., (t₀,t₁), Δt, ics...; invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity())

    @test prob1 == prob
    @test prob2 == prob
    @test prob3 == prob
    @test prob4 == prob

    # @test initial_conditions(dae) == (t₀, [x₀], [λ₀], [λ₀])
    # @test eltype(dae)
    # @test arrtype(dae)
    # @test axes(dae) == axes(x₀)
    # @test ndims(dae) == 2
    # @test nsamples(dae) == 1
    # @test nconstraints(dae) == 1


    eqs  = (pdae_v, pdae_f, pdae_u, pdae_g, pdae_ϕ, pdae_u, pdae_g, pdae_ψ, pdae_h)
    ics  = (q=q₀, p=p₀, λ=λ₀, μ=λ₀)
    hdae = HDAE(eqs...)
    prob = GeometricProblem(hdae, (t₀,t₁), Δt, ics)

    @test typeof(prob) <: GeometricProblem
    @test typeof(prob) <: HDAEProblem
    @test equtype(prob) == HDAE

    prob1 = HDAEProblem(eqs..., (t₀,t₁), Δt, ics)
    prob2 = HDAEProblem(eqs..., (t₀,t₁), Δt, ics; invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity())
    prob3 = HDAEProblem(eqs..., (t₀,t₁), Δt, ics...)
    prob4 = HDAEProblem(eqs..., (t₀,t₁), Δt, ics...; invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity())

    @test prob1 == prob
    @test prob2 == prob
    @test prob3 == prob
    @test prob4 == prob
    
end


@testset "$(rpad("IDAE Problem",80))" begin

    eqs  = (idae_ϑ, idae_f, pdae_u, pdae_g, pdae_ϕ)
    ics  = (q=q₀, p=p₀, λ=λ₀)
    idae = IDAE(eqs...; v̄=idae_v, f̄=idae_f)
    prob = GeometricProblem(idae, (t₀,t₁), Δt, ics)

    @test typeof(prob) <: GeometricProblem
    @test typeof(prob) <: IDAEProblem
    @test equtype(prob) == IDAE

    prob1 = IDAEProblem(eqs..., (t₀,t₁), Δt, ics)
    prob2 = IDAEProblem(eqs..., (t₀,t₁), Δt, ics; invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity())
    prob3 = IDAEProblem(eqs..., (t₀,t₁), Δt, ics...)
    prob4 = IDAEProblem(eqs..., (t₀,t₁), Δt, ics...; invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity())

    @test prob1 == prob
    @test prob2 == prob
    @test prob3 == prob
    @test prob4 == prob

    # @test initial_conditions(dae) == (t₀, [x₀], [λ₀], [λ₀])
    # @test eltype(dae)
    # @test arrtype(dae)
    # @test axes(dae) == axes(x₀)
    # @test ndims(dae) == 2
    # @test nsamples(dae) == 1
    # @test nconstraints(dae) == 1


    eqs  = (idae_ϑ, idae_f, pdae_u, pdae_g, pdae_ϕ, pdae_u, pdae_g, pdae_ϕ)
    ics  = (q=q₀, p=p₀, λ=λ₀, μ=λ₀)
    idae = IDAE(eqs...)
    prob = GeometricProblem(idae, (t₀,t₁), Δt, ics)

    @test typeof(prob) <: GeometricProblem
    @test typeof(prob) <: IDAEProblem
    @test equtype(prob) == IDAE

    prob1 = IDAEProblem(eqs..., (t₀,t₁), Δt, ics)
    prob2 = IDAEProblem(eqs..., (t₀,t₁), Δt, ics; invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity())
    prob3 = IDAEProblem(eqs..., (t₀,t₁), Δt, ics...)
    prob4 = IDAEProblem(eqs..., (t₀,t₁), Δt, ics...; invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity())

    @test prob1 == prob
    @test prob2 == prob
    @test prob3 == prob
    @test prob4 == prob
    
end


@testset "$(rpad("LDAE Problem",80))" begin

    eqs  = (idae_ϑ, idae_f, pdae_u, pdae_g, pdae_ϕ, ldae_ω, ldae_l)
    ics  = (q=q₀, p=p₀, λ=λ₀)
    ldae = LDAE(eqs...)
    prob = GeometricProblem(ldae, (t₀,t₁), Δt, ics)

    @test typeof(prob) <: GeometricProblem
    @test typeof(prob) <: LDAEProblem
    @test equtype(prob) == LDAE

    prob1 = LDAEProblem(eqs..., (t₀,t₁), Δt, ics)
    prob2 = LDAEProblem(eqs..., (t₀,t₁), Δt, ics; invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity())
    prob3 = LDAEProblem(eqs..., (t₀,t₁), Δt, ics...)
    prob4 = LDAEProblem(eqs..., (t₀,t₁), Δt, ics...; invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity())

    @test prob1 == prob
    @test prob2 == prob
    @test prob3 == prob
    @test prob4 == prob

    # @test initial_conditions(dae) == (t₀, [x₀], [λ₀], [λ₀])
    # @test eltype(dae)
    # @test arrtype(dae)
    # @test axes(dae) == axes(x₀)
    # @test ndims(dae) == 2
    # @test nsamples(dae) == 1
    # @test nconstraints(dae) == 1


    eqs  = (idae_ϑ, idae_f, pdae_u, pdae_g, pdae_ϕ, pdae_u, pdae_g, pdae_ϕ, ldae_ω, ldae_l)
    ics  = (q=q₀, p=p₀, λ=λ₀, μ=λ₀)
    ldae = LDAE(eqs...)
    prob = GeometricProblem(ldae, (t₀,t₁), Δt, ics)

    @test typeof(prob) <: GeometricProblem
    @test typeof(prob) <: LDAEProblem
    @test equtype(prob) == LDAE

    prob1 = LDAEProblem(eqs..., (t₀,t₁), Δt, ics)
    prob2 = LDAEProblem(eqs..., (t₀,t₁), Δt, ics; invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity())
    prob3 = LDAEProblem(eqs..., (t₀,t₁), Δt, ics...)
    prob4 = LDAEProblem(eqs..., (t₀,t₁), Δt, ics...; invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity())

    @test prob1 == prob
    @test prob2 == prob
    @test prob3 == prob
    @test prob4 == prob
    
end


@testset "$(rpad("SDE Problem",80))" begin

    ics  = (q=x₀,)
    sde  = SDE(sde_v, sde_B)
    prob = GeometricProblem(sde, (t₀,t₁), Δt, ics)

    @test typeof(prob) <: GeometricProblem
    @test typeof(prob) <: SDEProblem
    @test equtype(prob) == SDE

    @test periodicity(prob).q == periodicity(equation(prob))

    prob1 = SDEProblem(sde_v, sde_B, (t₀,t₁), Δt, ics)
    prob2 = SDEProblem(sde_v, sde_B, (t₀,t₁), Δt, ics; invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity())
    prob3 = SDEProblem(sde_v, sde_B, (t₀,t₁), Δt, ics...)
    prob4 = SDEProblem(sde_v, sde_B, (t₀,t₁), Δt, ics...; invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity())

    @test prob1 == prob
    @test prob2 == prob
    @test prob3 == prob
    @test prob4 == prob
    
end


@testset "$(rpad("PSDE Problem",80))" begin

    ics  = (q=q₀,p=p₀)
    psde = PSDE(psde_v, psde_f, psde_B, psde_G)
    prob = GeometricProblem(psde, (t₀,t₁), Δt, ics)

    @test typeof(prob) <: GeometricProblem
    @test typeof(prob) <: PSDEProblem
    @test equtype(prob) == PSDE

    @test periodicity(prob).q == periodicity(equation(prob))
    @test periodicity(prob).p == NullPeriodicity()

    prob1 = PSDEProblem(psde_v, psde_f, psde_B, psde_G, (t₀,t₁), Δt, ics)
    prob2 = PSDEProblem(psde_v, psde_f, psde_B, psde_G, (t₀,t₁), Δt, ics; invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity())
    prob3 = PSDEProblem(psde_v, psde_f, psde_B, psde_G, (t₀,t₁), Δt, ics...)
    prob4 = PSDEProblem(psde_v, psde_f, psde_B, psde_G, (t₀,t₁), Δt, ics...; invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity())

    @test prob1 == prob
    @test prob2 == prob
    @test prob3 == prob
    @test prob4 == prob
    
end


@testset "$(rpad("SPSDE Problem",80))" begin

    ics  = (q=q₀,p=p₀)
    psde = SPSDE(spsde_v, spsde_f1, spsde_f2, spsde_B, spsde_G1, spsde_G2)
    prob = GeometricProblem(psde, (t₀,t₁), Δt, ics)

    @test typeof(prob) <: GeometricProblem
    @test typeof(prob) <: SPSDEProblem
    @test equtype(prob) == SPSDE

    @test periodicity(prob).q == periodicity(equation(prob))
    @test periodicity(prob).p == NullPeriodicity()

    prob1 = SPSDEProblem(spsde_v, spsde_f1, spsde_f2, spsde_B, spsde_G1, spsde_G2, (t₀,t₁), Δt, ics)
    prob2 = SPSDEProblem(spsde_v, spsde_f1, spsde_f2, spsde_B, spsde_G1, spsde_G2, (t₀,t₁), Δt, ics; invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity())
    prob3 = SPSDEProblem(spsde_v, spsde_f1, spsde_f2, spsde_B, spsde_G1, spsde_G2, (t₀,t₁), Δt, ics...)
    prob4 = SPSDEProblem(spsde_v, spsde_f1, spsde_f2, spsde_B, spsde_G1, spsde_G2, (t₀,t₁), Δt, ics...; invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity())

    @test prob1 == prob
    @test prob2 == prob
    @test prob3 == prob
    @test prob4 == prob
    
end