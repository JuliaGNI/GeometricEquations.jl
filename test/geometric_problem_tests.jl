using GeometricEquations
using GeometricEquations: datatype, timetype, arrtype, equtype
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

    @test initial_conditions(prob) == (t₀, x₀)
    @test parameters(prob) == NullParameters()
    @test equation(prob) == ode
    @test nsamples(prob) == 1

    @test functions(prob) == functions(ode)
    @test functions(prob) == NamedTuple{(:v,)}((ode_v,))

    prob1 = GeometricProblem(ode, (t₀,t₁), Δt, ics, NullParameters())
    prob2 = GeometricProblem(ode, (t₀,t₁), Δt, ics, nothing)

    @test prob == prob1
    @test prob == prob2

    # @test ode == similar(ode, t₀, x₀)
    # @test ode == similar(ode, x₀)

end


@testset "$(rpad("ODE Problem",80))" begin

    ics  = (q=x₀,)
    ode  = ODE(ode_v)
    prob = GeometricProblem(ode, (t₀,t₁), Δt, ics)

    @test typeof(prob) <: GeometricProblem
    @test typeof(prob) <: ODEProblem
    @test equtype(prob) == ODE

    prob1 = ODEProblem(ode_v, (t₀,t₁), Δt, ics)
    prob2 = ODEProblem(ode_v, (t₀,t₁), Δt, ics, NullParameters())
    prob3 = ODEProblem(ode_v, (t₀,t₁), Δt, ics, NullParameters(); invariants=NullInvariants(), periodicity=NullPeriodicity())
    prob4 = ODEProblem(ode_v, (t₀,t₁), Δt, ics...)
    prob5 = ODEProblem(ode_v, (t₀,t₁), Δt, ics..., NullParameters())
    prob6 = ODEProblem(ode_v, (t₀,t₁), Δt, ics..., NullParameters(); invariants=NullInvariants(), periodicity=NullPeriodicity())

    @test prob1 == prob
    @test prob2 == prob
    @test prob3 == prob
    @test prob4 == prob
    @test prob5 == prob
    @test prob6 == prob
    
end


@testset "$(rpad("SODE Problem",80))" begin

    ics  = (q=x₀,)
    sode = SODE(sode_v)
    prob = GeometricProblem(sode, (t₀,t₁), Δt, ics)

    @test typeof(prob) <: GeometricProblem
    @test typeof(prob) <: SODEProblem
    @test equtype(prob) == SODE

    prob1 = SODEProblem(sode_v, (t₀,t₁), Δt, ics)
    prob2 = SODEProblem(sode_v, (t₀,t₁), Δt, ics, NullParameters())
    prob3 = SODEProblem(sode_v, (t₀,t₁), Δt, ics, NullParameters(); invariants=NullInvariants(), periodicity=NullPeriodicity())
    prob4 = SODEProblem(sode_v, (t₀,t₁), Δt, ics...)
    prob5 = SODEProblem(sode_v, (t₀,t₁), Δt, ics..., NullParameters())
    prob6 = SODEProblem(sode_v, (t₀,t₁), Δt, ics..., NullParameters(); invariants=NullInvariants(), periodicity=NullPeriodicity())

    @test prob1 == prob
    @test prob2 == prob
    @test prob3 == prob
    @test prob4 == prob
    @test prob5 == prob
    @test prob6 == prob
    

    sode = SODE(sode_v, sode_q)
    prob = GeometricProblem(sode, (t₀,t₁), Δt, ics)

    @test typeof(prob) <: GeometricProblem
    @test typeof(prob) <: SODEProblem
    @test equtype(prob) == SODE

    prob1 = SODEProblem(sode_v, sode_q, (t₀,t₁), Δt, ics)
    prob2 = SODEProblem(sode_v, sode_q, (t₀,t₁), Δt, ics, NullParameters())
    prob3 = SODEProblem(sode_v, sode_q, (t₀,t₁), Δt, ics, NullParameters(); invariants=NullInvariants(), periodicity=NullPeriodicity())
    prob4 = SODEProblem(sode_v, sode_q, (t₀,t₁), Δt, ics...)
    prob5 = SODEProblem(sode_v, sode_q, (t₀,t₁), Δt, ics..., NullParameters())
    prob6 = SODEProblem(sode_v, sode_q, (t₀,t₁), Δt, ics..., NullParameters(); invariants=NullInvariants(), periodicity=NullPeriodicity())

    @test prob1 == prob
    @test prob2 == prob
    @test prob3 == prob
    @test prob4 == prob
    @test prob5 == prob
    @test prob6 == prob
    
end


@testset "$(rpad("PODE Problem",80))" begin

    ics  = (q=q₀, p=p₀)
    pode = PODE(pode_eqs...)
    prob = GeometricProblem(pode, (t₀,t₁), Δt, ics)

    @test typeof(prob) <: GeometricProblem
    @test typeof(prob) <: PODEProblem
    @test equtype(prob) == PODE

    prob1 = PODEProblem(pode_eqs..., (t₀,t₁), Δt, ics)
    prob2 = PODEProblem(pode_eqs..., (t₀,t₁), Δt, ics, NullParameters())
    prob3 = PODEProblem(pode_eqs..., (t₀,t₁), Δt, ics, NullParameters(); invariants=NullInvariants(), periodicity=NullPeriodicity())
    prob4 = PODEProblem(pode_eqs..., (t₀,t₁), Δt, ics...)
    prob5 = PODEProblem(pode_eqs..., (t₀,t₁), Δt, ics..., NullParameters())
    prob6 = PODEProblem(pode_eqs..., (t₀,t₁), Δt, ics..., NullParameters(); invariants=NullInvariants(), periodicity=NullPeriodicity())

    @test prob1 == prob
    @test prob2 == prob
    @test prob3 == prob
    @test prob4 == prob
    @test prob5 == prob
    @test prob6 == prob
    
end


@testset "$(rpad("IODE Problem",80))" begin

    ics  = (q=q₀, p=p₀, λ=λ₀)
    iode = IODE(iode_eqs...)
    prob = GeometricProblem(iode, (t₀,t₁), Δt, ics)

    @test typeof(prob) <: GeometricProblem
    @test typeof(prob) <: IODEProblem
    @test equtype(prob) == IODE

    prob1 = IODEProblem(iode_eqs..., (t₀,t₁), Δt, ics)
    prob2 = IODEProblem(iode_eqs..., (t₀,t₁), Δt, ics, NullParameters())
    prob3 = IODEProblem(iode_eqs..., (t₀,t₁), Δt, ics, NullParameters(); invariants=NullInvariants(), periodicity=NullPeriodicity())
    prob4 = IODEProblem(iode_eqs..., (t₀,t₁), Δt, ics...)
    prob5 = IODEProblem(iode_eqs..., (t₀,t₁), Δt, ics..., NullParameters())
    prob6 = IODEProblem(iode_eqs..., (t₀,t₁), Δt, ics..., NullParameters(); invariants=NullInvariants(), periodicity=NullPeriodicity())

    @test prob1 == prob
    @test prob2 == prob
    @test prob3 == prob
    @test prob4 == prob
    @test prob5 == prob
    @test prob6 == prob
    
end


@testset "$(rpad("HODE Problem",80))" begin

    ics  = (q=q₀, p=p₀)
    hode = HODE(hode_eqs...)
    prob = GeometricProblem(hode, (t₀,t₁), Δt, ics)

    @test typeof(prob) <: GeometricProblem
    @test typeof(prob) <: HODEProblem
    @test equtype(prob) == HODE

    prob1 = HODEProblem(hode_eqs..., (t₀,t₁), Δt, ics)
    prob2 = HODEProblem(hode_eqs..., (t₀,t₁), Δt, ics, NullParameters())
    prob3 = HODEProblem(hode_eqs..., (t₀,t₁), Δt, ics, NullParameters(); invariants=NullInvariants(), periodicity=NullPeriodicity())
    prob4 = HODEProblem(hode_eqs..., (t₀,t₁), Δt, ics...)
    prob5 = HODEProblem(hode_eqs..., (t₀,t₁), Δt, ics..., NullParameters())
    prob6 = HODEProblem(hode_eqs..., (t₀,t₁), Δt, ics..., NullParameters(); invariants=NullInvariants(), periodicity=NullPeriodicity())

    @test prob1 == prob
    @test prob2 == prob
    @test prob3 == prob
    @test prob4 == prob
    @test prob5 == prob
    @test prob6 == prob
    
end


@testset "$(rpad("LODE Problem",80))" begin

    ics  = (q=q₀, p=p₀, λ=λ₀)
    lode = LODE(lode_eqs...)
    prob = GeometricProblem(lode, (t₀,t₁), Δt, ics)

    @test typeof(prob) <: GeometricProblem
    @test typeof(prob) <: LODEProblem
    @test equtype(prob) == LODE

    prob1 = LODEProblem(lode_eqs..., (t₀,t₁), Δt, ics)
    prob2 = LODEProblem(lode_eqs..., (t₀,t₁), Δt, ics, NullParameters())
    prob3 = LODEProblem(lode_eqs..., (t₀,t₁), Δt, ics, NullParameters(); invariants=NullInvariants(), periodicity=NullPeriodicity())
    prob4 = LODEProblem(lode_eqs..., (t₀,t₁), Δt, ics...)
    prob5 = LODEProblem(lode_eqs..., (t₀,t₁), Δt, ics..., NullParameters())
    prob6 = LODEProblem(lode_eqs..., (t₀,t₁), Δt, ics..., NullParameters(); invariants=NullInvariants(), periodicity=NullPeriodicity())

    @test prob1 == prob
    @test prob2 == prob
    @test prob3 == prob
    @test prob4 == prob
    @test prob5 == prob
    @test prob6 == prob
    
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
    prob2 = DAEProblem(eqs..., (t₀,t₁), Δt, ics; parameters=NullParameters())
    prob3 = DAEProblem(eqs..., (t₀,t₁), Δt, ics; invariants=NullInvariants())
    prob4 = DAEProblem(eqs..., (t₀,t₁), Δt, ics; periodicity=NullPeriodicity())

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

    prob1 = DAEProblem(eqs..., (t₀,t₁), Δt, ics...)
    prob2 = DAEProblem(eqs..., (t₀,t₁), Δt, ics...; parameters=NullParameters())
    prob3 = DAEProblem(eqs..., (t₀,t₁), Δt, ics...; invariants=NullInvariants())
    prob4 = DAEProblem(eqs..., (t₀,t₁), Δt, ics...; periodicity=NullPeriodicity())

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
    prob2 = PDAEProblem(eqs..., (t₀,t₁), Δt, ics; parameters=NullParameters())
    prob3 = PDAEProblem(eqs..., (t₀,t₁), Δt, ics; invariants=NullInvariants())
    prob4 = PDAEProblem(eqs..., (t₀,t₁), Δt, ics; periodicity=NullPeriodicity())

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

    prob1 = PDAEProblem(eqs..., (t₀,t₁), Δt, ics...)
    prob2 = PDAEProblem(eqs..., (t₀,t₁), Δt, ics...; parameters=NullParameters())
    prob3 = PDAEProblem(eqs..., (t₀,t₁), Δt, ics...; invariants=NullInvariants())
    prob4 = PDAEProblem(eqs..., (t₀,t₁), Δt, ics...; periodicity=NullPeriodicity())

    @test prob1 == prob
    @test prob2 == prob
    @test prob3 == prob
    @test prob4 == prob
    
end


@testset "$(rpad("HDAE Problem",80))" begin

    eqs  = (pdae_v, pdae_f, pdae_u, pdae_g, pdae_ϕ, nothing, nothing, nothing, hdae_ω, pdae_h)
    ics  = (q=q₀, p=p₀, λ=λ₀)
    hdae = HDAE(eqs...)
    prob = GeometricProblem(hdae, (t₀,t₁), Δt, ics)

    @test typeof(prob) <: GeometricProblem
    @test typeof(prob) <: HDAEProblem
    @test equtype(prob) == HDAE

    prob1 = HDAEProblem(eqs..., (t₀,t₁), Δt, ics)
    prob2 = HDAEProblem(eqs..., (t₀,t₁), Δt, ics; parameters=NullParameters())
    prob3 = HDAEProblem(eqs..., (t₀,t₁), Δt, ics; invariants=NullInvariants())
    prob4 = HDAEProblem(eqs..., (t₀,t₁), Δt, ics; periodicity=NullPeriodicity())

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


    eqs  = (pdae_v, pdae_f, pdae_u, pdae_g, pdae_ϕ, pdae_u, pdae_g, pdae_ψ, hdae_ω, pdae_h)
    ics  = (q=q₀, p=p₀, λ=λ₀, μ=λ₀)
    hdae = HDAE(eqs...)
    prob = GeometricProblem(hdae, (t₀,t₁), Δt, ics)

    @test typeof(prob) <: GeometricProblem
    @test typeof(prob) <: HDAEProblem
    @test equtype(prob) == HDAE

    prob1 = HDAEProblem(eqs..., (t₀,t₁), Δt, ics...)
    prob2 = HDAEProblem(eqs..., (t₀,t₁), Δt, ics...; parameters=NullParameters())
    prob3 = HDAEProblem(eqs..., (t₀,t₁), Δt, ics...; invariants=NullInvariants())
    prob4 = HDAEProblem(eqs..., (t₀,t₁), Δt, ics...; periodicity=NullPeriodicity())

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
    prob2 = IDAEProblem(eqs..., (t₀,t₁), Δt, ics; parameters=NullParameters())
    prob3 = IDAEProblem(eqs..., (t₀,t₁), Δt, ics; invariants=NullInvariants())
    prob4 = IDAEProblem(eqs..., (t₀,t₁), Δt, ics; periodicity=NullPeriodicity())

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

    prob1 = IDAEProblem(eqs..., (t₀,t₁), Δt, ics...)
    prob2 = IDAEProblem(eqs..., (t₀,t₁), Δt, ics...; parameters=NullParameters())
    prob3 = IDAEProblem(eqs..., (t₀,t₁), Δt, ics...; invariants=NullInvariants())
    prob4 = IDAEProblem(eqs..., (t₀,t₁), Δt, ics...; periodicity=NullPeriodicity())

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
    prob2 = LDAEProblem(eqs..., (t₀,t₁), Δt, ics; parameters=NullParameters())
    prob3 = LDAEProblem(eqs..., (t₀,t₁), Δt, ics; invariants=NullInvariants())
    prob4 = LDAEProblem(eqs..., (t₀,t₁), Δt, ics; periodicity=NullPeriodicity())

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

    prob1 = LDAEProblem(eqs..., (t₀,t₁), Δt, ics...)
    prob2 = LDAEProblem(eqs..., (t₀,t₁), Δt, ics...; parameters=NullParameters())
    prob3 = LDAEProblem(eqs..., (t₀,t₁), Δt, ics...; invariants=NullInvariants())
    prob4 = LDAEProblem(eqs..., (t₀,t₁), Δt, ics...; periodicity=NullPeriodicity())

    @test prob1 == prob
    @test prob2 == prob
    @test prob3 == prob
    @test prob4 == prob
    
end
