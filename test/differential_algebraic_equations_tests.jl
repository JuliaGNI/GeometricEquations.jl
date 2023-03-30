
using GeometricEquations
using GeometricEquations: check_parameters, parameter_types, symplectic_matrix
using Test

include("functions.jl")
include("initial_conditions.jl")


@testset "$(rpad("Differential Algebraic Equations (DAE)",80))" begin

    dae_eqs  = (dae_v, dae_u, dae_ϕ, nothing, nothing, dae_v)
    dae_eqs1 = (dae_v, dae_u, dae_ϕ, nothing, nothing)
    dae_eqs2 = (dae_v, dae_u, dae_ϕ)

    dae = DAE(dae_eqs..., NullInvariants(), NullParameters(), NullPeriodicity())

    @test invariants(dae) == NullInvariants()
    @test parameters(dae) == NullParameters()
    @test periodicity(dae) == NullPeriodicity()

    @test hasvectorfield(dae) == true
    @test hassecondary(dae) == false
    @test hasinvariants(dae) == false
    @test hasparameters(dae) == false
    @test hasperiodicity(dae) == false

    funcs = functions(dae)
    @test funcs.v == dae_v == dae.v
    @test funcs.u == dae_u == dae.u
    @test funcs.ϕ == dae_ϕ == dae.ϕ
    @test funcs.v̄ == dae_v == dae.v̄
    @test !haskey(funcs, :ū)
    @test !haskey(funcs, :ψ)

    funcs = functions(dae, parameters(dae))
    @test funcs.v != dae_v
    @test funcs.u != dae_u
    @test funcs.ϕ != dae_ϕ
    @test funcs.v̄ != dae_v
    
    for eqs in (dae_eqs1, dae_eqs2)
        dae1 = DAE(eqs...)
        dae2 = DAE(eqs...; invariants=NullInvariants())
        dae3 = DAE(eqs...; parameters=NullParameters())
        dae4 = DAE(eqs...; periodicity=NullPeriodicity())
        
        @test dae == dae1
        @test dae == dae2
        @test dae == dae3
        @test dae == dae4
 
        @test hash(dae) == hash(dae1)
        @test hash(dae) == hash(dae2)
        @test hash(dae) == hash(dae3)
        @test hash(dae) == hash(dae4)
    end


    dae_eqs  = (dae_v, dae_u, dae_ϕ, dae_ū, dae_ψ)
    dae_args = (invariants=(h=pdae_h,), parameters=(a=1,), periodicity=π*ones(2))

    dae = DAE(dae_eqs..., dae_v, dae_args.invariants, parameter_types(dae_args.parameters), dae_args.periodicity)

    @test invariants(dae) == dae_args.invariants
    @test parameters(dae) == parameter_types(dae_args.parameters)
    @test periodicity(dae) == dae_args.periodicity

    @test hasvectorfield(dae) == true
    @test hassecondary(dae) == true
    @test hasinvariants(dae) == true
    @test hasparameters(dae) == true
    @test hasperiodicity(dae) == true

    @test check_parameters(dae, dae_args.parameters) == true

    funcs = functions(dae)
    @test funcs.v == dae_v == dae.v
    @test funcs.u == dae_u == dae.u
    @test funcs.ϕ == dae_ϕ == dae.ϕ
    @test funcs.ū == dae_ū == dae.ū
    @test funcs.ψ == dae_ψ == dae.ψ
    @test funcs.v̄ == dae_v == dae.v̄

    funcs = functions(dae, dae_args.parameters)
    @test funcs.v != dae_v
    @test funcs.u != dae_u
    @test funcs.ϕ != dae_ϕ
    @test funcs.ū != dae_ū
    @test funcs.ψ != dae_ψ
    @test funcs.v̄ != dae_v

    dae1 = DAE(dae_eqs...; v̄=dae_v, invariants=dae_args.invariants, parameters=parameter_types(dae_args.parameters), periodicity=dae_args.periodicity)
    dae2 = DAE(dae_eqs...; invariants=dae_args.invariants, parameters=parameter_types(dae_args.parameters), periodicity=dae_args.periodicity)

    @test dae == dae1
    @test dae == dae2

    @test hash(dae) == hash(dae1)
    @test hash(dae) == hash(dae2)

end


@testset "$(rpad("Partitioned Differential Algebraic Equations (PDAE)",80))" begin

    pdae_eqs  = (pdae_v, pdae_f, pdae_u, pdae_g, pdae_ϕ, nothing, nothing, nothing, pdae_v, pdae_f)
    pdae_eqs1 = (pdae_v, pdae_f, pdae_u, pdae_g, pdae_ϕ, nothing, nothing, nothing)
    pdae_eqs2 = (pdae_v, pdae_f, pdae_u, pdae_g, pdae_ϕ)

    pdae  = PDAE(pdae_eqs..., NullInvariants(), NullParameters(), NullPeriodicity())

    @test invariants(pdae) == NullInvariants()
    @test parameters(pdae) == NullParameters()
    @test periodicity(pdae) == NullPeriodicity()

    @test hasvectorfield(pdae) == true
    @test hassecondary(pdae) == false
    @test hasinvariants(pdae) == false
    @test hasparameters(pdae) == false
    @test hasperiodicity(pdae) == false

    funcs = functions(pdae)
    @test funcs.v == pdae_v == pdae.v
    @test funcs.f == pdae_f == pdae.f
    @test funcs.u == pdae_u == pdae.u
    @test funcs.g == pdae_g == pdae.g
    @test funcs.ϕ == pdae_ϕ == pdae.ϕ
    @test funcs.v̄ == pdae_v == pdae.v̄
    @test funcs.f̄ == pdae_f == pdae.f̄
    
    # @test pdae == similar(pdae, t₀, q₀, p₀, λ₀)
    # @test pdae == similar(pdae, t₀, q₀, p₀)
    # @test pdae == similar(pdae, q₀, p₀)

    for eqs in (pdae_eqs1, pdae_eqs2)
        pdae1 = PDAE(eqs...)
        pdae2 = PDAE(eqs...; invariants=NullInvariants())
        pdae3 = PDAE(eqs...; parameters=NullParameters())
        pdae4 = PDAE(eqs...; periodicity=NullPeriodicity())

        @test pdae == pdae1
        @test pdae == pdae2
        @test pdae == pdae3
        @test pdae == pdae4

        @test hash(pdae) == hash(pdae1)
        @test hash(pdae) == hash(pdae2)
        @test hash(pdae) == hash(pdae3)
        @test hash(pdae) == hash(pdae4)
    end


    pdae_eqs  = (pdae_v, pdae_f, pdae_u, pdae_g, pdae_ϕ, pdae_u, pdae_g, pdae_ψ)
    pdae_args = (invariants=(h=pdae_h,), parameters=(a=1,), periodicity=π*ones(1))

    pdae  = PDAE(pdae_eqs..., pdae_v, pdae_f, pdae_args.invariants, parameter_types(pdae_args.parameters), pdae_args.periodicity)

    @test invariants(pdae) == pdae_args.invariants
    @test parameters(pdae) == parameter_types(pdae_args.parameters)
    @test periodicity(pdae) == pdae_args.periodicity

    @test hasvectorfield(pdae) == true
    @test hassecondary(pdae) == true
    @test hasinvariants(pdae) == true
    @test hasparameters(pdae) == true
    @test hasperiodicity(pdae) == true

    funcs = functions(pdae)
    @test funcs.v == pdae_v == pdae.v
    @test funcs.f == pdae_f == pdae.f
    @test funcs.u == pdae_u == pdae.u
    @test funcs.g == pdae_g == pdae.g
    @test funcs.ϕ == pdae_ϕ == pdae.ϕ
    @test funcs.ū == pdae_u == pdae.ū
    @test funcs.ḡ == pdae_g == pdae.ḡ
    @test funcs.ψ == pdae_ψ == pdae.ψ
    @test funcs.v̄ == pdae_v == pdae.v̄
    @test funcs.f̄ == pdae_f == pdae.f̄
    
    funcs = functions(pdae, pdae_args.parameters)
    @test funcs.v != pdae_v
    @test funcs.f != pdae_f
    @test funcs.u != pdae_u
    @test funcs.g != pdae_g
    @test funcs.ϕ != pdae_ϕ
    @test funcs.ū != pdae_u
    @test funcs.ḡ != pdae_g
    @test funcs.ψ != pdae_ψ
    @test funcs.v̄ != pdae_v
    @test funcs.f̄ != pdae_f
    
    pdae1 = PDAE(pdae_eqs...; v̄=pdae_v, f̄=pdae_f, invariants=pdae_args.invariants, parameters=parameter_types(pdae_args.parameters), periodicity=pdae_args.periodicity)
    pdae2 = PDAE(pdae_eqs...; invariants=pdae_args.invariants, parameters=parameter_types(pdae_args.parameters), periodicity=pdae_args.periodicity)

    @test pdae == pdae1
    @test pdae == pdae2

    @test hash(pdae) == hash(pdae1)
    @test hash(pdae) == hash(pdae2)

end


@testset "$(rpad("Implicit Differential Algebraic Equations (IDAE)",80))" begin

    idae_eqs  = (idae_ϑ, pdae_f, pdae_u, pdae_g, pdae_ϕ, nothing, nothing, nothing, idae_v, pdae_f)
    idae_eqs1 = (idae_ϑ, pdae_f, pdae_u, pdae_g, pdae_ϕ, nothing, nothing, nothing)
    idae_eqs2 = (idae_ϑ, pdae_f, pdae_u, pdae_g, pdae_ϕ)

    idae  = IDAE(idae_eqs..., NullInvariants(), NullParameters(), NullPeriodicity())

    @test invariants(idae) == NullInvariants()
    @test parameters(idae) == NullParameters()
    @test periodicity(idae) == NullPeriodicity()

    @test hasvectorfield(idae) == true
    @test hassecondary(idae) == false
    @test hasinvariants(idae) == false
    @test hasparameters(idae) == false
    @test hasperiodicity(idae) == false

    funcs = functions(idae)
    @test funcs.ϑ == idae_ϑ == idae.ϑ
    @test funcs.f == pdae_f == idae.f
    @test funcs.u == pdae_u == idae.u
    @test funcs.g == pdae_g == idae.g
    @test funcs.ϕ == pdae_ϕ == idae.ϕ
    @test funcs.v̄ == idae_v == idae.v̄
    @test funcs.f̄ == pdae_f == idae.f̄
    

    for eqs in (idae_eqs1, idae_eqs2)
        idae1 = IDAE(eqs...)
        idae2 = IDAE(eqs...; invariants=NullInvariants())
        idae3 = IDAE(eqs...; parameters=NullParameters())
        idae4 = IDAE(eqs...; periodicity=NullPeriodicity())

        @test idae == idae1
        @test idae == idae2
        @test idae == idae3
        @test idae == idae4

        @test hash(idae) == hash(idae1)
        @test hash(idae) == hash(idae2)
        @test hash(idae) == hash(idae3)
        @test hash(idae) == hash(idae4)
    end


    idae_eqs  = (idae_ϑ, pdae_f, pdae_u, pdae_g, pdae_ϕ, pdae_u, pdae_g, pdae_ψ)
    idae_args = (invariants=(h=pdae_h,), parameters=(a=1,), periodicity=π*ones(1))

    idae  = IDAE(idae_eqs..., idae_v, pdae_f, idae_args.invariants, parameter_types(idae_args.parameters), idae_args.periodicity)

    @test invariants(idae) == idae_args.invariants
    @test parameters(idae) == parameter_types(idae_args.parameters)
    @test periodicity(idae) == idae_args.periodicity

    @test hasvectorfield(idae) == true
    @test hassecondary(idae) == true
    @test hasinvariants(idae) == true
    @test hasparameters(idae) == true
    @test hasperiodicity(idae) == true

    funcs = functions(idae)
    @test funcs.ϑ == idae_ϑ == idae.ϑ
    @test funcs.f == pdae_f == idae.f
    @test funcs.u == pdae_u == idae.u
    @test funcs.g == pdae_g == idae.g
    @test funcs.ϕ == pdae_ϕ == idae.ϕ
    @test funcs.ū == pdae_u == idae.ū
    @test funcs.ḡ == pdae_g == idae.ḡ
    @test funcs.ψ == pdae_ψ == idae.ψ
    @test funcs.v̄ == idae_v == idae.v̄
    @test funcs.f̄ == pdae_f == idae.f̄
    
    funcs = functions(idae, idae_args.parameters)
    @test funcs.ϑ != idae_ϑ
    @test funcs.f != pdae_f
    @test funcs.u != pdae_u
    @test funcs.g != pdae_g
    @test funcs.ϕ != pdae_ϕ
    @test funcs.ū != pdae_u
    @test funcs.ḡ != pdae_g
    @test funcs.ψ != pdae_ψ
    @test funcs.v̄ != idae_v
    @test funcs.f̄ != pdae_f
    
    idae1 = IDAE(idae_eqs...; v̄=idae_v, f̄=pdae_f, invariants=idae_args.invariants, parameters=parameter_types(idae_args.parameters), periodicity=idae_args.periodicity)
    idae2 = IDAE(idae_eqs...; invariants=idae_args.invariants, parameters=parameter_types(idae_args.parameters), periodicity=idae_args.periodicity)

    @test idae == idae1
    @test idae == idae2

    @test hash(idae) == hash(idae1)
    @test hash(idae) == hash(idae2)

end


@testset "$(rpad("Hamiltonian Differential Algebraic Equations (HDAE)",80))" begin

    hdae_eqs  = (pdae_v, pdae_f, pdae_u, pdae_g, pdae_ϕ, nothing, nothing, nothing, pdae_v, pdae_f, pdae_h)
    hdae_eqs1 = (pdae_v, pdae_f, pdae_u, pdae_g, pdae_ϕ, nothing, nothing, nothing, pdae_h)
    hdae_eqs2 = (pdae_v, pdae_f, pdae_u, pdae_g, pdae_ϕ, pdae_h)

    hdae = HDAE(hdae_eqs..., NullInvariants(), NullParameters(), NullPeriodicity())

    @test hasvectorfield(hdae) == true
    @test hashamiltonian(hdae) == true
    @test hassecondary(hdae) == false
    @test hasinvariants(hdae) == false
    @test hasparameters(hdae) == false
    @test hasperiodicity(hdae) == false

    funcs = functions(hdae)
    @test funcs.v == pdae_v == hdae.v
    @test funcs.f == pdae_f == hdae.f
    @test funcs.u == pdae_u == hdae.u
    @test funcs.g == pdae_g == hdae.g
    @test funcs.ϕ == pdae_ϕ == hdae.ϕ
    @test funcs.v̄ == pdae_v == hdae.v̄
    @test funcs.f̄ == pdae_f == hdae.f̄
    @test funcs.h == pdae_h == hdae.hamiltonian

    for eqs in (hdae_eqs1, hdae_eqs2)
        hdae1 = HDAE(eqs...)
        hdae2 = HDAE(eqs...; invariants=NullInvariants())
        hdae3 = HDAE(eqs...; parameters=NullParameters())
        hdae4 = HDAE(eqs...; periodicity=NullPeriodicity())
        
        @test hdae == hdae1
        @test hdae == hdae2
        @test hdae == hdae3
        @test hdae == hdae4
 
        @test hash(hdae) == hash(hdae1)
        @test hash(hdae) == hash(hdae2)
        @test hash(hdae) == hash(hdae3)
        @test hash(hdae) == hash(hdae4)
    end


    hdae_eqs  = (pdae_v, pdae_f, pdae_u, pdae_g, pdae_ϕ, pdae_u, pdae_g, pdae_ψ, pdae_v, pdae_f, pdae_h)
    hdae_eqs1 = (pdae_v, pdae_f, pdae_u, pdae_g, pdae_ϕ, pdae_u, pdae_g, pdae_ψ, pdae_h)
    hdae_args = (invariants=(h=pdae_h,), parameters=(a=1,), periodicity=π*ones(1))

    hdae = HDAE(hdae_eqs..., hdae_args.invariants, parameter_types(hdae_args.parameters), hdae_args.periodicity)

    @test hasvectorfield(hdae) == true
    @test hashamiltonian(hdae) == true
    @test hassecondary(hdae) == true
    @test hasinvariants(hdae) == true
    @test hasparameters(hdae) == true
    @test hasperiodicity(hdae) == true

    funcs = functions(hdae)
    @test funcs.v == pdae_v == hdae.v
    @test funcs.f == pdae_f == hdae.f
    @test funcs.u == pdae_u == hdae.u
    @test funcs.g == pdae_g == hdae.g
    @test funcs.ϕ == pdae_ϕ == hdae.ϕ
    @test funcs.ū == pdae_u == hdae.ū
    @test funcs.ḡ == pdae_g == hdae.ḡ
    @test funcs.ψ == pdae_ψ == hdae.ψ
    @test funcs.v̄ == pdae_v == hdae.v̄
    @test funcs.f̄ == pdae_f == hdae.f̄
    @test funcs.h == pdae_h == hdae.hamiltonian

    funcs = functions(hdae, hdae_args.parameters)
    @test funcs.v != pdae_v
    @test funcs.f != pdae_f
    @test funcs.u != pdae_u
    @test funcs.g != pdae_g
    @test funcs.ϕ != pdae_ϕ
    @test funcs.ū != pdae_u
    @test funcs.ḡ != pdae_g
    @test funcs.ψ != pdae_ψ
    @test funcs.v̄ != pdae_v
    @test funcs.f̄ != pdae_f
    @test funcs.h != pdae_h

    hdae1 = HDAE(hdae_eqs1...; v̄=pdae_v, f̄=pdae_f, invariants=hdae_args.invariants, parameters=parameter_types(hdae_args.parameters), periodicity=hdae_args.periodicity)
    hdae2 = HDAE(hdae_eqs1...; invariants=hdae_args.invariants, parameters=parameter_types(hdae_args.parameters), periodicity=hdae_args.periodicity)

    @test hdae == hdae1
    @test hdae == hdae2

    @test hash(hdae) == hash(hdae1)
    @test hash(hdae) == hash(hdae2)

end


@testset "$(rpad("Variational Differential Algebraic Equations (LDAE)",80))" begin

    ldae_eqs  = (iode_ϑ, iode_f, iode_u, iode_g, pdae_ϕ, nothing, nothing, nothing, lode_ω, ldae_v, iode_f, lode_l)
    ldae_eqs1 = (iode_ϑ, iode_f, iode_u, iode_g, pdae_ϕ, nothing, nothing, nothing, lode_ω, lode_l)
    ldae_eqs2 = (iode_ϑ, iode_f, iode_u, iode_g, pdae_ϕ, lode_ω, lode_l)

    ldae = LDAE(ldae_eqs..., NullInvariants(), NullParameters(), NullPeriodicity())

    @test hassecondary(ldae) == false
    @test hasinvariants(ldae) == false
    @test hasparameters(ldae) == false
    @test hasperiodicity(ldae) == false

    funcs = functions(ldae)
    @test funcs.ϑ == iode_ϑ == ldae.ϑ
    @test funcs.f == iode_f == ldae.f
    @test funcs.u == iode_u == ldae.u
    @test funcs.g == iode_g == ldae.g
    @test funcs.ϕ == pdae_ϕ == ldae.ϕ
    @test funcs.ω == lode_ω == ldae.ω
    @test funcs.v̄ == ldae_v == ldae.v̄
    @test funcs.f̄ == iode_f == ldae.f̄
    @test funcs.l == lode_l == ldae.lagrangian

    for eqs in (ldae_eqs1, ldae_eqs2)
        ldae1 = LDAE(eqs...)
        ldae2 = LDAE(eqs...; invariants=NullInvariants())
        ldae3 = LDAE(eqs...; parameters=NullParameters())
        ldae4 = LDAE(eqs...; periodicity=NullPeriodicity())

        @test ldae.ϑ == ldae1.ϑ
        @test ldae.f == ldae1.f
        @test ldae.u == ldae1.u
        @test ldae.g == ldae1.g
        @test ldae.ϕ == ldae1.ϕ
        @test ldae.ω == ldae1.ω
        @test ldae.v̄ == ldae1.v̄
        @test ldae.f̄ == ldae1.f̄

        @test ldae == ldae1
        @test ldae == ldae2
        @test ldae == ldae3
        @test ldae == ldae4

        @test hash(ldae) == hash(ldae1)
        @test hash(ldae) == hash(ldae2)
        @test hash(ldae) == hash(ldae3)
        @test hash(ldae) == hash(ldae4)
    end


    ldae_eqs  = (iode_ϑ, iode_f, iode_u, iode_g, pdae_ϕ, pdae_u, pdae_g, pdae_ψ)
    ldae_args = (invariants=(h=iode_h,), parameters=(a=1,), periodicity=π*ones(1))

    ldae = LDAE(ldae_eqs..., lode_ω, iode_v, iode_f, lode_l, ldae_args.invariants, parameter_types(ldae_args.parameters), ldae_args.periodicity)

    @test hassecondary(ldae) == true
    @test hasinvariants(ldae) == true
    @test hasparameters(ldae) == true
    @test hasperiodicity(ldae) == true

    funcs = functions(ldae)
    @test funcs.ϑ == iode_ϑ == ldae.ϑ
    @test funcs.f == iode_f == ldae.f
    @test funcs.u == iode_u == ldae.u
    @test funcs.g == iode_g == ldae.g
    @test funcs.ϕ == pdae_ϕ == ldae.ϕ
    @test funcs.ω == lode_ω == ldae.ω
    @test funcs.v̄ == iode_v == ldae.v̄
    @test funcs.f̄ == iode_f == ldae.f̄
    @test funcs.l == lode_l == ldae.lagrangian

    funcs = functions(ldae, ldae_args.parameters)
    @test funcs.ϑ != iode_ϑ
    @test funcs.f != iode_f
    @test funcs.u != iode_u
    @test funcs.g != iode_g
    @test funcs.ϕ != pdae_ϕ
    @test funcs.ω != lode_ω
    @test funcs.v̄ != iode_v
    @test funcs.f̄ != iode_f
    @test funcs.l != lode_l

    # ldae1 = LDAE(ldae_eqs..., lode_l, lode_ω, t₀, [q₀], [p₀], [λ₀], [λ₀]; v̄=iode_v, ldae_args...)
    # ldae2 = LDAE(ldae_eqs..., lode_l, lode_ω, t₀, [q₀], [p₀], [λ₀]; v̄=iode_v, ldae_args...)
    # ldae3 = LDAE(ldae_eqs..., lode_l, lode_ω, [q₀], [p₀], [λ₀], [λ₀]; v̄=iode_v, ldae_args...)
    # ldae4 = LDAE(ldae_eqs..., lode_l, lode_ω, [q₀], [p₀], [λ₀]; v̄=iode_v, ldae_args...)
    # ldae5 = LDAE(ldae_eqs..., lode_l, lode_ω, t₀, q₀, p₀, λ₀, λ₀; v̄=iode_v, ldae_args...)
    # ldae6 = LDAE(ldae_eqs..., lode_l, lode_ω, t₀, q₀, p₀, λ₀; v̄=iode_v, ldae_args...)
    # ldae7 = LDAE(ldae_eqs..., lode_l, lode_ω, q₀, p₀, λ₀, λ₀; v̄=iode_v, ldae_args...)
    # ldae8 = LDAE(ldae_eqs..., lode_l, lode_ω, q₀, p₀, λ₀; v̄=iode_v, ldae_args...)

    # @test ldae == ldae1
    # @test ldae == ldae2
    # @test ldae == ldae3
    # @test ldae == ldae4
    # @test ldae == ldae5
    # @test ldae == ldae6
    # @test ldae == ldae7
    # @test ldae == ldae8

    # @test hash(ldae) == hash(ldae1)
    # @test hash(ldae) == hash(ldae2)
    # @test hash(ldae) == hash(ldae3)
    # @test hash(ldae) == hash(ldae4)
    # @test hash(ldae) == hash(ldae5)
    # @test hash(ldae) == hash(ldae6)
    # @test hash(ldae) == hash(ldae7)
    # @test hash(ldae) == hash(ldae8)
    
end


@testset "$(rpad("Split Partitioned Differential Algebraic Equations (SPDAE)",80))" begin

    spdae_v = (pdae_v, pdae_u, pdae_u)
    spdae_f = (pdae_f, pdae_g, pdae_g)
    spdae_eqs = (spdae_v, spdae_f, pdae_ϕ, pdae_ψ)

    spdae = SPDAE(spdae_eqs..., t₀, [q₀], [p₀], [λ₀], [λ₀], NullInvariants(), NullParameters(), NullPeriodicity())

    @test ndims(spdae) == 1
    @test nsamples(spdae) == 1
    @test nconstraints(spdae) == 1

    @test periodicity(spdae) == zero(q₀)
    @test initial_conditions(spdae) == (t₀, [q₀], [p₀], [λ₀], [λ₀])

    @test hasinvariants(spdae) == false
    @test hasparameters(spdae) == false
    @test hasperiodicity(spdae) == false

    funcs = functions(spdae)
    @test funcs.v == spdae_v == spdae.v
    @test funcs.f == spdae_f == spdae.f
    @test funcs.ϕ ==  pdae_ϕ == spdae.ϕ
    @test funcs.ψ ==  pdae_ψ == spdae.ψ

    @test spdae == similar(spdae, t₀, q₀, p₀, λ₀)
    @test spdae == similar(spdae, t₀, q₀, p₀)
    @test spdae == similar(spdae, q₀, p₀)

    # spdae1 = SPDAE(spdae_eqs..., t₀, [q₀], [p₀], [λ₀], [λ₀])
    # spdae2 = SPDAE(spdae_eqs..., t₀, [q₀], [p₀], [λ₀])
    # spdae3 = SPDAE(spdae_eqs..., [q₀], [p₀], [λ₀], [λ₀])
    # spdae4 = SPDAE(spdae_eqs..., [q₀], [p₀], [λ₀])
    spdae5 = SPDAE(spdae_eqs..., t₀, q₀, p₀, λ₀, λ₀)
    spdae6 = SPDAE(spdae_eqs..., t₀, q₀, p₀, λ₀)
    spdae7 = SPDAE(spdae_eqs..., q₀, p₀, λ₀, λ₀)
    spdae8 = SPDAE(spdae_eqs..., q₀, p₀, λ₀)

    # @test spdae == spdae1
    # @test spdae == spdae2
    # @test spdae == spdae3
    # @test spdae == spdae4
    @test spdae == spdae5
    @test spdae == spdae6
    @test spdae == spdae7
    @test spdae == spdae8

    # @test hash(spdae) == hash(spdae1)
    # @test hash(spdae) == hash(spdae2)
    # @test hash(spdae) == hash(spdae3)
    # @test hash(spdae) == hash(spdae4)
    @test hash(spdae) == hash(spdae5)
    @test hash(spdae) == hash(spdae6)
    @test hash(spdae) == hash(spdae7)
    @test hash(spdae) == hash(spdae8)

end
