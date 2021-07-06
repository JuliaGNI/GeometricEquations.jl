
using GeometricEquations
using GeometricEquations: symplectic_matrix
using Test

include("functions.jl")
include("initial_conditions.jl")


@testset "$(rpad("Differential Algebraic Equations (DAE)",80))" begin

    dae_eqs  = (dae_v, dae_u, nothing, dae_ϕ, nothing, dae_v)
    dae_eqs1 = (dae_v, dae_u, nothing, dae_ϕ, nothing)
    dae_eqs2 = (dae_v, dae_u, dae_ϕ)

    dae = DAE(dae_eqs..., t₀, [x₀], [λ₀], [λ₀], NullInvariants(), NullParameters(), nothing)

    @test axes(dae) == axes(x₀)
    @test ndims(dae) == 2
    @test nsamples(dae) == 1
    @test nconstraints(dae) == 1

    @test periodicity(dae) == zero(x₀)
    @test initial_conditions(dae) == (t₀, [x₀], [λ₀], [λ₀])

    @test hassecondary(dae) == false
    @test hasinvariants(dae) == false
    @test hasparameters(dae) == false
    @test hasperiodicity(dae) == false

    functions = get_functions(dae)
    @test functions.v == dae_v == dae.v
    @test functions.u == dae_u == dae.u
    @test functions.ϕ == dae_ϕ == dae.ϕ
    @test functions.v̄ == dae_v == dae.v̄

    @test dae == similar(dae, t₀, x₀, λ₀)
    @test dae == similar(dae, x₀, λ₀)

    for eqs in (dae_eqs1, dae_eqs2)
        dae1 = DAE(eqs..., t₀, [x₀], [λ₀], [λ₀])
        dae2 = DAE(eqs..., t₀, [x₀], [λ₀])
        dae3 = DAE(eqs..., [x₀], [λ₀], [λ₀])
        dae4 = DAE(eqs..., [x₀], [λ₀])
        dae5 = DAE(eqs..., t₀, x₀, λ₀, λ₀)
        dae6 = DAE(eqs..., t₀, x₀, λ₀)
        dae7 = DAE(eqs..., x₀, λ₀, λ₀)
        dae8 = DAE(eqs..., x₀, λ₀)
        
        @test dae == dae1
        @test dae == dae2
        @test dae == dae3
        @test dae == dae4
        @test dae == dae5
        @test dae == dae6
        @test dae == dae7
        @test dae == dae8

        @test hash(dae) == hash(dae1)
        @test hash(dae) == hash(dae2)
        @test hash(dae) == hash(dae3)
        @test hash(dae) == hash(dae4)
        @test hash(dae) == hash(dae5)
        @test hash(dae) == hash(dae6)
        @test hash(dae) == hash(dae7)
        @test hash(dae) == hash(dae8)
    end


    dae_eqs  = (dae_v, dae_u, dae_ū, dae_ϕ, dae_ψ)
    dae_args = (invariants=(h=pdae_h,), parameters=(a=1,), periodicity=π*ones(2))

    dae = DAE(dae_eqs..., dae_v, t₀, [x₀], [λ₀], [λ₀], dae_args.invariants, dae_args.parameters, dae_args.periodicity)

    @test axes(dae) == axes(x₀)
    @test ndims(dae) == 2
    @test nsamples(dae) == 1
    @test nconstraints(dae) == 1

    @test periodicity(dae) == dae_args.periodicity
    @test initial_conditions(dae) == (t₀, [x₀], [λ₀], [λ₀])

    @test hassecondary(dae) == true
    @test hasinvariants(dae) == true
    @test hasparameters(dae) == true
    @test hasperiodicity(dae) == true

    functions = get_functions(dae)
    @test functions.v != dae_v == dae.v
    @test functions.u != dae_u == dae.u
    @test functions.ū != dae_ū == dae.ū
    @test functions.ϕ != dae_ϕ == dae.ϕ
    @test functions.ψ != dae_ψ == dae.ψ
    @test functions.v̄ != dae_v == dae.v̄

    @test dae == similar(dae, t₀, x₀, λ₀)
    @test dae == similar(dae, x₀, λ₀)

    dae1 = DAE(dae_eqs..., [x₀], [λ₀]; dae_args...)
    dae2 = DAE(dae_eqs..., [x₀], [λ₀]; dae_args...)
    dae3 = DAE(dae_eqs..., t₀, x₀, λ₀; dae_args...)
    dae4 = DAE(dae_eqs..., x₀, λ₀; dae_args...)

    @test dae == dae1
    @test dae == dae2
    @test dae == dae3
    @test dae == dae4

    @test hash(dae) == hash(dae1)
    @test hash(dae) == hash(dae2)
    @test hash(dae) == hash(dae3)
    @test hash(dae) == hash(dae4)

end


@testset "$(rpad("Partitioned Differential Algebraic Equations (PDAE)",80))" begin

    pdae_eqs  = (pdae_v, pdae_f, pdae_u, pdae_g, pdae_ϕ, nothing, nothing, nothing, pdae_v, pdae_f)
    pdae_eqs1 = (pdae_v, pdae_f, pdae_u, pdae_g, pdae_ϕ, nothing, nothing, nothing)
    pdae_eqs2 = (pdae_v, pdae_f, pdae_u, pdae_g, pdae_ϕ)

    pdae  = PDAE(pdae_eqs..., t₀, [q₀], [p₀], [λ₀], [λ₀], NullInvariants(), NullParameters(), nothing)

    @test axes(pdae) == axes(q₀)
    @test ndims(pdae) == 1
    @test nsamples(pdae) == 1
    @test nconstraints(pdae) == 1

    @test periodicity(pdae) == zero(q₀)
    @test initial_conditions(pdae) == (t₀, [q₀], [p₀], [λ₀], [λ₀])

    @test hassecondary(pdae) == false
    @test hasinvariants(pdae) == false
    @test hasparameters(pdae) == false
    @test hasperiodicity(pdae) == false

    functions = get_functions(pdae)
    @test functions.v == pdae_v == pdae.v
    @test functions.f == pdae_f == pdae.f
    @test functions.u == pdae_u == pdae.u
    @test functions.g == pdae_g == pdae.g
    @test functions.ϕ == pdae_ϕ == pdae.ϕ
    @test functions.v̄ == pdae_v == pdae.v̄
    @test functions.f̄ == pdae_f == pdae.f̄
    
    @test pdae == similar(pdae, t₀, q₀, p₀, λ₀)
    @test pdae == similar(pdae, t₀, q₀, p₀)
    @test pdae == similar(pdae, q₀, p₀)

    for eqs in (pdae_eqs1, pdae_eqs2)
        pdae1 = PDAE(eqs..., t₀, [q₀], [p₀], [λ₀], [λ₀])
        pdae2 = PDAE(eqs..., t₀, [q₀], [p₀], [λ₀])
        pdae3 = PDAE(eqs..., [q₀], [p₀], [λ₀], [λ₀])
        pdae4 = PDAE(eqs..., [q₀], [p₀], [λ₀])
        pdae5 = PDAE(eqs..., t₀, q₀, p₀, λ₀, λ₀)
        pdae6 = PDAE(eqs..., t₀, q₀, p₀, λ₀)
        pdae7 = PDAE(eqs..., q₀, p₀, λ₀, λ₀)
        pdae8 = PDAE(eqs..., q₀, p₀, λ₀)

        @test pdae == pdae1
        @test pdae == pdae2
        @test pdae == pdae3
        @test pdae == pdae4
        @test pdae == pdae5
        @test pdae == pdae6
        @test pdae == pdae7
        @test pdae == pdae8

        @test hash(pdae) == hash(pdae1)
        @test hash(pdae) == hash(pdae2)
        @test hash(pdae) == hash(pdae3)
        @test hash(pdae) == hash(pdae4)
        @test hash(pdae) == hash(pdae5)
        @test hash(pdae) == hash(pdae6)
        @test hash(pdae) == hash(pdae7)
        @test hash(pdae) == hash(pdae8)
    end


    pdae_eqs  = (pdae_v, pdae_f, pdae_u, pdae_g, pdae_ϕ, pdae_u, pdae_g, pdae_ψ)
    pdae_args = (invariants=(h=pdae_h,), parameters=(a=1,), periodicity=π*ones(1))

    pdae  = PDAE(pdae_eqs..., pdae_v, pdae_f, t₀, [q₀], [p₀], [λ₀], [λ₀], pdae_args.invariants, pdae_args.parameters, pdae_args.periodicity)

    @test axes(pdae) == axes(q₀)
    @test ndims(pdae) == 1
    @test nsamples(pdae) == 1
    @test nconstraints(pdae) == 1

    @test periodicity(pdae) == pdae_args.periodicity
    @test initial_conditions(pdae) == (t₀, [q₀], [p₀], [λ₀], [λ₀])

    @test hassecondary(pdae) == true
    @test hasinvariants(pdae) == true
    @test hasparameters(pdae) == true
    @test hasperiodicity(pdae) == true

    functions = get_functions(pdae)
    @test functions.v != pdae_v == pdae.v
    @test functions.f != pdae_f == pdae.f
    @test functions.u != pdae_u == pdae.u
    @test functions.g != pdae_g == pdae.g
    @test functions.ϕ != pdae_ϕ == pdae.ϕ
    @test functions.ū != pdae_u == pdae.ū
    @test functions.ḡ != pdae_g == pdae.ḡ
    @test functions.ψ != pdae_ψ == pdae.ψ
    @test functions.v̄ != pdae_v == pdae.v̄
    @test functions.f̄ != pdae_f == pdae.f̄
    
    @test pdae == similar(pdae, t₀, q₀, p₀, λ₀)
    @test pdae == similar(pdae, t₀, q₀, p₀)
    @test pdae == similar(pdae, q₀, p₀)

    pdae1 = PDAE(pdae_eqs..., t₀, [q₀], [p₀], [λ₀], [λ₀]; pdae_args...)
    pdae2 = PDAE(pdae_eqs..., t₀, [q₀], [p₀], [λ₀]; pdae_args...)
    pdae3 = PDAE(pdae_eqs..., [q₀], [p₀], [λ₀], [λ₀]; pdae_args...)
    pdae4 = PDAE(pdae_eqs..., [q₀], [p₀], [λ₀]; pdae_args...)
    pdae5 = PDAE(pdae_eqs..., t₀, q₀, p₀, λ₀, λ₀; pdae_args...)
    pdae6 = PDAE(pdae_eqs..., t₀, q₀, p₀, λ₀; pdae_args...)
    pdae7 = PDAE(pdae_eqs..., q₀, p₀, λ₀, λ₀; pdae_args...)
    pdae8 = PDAE(pdae_eqs..., q₀, p₀, λ₀; pdae_args...)

    @test pdae == pdae1
    @test pdae == pdae2
    @test pdae == pdae3
    @test pdae == pdae4
    @test pdae == pdae5
    @test pdae == pdae6
    @test pdae == pdae7
    @test pdae == pdae8

    @test hash(pdae) == hash(pdae1)
    @test hash(pdae) == hash(pdae2)
    @test hash(pdae) == hash(pdae3)
    @test hash(pdae) == hash(pdae4)
    @test hash(pdae) == hash(pdae5)
    @test hash(pdae) == hash(pdae6)
    @test hash(pdae) == hash(pdae7)
    @test hash(pdae) == hash(pdae8)

end


@testset "$(rpad("Implicit Differential Algebraic Equations (IDAE)",80))" begin

    idae_eqs  = (pdae_p, pdae_f, pdae_u, pdae_g, pdae_ϕ, nothing, nothing, nothing, pdae_v, pdae_f)
    idae_eqs1 = (pdae_p, pdae_f, pdae_u, pdae_g, pdae_ϕ, nothing, nothing, nothing)
    idae_eqs2 = (pdae_p, pdae_f, pdae_u, pdae_g, pdae_ϕ)

    idae  = IDAE(idae_eqs..., t₀, [q₀], [p₀], [λ₀], [λ₀], NullInvariants(), NullParameters(), nothing)

    @test axes(idae) == axes(q₀)
    @test ndims(idae) == 1
    @test nsamples(idae) == 1
    @test nconstraints(idae) == 1

    @test periodicity(idae) == zero(q₀)
    @test initial_conditions(idae) == (t₀, [q₀], [p₀], [λ₀], [λ₀])

    @test hassecondary(idae) == false
    @test hasinvariants(idae) == false
    @test hasparameters(idae) == false
    @test hasperiodicity(idae) == false

    functions = get_functions(idae)
    @test functions.ϑ == pdae_p == idae.ϑ
    @test functions.f == pdae_f == idae.f
    @test functions.u == pdae_u == idae.u
    @test functions.g == pdae_g == idae.g
    @test functions.ϕ == pdae_ϕ == idae.ϕ
    @test functions.v̄ == pdae_v == idae.v̄
    @test functions.f̄ == pdae_f == idae.f̄
    
    @test idae == similar(idae, t₀, q₀, p₀, λ₀)
    @test idae == similar(idae, t₀, q₀, p₀)
    @test idae == similar(idae, q₀, p₀)

    for eqs in (idae_eqs1, idae_eqs2)
        idae1 = IDAE(eqs..., t₀, [q₀], [p₀], [λ₀], [λ₀]; v̄=pdae_v)
        idae2 = IDAE(eqs..., t₀, [q₀], [p₀], [λ₀]; v̄=pdae_v)
        idae3 = IDAE(eqs..., [q₀], [p₀], [λ₀], [λ₀]; v̄=pdae_v)
        idae4 = IDAE(eqs..., [q₀], [p₀], [λ₀]; v̄=pdae_v)
        idae5 = IDAE(eqs..., t₀, q₀, p₀, λ₀, λ₀; v̄=pdae_v)
        idae6 = IDAE(eqs..., t₀, q₀, p₀, λ₀; v̄=pdae_v)
        idae7 = IDAE(eqs..., q₀, p₀, λ₀, λ₀; v̄=pdae_v)
        idae8 = IDAE(eqs..., q₀, p₀, λ₀; v̄=pdae_v)

        @test idae == idae1
        @test idae == idae2
        @test idae == idae3
        @test idae == idae4
        @test idae == idae5
        @test idae == idae6
        @test idae == idae7
        @test idae == idae8

        @test hash(idae) == hash(idae1)
        @test hash(idae) == hash(idae2)
        @test hash(idae) == hash(idae3)
        @test hash(idae) == hash(idae4)
        @test hash(idae) == hash(idae5)
        @test hash(idae) == hash(idae6)
        @test hash(idae) == hash(idae7)
        @test hash(idae) == hash(idae8)
    end


    idae_eqs  = (pdae_p, pdae_f, pdae_u, pdae_g, pdae_ϕ, pdae_u, pdae_g, pdae_ψ)
    idae_args = (invariants=(h=pdae_h,), parameters=(a=1,), periodicity=π*ones(1))

    idae  = IDAE(idae_eqs..., pdae_v, pdae_f, t₀, [q₀], [p₀], [λ₀], [λ₀], idae_args.invariants, idae_args.parameters, idae_args.periodicity)

    @test axes(idae) == axes(q₀)
    @test ndims(idae) == 1
    @test nsamples(idae) == 1
    @test nconstraints(idae) == 1

    @test periodicity(idae) == idae_args.periodicity
    @test initial_conditions(idae) == (t₀, [q₀], [p₀], [λ₀], [λ₀])

    @test hassecondary(idae) == true
    @test hasinvariants(idae) == true
    @test hasparameters(idae) == true
    @test hasperiodicity(idae) == true

    functions = get_functions(idae)
    @test functions.ϑ != pdae_p == idae.ϑ
    @test functions.f != pdae_f == idae.f
    @test functions.u != pdae_u == idae.u
    @test functions.g != pdae_g == idae.g
    @test functions.ϕ != pdae_ϕ == idae.ϕ
    @test functions.ū != pdae_u == idae.ū
    @test functions.ḡ != pdae_g == idae.ḡ
    @test functions.ψ != pdae_ψ == idae.ψ
    @test functions.v̄ != pdae_v == idae.v̄
    @test functions.f̄ != pdae_f == idae.f̄
    
    @test idae == similar(idae, t₀, q₀, p₀, λ₀)
    @test idae == similar(idae, t₀, q₀, p₀)
    @test idae == similar(idae, q₀, p₀)

    idae1 = IDAE(idae_eqs..., t₀, [q₀], [p₀], [λ₀], [λ₀]; v̄=pdae_v, idae_args...)
    idae2 = IDAE(idae_eqs..., t₀, [q₀], [p₀], [λ₀]; v̄=pdae_v, idae_args...)
    idae3 = IDAE(idae_eqs..., [q₀], [p₀], [λ₀], [λ₀]; v̄=pdae_v, idae_args...)
    idae4 = IDAE(idae_eqs..., [q₀], [p₀], [λ₀]; v̄=pdae_v, idae_args...)
    idae5 = IDAE(idae_eqs..., t₀, q₀, p₀, λ₀, λ₀; v̄=pdae_v, idae_args...)
    idae6 = IDAE(idae_eqs..., t₀, q₀, p₀, λ₀; v̄=pdae_v, idae_args...)
    idae7 = IDAE(idae_eqs..., q₀, p₀, λ₀, λ₀; v̄=pdae_v, idae_args...)
    idae8 = IDAE(idae_eqs..., q₀, p₀, λ₀; v̄=pdae_v, idae_args...)

    @test idae == idae1
    @test idae == idae2
    @test idae == idae3
    @test idae == idae4
    @test idae == idae5
    @test idae == idae6
    @test idae == idae7
    @test idae == idae8

    @test hash(idae) == hash(idae1)
    @test hash(idae) == hash(idae2)
    @test hash(idae) == hash(idae3)
    @test hash(idae) == hash(idae4)
    @test hash(idae) == hash(idae5)
    @test hash(idae) == hash(idae6)
    @test hash(idae) == hash(idae7)
    @test hash(idae) == hash(idae8)

end


@testset "$(rpad("Hamiltonian Differential Algebraic Equations (HDAE)",80))" begin

    hdae_eqs  = (pdae_v, pdae_f, pdae_u, pdae_g, pdae_ϕ, nothing, nothing, nothing, pdae_v, pdae_f, symplectic_matrix)
    hdae_eqs1 = (pdae_v, pdae_f, pdae_u, pdae_g, pdae_ϕ, nothing, nothing, nothing, pdae_h)
    hdae_eqs2 = (pdae_v, pdae_f, pdae_u, pdae_g, pdae_ϕ, pdae_h)

    hdae  = HDAE(hdae_eqs..., t₀, [q₀], [p₀], [λ₀], [λ₀], pdae_h, NullInvariants(), NullParameters(), nothing)

    @test axes(hdae) == axes(q₀)
    @test ndims(hdae) == 1
    @test nsamples(hdae) == 1
    @test nconstraints(hdae) == 1

    @test periodicity(hdae) == zero(q₀)
    @test initial_conditions(hdae) == (t₀, [q₀], [p₀], [λ₀], [λ₀])

    @test hassecondary(hdae) == false
    @test hasinvariants(hdae) == false
    @test hasparameters(hdae) == false
    @test hasperiodicity(hdae) == false

    functions = get_functions(hdae)
    @test functions.v == pdae_v == hdae.v
    @test functions.f == pdae_f == hdae.f
    @test functions.u == pdae_u == hdae.u
    @test functions.g == pdae_g == hdae.g
    @test functions.ϕ == pdae_ϕ == hdae.ϕ
    @test functions.v̄ == pdae_v == hdae.v̄
    @test functions.f̄ == pdae_f == hdae.f̄
    @test functions.h == pdae_h == hdae.hamiltonian

    @test hdae == similar(hdae, t₀, q₀, p₀, λ₀)
    @test hdae == similar(hdae, t₀, q₀, p₀)
    @test hdae == similar(hdae, q₀, p₀)

    for eqs in (hdae_eqs1, hdae_eqs2)
        hdae1 = HDAE(eqs..., t₀, [q₀], [p₀], [λ₀], [λ₀])
        hdae2 = HDAE(eqs..., t₀, [q₀], [p₀], [λ₀])
        hdae3 = HDAE(eqs..., [q₀], [p₀], [λ₀], [λ₀])
        hdae4 = HDAE(eqs..., [q₀], [p₀], [λ₀])
        hdae5 = HDAE(eqs..., t₀, q₀, p₀, λ₀, λ₀)
        hdae6 = HDAE(eqs..., t₀, q₀, p₀, λ₀)
        hdae7 = HDAE(eqs..., q₀, p₀, λ₀, λ₀)
        hdae8 = HDAE(eqs..., q₀, p₀, λ₀)

        @test hdae == hdae1
        @test hdae == hdae2
        @test hdae == hdae3
        @test hdae == hdae4
        @test hdae == hdae5
        @test hdae == hdae6
        @test hdae == hdae7
        @test hdae == hdae8

        @test hash(hdae) == hash(hdae1)
        @test hash(hdae) == hash(hdae2)
        @test hash(hdae) == hash(hdae3)
        @test hash(hdae) == hash(hdae4)
        @test hash(hdae) == hash(hdae5)
        @test hash(hdae) == hash(hdae6)
        @test hash(hdae) == hash(hdae7)
        @test hash(hdae) == hash(hdae8)
    end


    hdae_eqs  = (pdae_v, pdae_f, pdae_u, pdae_g, pdae_ϕ, pdae_u, pdae_g, pdae_ψ)
    hdae_args = (invariants=(h=pdae_h,), parameters=(a=1,), periodicity=π*ones(1))

    hdae  = HDAE(hdae_eqs..., pdae_v, pdae_f, symplectic_matrix, t₀, [q₀], [p₀], [λ₀], [λ₀], pdae_h, hdae_args.invariants, hdae_args.parameters, hdae_args.periodicity)

    @test axes(hdae) == axes(q₀)
    @test ndims(hdae) == 1
    @test nsamples(hdae) == 1
    @test nconstraints(hdae) == 1

    @test periodicity(hdae) == hdae_args.periodicity
    @test initial_conditions(hdae) == (t₀, [q₀], [p₀], [λ₀], [λ₀])

    @test hassecondary(hdae) == true
    @test hasinvariants(hdae) == true
    @test hasparameters(hdae) == true
    @test hasperiodicity(hdae) == true

    functions = get_functions(hdae)
    @test functions.v != pdae_v == hdae.v
    @test functions.f != pdae_f == hdae.f
    @test functions.u != pdae_u == hdae.u
    @test functions.g != pdae_g == hdae.g
    @test functions.ϕ != pdae_ϕ == hdae.ϕ
    @test functions.ū != pdae_u == hdae.ū
    @test functions.ḡ != pdae_g == hdae.ḡ
    @test functions.ψ != pdae_ψ == hdae.ψ
    @test functions.v̄ != pdae_v == hdae.v̄
    @test functions.f̄ != pdae_f == hdae.f̄
    @test functions.h != pdae_h == hdae.hamiltonian

    @test hdae == similar(hdae, t₀, q₀, p₀, λ₀)
    @test hdae == similar(hdae, t₀, q₀, p₀)
    @test hdae == similar(hdae, q₀, p₀)

    hdae1 = HDAE(hdae_eqs..., pdae_h, t₀, [q₀], [p₀], [λ₀], [λ₀]; hdae_args...)
    hdae2 = HDAE(hdae_eqs..., pdae_h, t₀, [q₀], [p₀], [λ₀]; hdae_args...)
    hdae3 = HDAE(hdae_eqs..., pdae_h, [q₀], [p₀], [λ₀], [λ₀]; hdae_args...)
    hdae4 = HDAE(hdae_eqs..., pdae_h, [q₀], [p₀], [λ₀]; hdae_args...)
    hdae5 = HDAE(hdae_eqs..., pdae_h, t₀, q₀, p₀, λ₀, λ₀; hdae_args...)
    hdae6 = HDAE(hdae_eqs..., pdae_h, t₀, q₀, p₀, λ₀; hdae_args...)
    hdae7 = HDAE(hdae_eqs..., pdae_h, q₀, p₀, λ₀, λ₀; hdae_args...)
    hdae8 = HDAE(hdae_eqs..., pdae_h, q₀, p₀, λ₀; hdae_args...)

    @test hdae == hdae1
    @test hdae == hdae2
    @test hdae == hdae3
    @test hdae == hdae4
    @test hdae == hdae5
    @test hdae == hdae6
    @test hdae == hdae7
    @test hdae == hdae8

    @test hash(hdae) == hash(hdae1)
    @test hash(hdae) == hash(hdae2)
    @test hash(hdae) == hash(hdae3)
    @test hash(hdae) == hash(hdae4)
    @test hash(hdae) == hash(hdae5)
    @test hash(hdae) == hash(hdae6)
    @test hash(hdae) == hash(hdae7)
    @test hash(hdae) == hash(hdae8)

end


@testset "$(rpad("Variational Differential Algebraic Equations (LDAE)",80))" begin

    ldae_eqs  = (iode_ϑ, iode_f, iode_u, iode_g, pdae_ϕ, nothing, nothing, nothing, lode_ω)
    ldae_eqs1 = (iode_ϑ, iode_f, iode_u, iode_g, pdae_ϕ, nothing, nothing, nothing, lode_l, lode_ω)
    ldae_eqs2 = (iode_ϑ, iode_f, iode_u, iode_g, pdae_ϕ, lode_l, lode_ω)

    ldae = LDAE(ldae_eqs..., iode_v, iode_f, t₀, [q₀], [p₀], [λ₀], [λ₀], lode_l, NullInvariants(), NullParameters(), nothing)

    @test ndims(ldae) == 1
    @test nsamples(ldae) == 1
    @test nconstraints(ldae) == 1

    @test periodicity(ldae) == zero(q₀)
    @test initial_conditions(ldae) == (t₀, [q₀], [p₀], [λ₀], [λ₀])

    @test hassecondary(ldae) == false
    @test hasinvariants(ldae) == false
    @test hasparameters(ldae) == false
    @test hasperiodicity(ldae) == false

    functions = get_functions(ldae)
    @test functions.ϑ == iode_ϑ == ldae.ϑ
    @test functions.f == iode_f == ldae.f
    @test functions.u == iode_u == ldae.u
    @test functions.g == iode_g == ldae.g
    @test functions.ϕ == pdae_ϕ == ldae.ϕ
    @test functions.v̄ == iode_v == ldae.v̄
    @test functions.f̄ == iode_f == ldae.f̄
    @test functions.ω == lode_ω == ldae.ω
    @test functions.l == lode_l == ldae.lagrangian

    @test ldae == similar(ldae, t₀, q₀, p₀, λ₀, λ₀)
    @test ldae == similar(ldae, t₀, q₀, p₀, λ₀)
    @test ldae == similar(ldae, t₀, q₀, p₀)
    @test ldae == similar(ldae, q₀, p₀)
    
    for eqs in (ldae_eqs1, ldae_eqs2)
        ldae1 = LDAE(eqs..., t₀, [q₀], [p₀], [λ₀], [λ₀]; v̄=iode_v)
        ldae2 = LDAE(eqs..., t₀, [q₀], [p₀], [λ₀]; v̄=iode_v)
        ldae3 = LDAE(eqs..., [q₀], [p₀], [λ₀], [λ₀]; v̄=iode_v)
        ldae4 = LDAE(eqs..., [q₀], [p₀], [λ₀]; v̄=iode_v)
        ldae5 = LDAE(eqs..., t₀, q₀, p₀, λ₀, λ₀; v̄=iode_v)
        ldae6 = LDAE(eqs..., t₀, q₀, p₀, λ₀; v̄=iode_v)
        ldae7 = LDAE(eqs..., q₀, p₀, λ₀, λ₀; v̄=iode_v)
        ldae8 = LDAE(eqs..., q₀, p₀, λ₀; v̄=iode_v)

        @test ldae == ldae1
        @test ldae == ldae2
        @test ldae == ldae3
        @test ldae == ldae4
        @test ldae == ldae5
        @test ldae == ldae6
        @test ldae == ldae7
        @test ldae == ldae8

        @test hash(ldae) == hash(ldae1)
        @test hash(ldae) == hash(ldae2)
        @test hash(ldae) == hash(ldae3)
        @test hash(ldae) == hash(ldae4)
        @test hash(ldae) == hash(ldae5)
        @test hash(ldae) == hash(ldae6)
        @test hash(ldae) == hash(ldae7)
        @test hash(ldae) == hash(ldae8)
    end


    ldae_eqs  = (iode_ϑ, iode_f, iode_u, iode_g, pdae_ϕ, pdae_u, pdae_g, pdae_ψ)
    ldae_args = (invariants=(h=iode_h,), parameters=(a=1,), periodicity=π*ones(1))

    ldae = LDAE(ldae_eqs..., lode_ω, iode_v, iode_f, t₀, [q₀], [p₀], [λ₀], [λ₀], lode_l, ldae_args.invariants, ldae_args.parameters, ldae_args.periodicity)

    @test ndims(ldae) == 1
    @test nsamples(ldae) == 1
    @test nconstraints(ldae) == 1

    @test periodicity(ldae) == ldae_args.periodicity
    @test initial_conditions(ldae) == (t₀, [q₀], [p₀], [λ₀], [λ₀])

    @test hassecondary(ldae) == true
    @test hasinvariants(ldae) == true
    @test hasparameters(ldae) == true
    @test hasperiodicity(ldae) == true

    functions = get_functions(ldae)
    @test functions.ϑ != iode_ϑ == ldae.ϑ
    @test functions.f != iode_f == ldae.f
    @test functions.u != iode_u == ldae.u
    @test functions.g != iode_g == ldae.g
    @test functions.ϕ != pdae_ϕ == ldae.ϕ
    @test functions.v̄ != iode_v == ldae.v̄
    @test functions.f̄ != iode_f == ldae.f̄
    @test functions.ω != lode_ω == ldae.ω
    @test functions.l != lode_l == ldae.lagrangian

    @test ldae == similar(ldae, t₀, q₀, p₀, λ₀, λ₀)
    @test ldae == similar(ldae, t₀, q₀, p₀, λ₀)
    @test ldae == similar(ldae, t₀, q₀, p₀)
    @test ldae == similar(ldae, q₀, p₀)

    ldae1 = LDAE(ldae_eqs..., lode_l, lode_ω, t₀, [q₀], [p₀], [λ₀], [λ₀]; v̄=iode_v, ldae_args...)
    ldae2 = LDAE(ldae_eqs..., lode_l, lode_ω, t₀, [q₀], [p₀], [λ₀]; v̄=iode_v, ldae_args...)
    ldae3 = LDAE(ldae_eqs..., lode_l, lode_ω, [q₀], [p₀], [λ₀], [λ₀]; v̄=iode_v, ldae_args...)
    ldae4 = LDAE(ldae_eqs..., lode_l, lode_ω, [q₀], [p₀], [λ₀]; v̄=iode_v, ldae_args...)
    ldae5 = LDAE(ldae_eqs..., lode_l, lode_ω, t₀, q₀, p₀, λ₀, λ₀; v̄=iode_v, ldae_args...)
    ldae6 = LDAE(ldae_eqs..., lode_l, lode_ω, t₀, q₀, p₀, λ₀; v̄=iode_v, ldae_args...)
    ldae7 = LDAE(ldae_eqs..., lode_l, lode_ω, q₀, p₀, λ₀, λ₀; v̄=iode_v, ldae_args...)
    ldae8 = LDAE(ldae_eqs..., lode_l, lode_ω, q₀, p₀, λ₀; v̄=iode_v, ldae_args...)

    @test ldae == ldae1
    @test ldae == ldae2
    @test ldae == ldae3
    @test ldae == ldae4
    @test ldae == ldae5
    @test ldae == ldae6
    @test ldae == ldae7
    @test ldae == ldae8

    @test hash(ldae) == hash(ldae1)
    @test hash(ldae) == hash(ldae2)
    @test hash(ldae) == hash(ldae3)
    @test hash(ldae) == hash(ldae4)
    @test hash(ldae) == hash(ldae5)
    @test hash(ldae) == hash(ldae6)
    @test hash(ldae) == hash(ldae7)
    @test hash(ldae) == hash(ldae8)
    
end


@testset "$(rpad("Split Partitioned Differential Algebraic Equations (SPDAE)",80))" begin

    spdae_v = (pdae_v, pdae_u, pdae_u)
    spdae_f = (pdae_f, pdae_g, pdae_g)
    spdae_eqs = (spdae_v, spdae_f, pdae_ϕ, pdae_ψ)

    spdae = SPDAE(spdae_eqs..., t₀, [q₀], [p₀], [λ₀], [λ₀], NullInvariants(), NullParameters(), nothing)

    @test ndims(spdae) == 1
    @test nsamples(spdae) == 1
    @test nconstraints(spdae) == 1

    @test periodicity(spdae) == zero(q₀)
    @test initial_conditions(spdae) == (t₀, [q₀], [p₀], [λ₀], [λ₀])

    @test hasinvariants(spdae) == false
    @test hasparameters(spdae) == false
    @test hasperiodicity(spdae) == false

    functions = get_functions(spdae)
    @test functions.v == spdae_v == spdae.v
    @test functions.f == spdae_f == spdae.f
    @test functions.ϕ ==  pdae_ϕ == spdae.ϕ
    @test functions.ψ ==  pdae_ψ == spdae.ψ

    @test spdae == similar(spdae, t₀, q₀, p₀, λ₀)
    @test spdae == similar(spdae, t₀, q₀, p₀)
    @test spdae == similar(spdae, q₀, p₀)

    spdae1 = SPDAE(spdae_eqs..., t₀, [q₀], [p₀], [λ₀], [λ₀])
    spdae2 = SPDAE(spdae_eqs..., t₀, [q₀], [p₀], [λ₀])
    spdae3 = SPDAE(spdae_eqs..., [q₀], [p₀], [λ₀], [λ₀])
    spdae4 = SPDAE(spdae_eqs..., [q₀], [p₀], [λ₀])
    spdae5 = SPDAE(spdae_eqs..., t₀, q₀, p₀, λ₀, λ₀)
    spdae6 = SPDAE(spdae_eqs..., t₀, q₀, p₀, λ₀)
    spdae7 = SPDAE(spdae_eqs..., q₀, p₀, λ₀, λ₀)
    spdae8 = SPDAE(spdae_eqs..., q₀, p₀, λ₀)

    @test spdae == spdae1
    @test spdae == spdae2
    @test spdae == spdae3
    @test spdae == spdae4
    @test spdae == spdae5
    @test spdae == spdae6
    @test spdae == spdae7
    @test spdae == spdae8

    @test hash(spdae) == hash(spdae1)
    @test hash(spdae) == hash(spdae2)
    @test hash(spdae) == hash(spdae3)
    @test hash(spdae) == hash(spdae4)
    @test hash(spdae) == hash(spdae5)
    @test hash(spdae) == hash(spdae6)
    @test hash(spdae) == hash(spdae7)
    @test hash(spdae) == hash(spdae8)

end
