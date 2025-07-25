
using GeometricEquations
using Test

include("functions.jl")
include("initial_conditions.jl")


@testset "$(rpad("Ordinary Differential Equations (ODE)",80))" begin

    ode  = ODE(ode_eqs..., ode_igs..., NullInvariants(), NullParameters(), NullPeriodicity())
    ode1 = ODE(ode_eqs...)
    ode2 = ODE(ode_eqs...; invariants=NullInvariants())
    ode3 = ODE(ode_eqs...; parameters=NullParameters())
    ode4 = ODE(ode_eqs...; periodicity=NullPeriodicity())

    @test ode == ode1
    @test ode == ode2
    @test ode == ode3
    @test ode == ode4

    @test hash(ode) == hash(ode1)
    @test hash(ode) == hash(ode2)
    @test hash(ode) == hash(ode3)
    @test hash(ode) == hash(ode4)

    @test functions(ode) == NamedTuple{(:v,)}(ode_eqs)
    @test solutions(ode) == NamedTuple()
    @test initialguess(ode) == NamedTuple{(:v,)}(ode_igs)

    @test parameters(ode) == NullParameters()
    @test invariants(ode) == NullInvariants()
    @test periodicity(ode) == NullPeriodicity()

    @test hasvectorfield(ode) == true
    @test hassolution(ode) == false
    @test hasinitialguess(ode) == true
    @test hasprimary(ode) == false
    @test hassecondary(ode) == false

    @test hasinvariants(ode) == false
    @test hasparameters(ode) == false
    @test hasperiodicity(ode) == false

    @test hashamiltonian(ode) == false
    @test haslagrangian(ode) == false

    funcs = functions(ode)

    @test_nowarn funcs.v(zero(x₀), t₀, x₀, NullParameters())

    funcs = functions(ode, NullParameters())

    @test_nowarn funcs.v(zero(x₀), t₀, x₀)


    # Test for periodicity
    odep = ODE(ode_eqs...; periodicity=([-π,0],[+π,2π]))

    @test periodicity(odep) == (Float64[-π,0],Float64[+π,2π])
    @test getperiodicity(odep) == BitArray([true,true])
    @test hasperiodicity(odep) == true

    odep = ODE(ode_eqs...; periodicity=([-Inf,0],[+Inf,2π]))

    @test periodicity(odep) == (Float64[-Inf,0],Float64[+Inf,2π])
    @test getperiodicity(odep) == BitArray([false,true])
    @test hasperiodicity(odep) == true

    odep = ODE(ode_eqs...; periodicity=([-Inf,-Inf],[+Inf,+Inf]))

    @test periodicity(odep) == NullPeriodicity()
    @test ismissing(getperiodicity(odep))
    @test hasperiodicity(odep) == false

end


@testset "$(rpad("Split Ordinary Differential Equations (SODE)",80))" begin

    @test_throws AssertionError SODE(nothing, nothing)
    @test_throws AssertionError SODE(sode_eqs, (sode_sols..., sode_sols...))

    @test_nowarn SODE(sode_eqs, nothing)
    @test_nowarn SODE(nothing, sode_sols)
    @test_nowarn SODE(sode_eqs, sode_sols)


    sode  = SODE(sode_eqs, nothing, sode_igs..., NullInvariants(), NullParameters(), NullPeriodicity())
    sode0 = SODE(sode_eqs)
    sode1 = SODE(sode_eqs; v̄ = ode_v)
    sode2 = SODE(sode_eqs; v̄ = ode_v, invariants=NullInvariants())
    sode3 = SODE(sode_eqs; v̄ = ode_v, parameters=NullParameters())
    sode4 = SODE(sode_eqs; v̄ = ode_v, periodicity=NullPeriodicity())

    @test sode != sode0
    @test sode == sode1
    @test sode == sode2
    @test sode == sode3
    @test sode == sode4

    @test hash(sode) == hash(sode1)
    @test hash(sode) == hash(sode2)
    @test hash(sode) == hash(sode3)
    @test hash(sode) == hash(sode4)

    @test nsteps(sode) == 2
    @test functions(sode) == (v = sode_eqs,)
    @test solutions(sode) == NamedTuple()
    @test initialguess(sode) == (v = ode_v,)

    @test parameters(sode) == NullParameters()
    @test invariants(sode) == NullInvariants()
    @test periodicity(sode) == NullPeriodicity()

    @test hasvectorfield(sode) == true
    @test hassolution(sode) == false
    @test hasinitialguess(sode) == true
    @test hasprimary(sode) == false
    @test hassecondary(sode) == false

    @test hasinvariants(sode) == false
    @test hasparameters(sode) == false
    @test hasperiodicity(sode) == false

    @test hashamiltonian(sode) == false
    @test haslagrangian(sode) == false


    sode  = SODE(sode_eqs, sode_sols, ode_v, NullInvariants(), NullParameters(), NullPeriodicity())
    sode0 = SODE(sode_eqs, sode_sols)
    sode1 = SODE(sode_eqs, sode_sols; v̄ = ode_v)
    sode2 = SODE(sode_eqs, sode_sols; v̄ = ode_v, invariants=NullInvariants())
    sode3 = SODE(sode_eqs, sode_sols; v̄ = ode_v, parameters=NullParameters())
    sode4 = SODE(sode_eqs, sode_sols; v̄ = ode_v, periodicity=NullPeriodicity())

    @test sode == sode1
    @test sode == sode2
    @test sode == sode3
    @test sode == sode4

    @test hash(sode) == hash(sode1)
    @test hash(sode) == hash(sode2)
    @test hash(sode) == hash(sode3)
    @test hash(sode) == hash(sode4)

    @test nsteps(sode) == 2
    @test functions(sode) == (v = sode_eqs,)
    @test solutions(sode) == (q = sode_sols,)
    @test initialguess(sode) == (v = ode_v,)

    @test parameters(sode) == NullParameters()
    @test invariants(sode) == NullInvariants()
    @test periodicity(sode) == NullPeriodicity()

    @test hasvectorfield(sode) == true
    @test hassolution(sode) == true
    @test hasinitialguess(sode) == true
    @test hasprimary(sode) == false
    @test hassecondary(sode) == false

    @test hasinvariants(sode) == false
    @test hasparameters(sode) == false
    @test hasperiodicity(sode) == false

    @test hashamiltonian(sode) == false
    @test haslagrangian(sode) == false


    sode  = SODE(nothing, sode_sols, ode_v, NullInvariants(), NullParameters(), NullPeriodicity())
    sode0 = SODE(nothing, sode_sols)
    sode1 = SODE(nothing, sode_sols, v̄ = ode_v)
    sode2 = SODE(nothing, sode_sols; v̄ = ode_v, invariants=NullInvariants())
    sode3 = SODE(nothing, sode_sols; v̄ = ode_v, parameters=NullParameters())
    sode4 = SODE(nothing, sode_sols; v̄ = ode_v, periodicity=NullPeriodicity())

    @test sode == sode1
    @test sode == sode2
    @test sode == sode3
    @test sode == sode4

    @test hash(sode) == hash(sode1)
    @test hash(sode) == hash(sode2)
    @test hash(sode) == hash(sode3)
    @test hash(sode) == hash(sode4)

    @test nsteps(sode) == 2
    @test functions(sode) == NamedTuple()
    @test solutions(sode) == (q = sode_sols,)
    @test initialguess(sode) == (v = ode_v,)

    @test parameters(sode) == NullParameters()
    @test invariants(sode) == NullInvariants()
    @test periodicity(sode) == NullPeriodicity()

    @test hasvectorfield(sode) == false
    @test hassolution(sode) == true
    @test hasinitialguess(sode) == true
    @test hasprimary(sode) == false
    @test hassecondary(sode) == false

    @test hasinvariants(sode) == false
    @test hasparameters(sode) == false
    @test hasperiodicity(sode) == false

    @test hashamiltonian(sode) == false
    @test haslagrangian(sode) == false

    # Test for periodicity
    sodep = SODE(sode_eqs; periodicity=([-π,0],[+π,2π]))

    @test periodicity(sodep) == (Float64[-π,0],Float64[+π,2π])
    @test getperiodicity(sodep) == BitArray([true,true])
    @test hasperiodicity(sodep) == true

    sodep = SODE(sode_eqs; periodicity=([-Inf,0],[+Inf,2π]))

    @test periodicity(sodep) == (Float64[-Inf,0],Float64[+Inf,2π])
    @test getperiodicity(sodep) == BitArray([false,true])
    @test hasperiodicity(sodep) == true

    sodep = SODE(sode_eqs; periodicity=([-Inf,-Inf],[+Inf,+Inf]))

    @test periodicity(sodep) == NullPeriodicity()
    @test ismissing(getperiodicity(sodep))
    @test hasperiodicity(sodep) == false

end


@testset "$(rpad("Partitioned Ordinary Differential Equations (PODE)",80))" begin

    pode  = PODE(pode_eqs..., pode_igs..., NullInvariants(), NullParameters(), NullPeriodicity())
    pode1 = PODE(pode_eqs...)
    pode2 = PODE(pode_eqs...; invariants=NullInvariants())
    pode3 = PODE(pode_eqs...; parameters=NullParameters())
    pode4 = PODE(pode_eqs...; periodicity=NullPeriodicity())

    @test pode == pode1
    @test pode == pode2
    @test pode == pode3
    @test pode == pode4

    @test hash(pode) == hash(pode1)
    @test hash(pode) == hash(pode2)
    @test hash(pode) == hash(pode3)
    @test hash(pode) == hash(pode4)

    @test functions(pode) == NamedTuple{(:v,:f)}(pode_eqs)
    @test solutions(pode) == NamedTuple()
    @test initialguess(pode) == NamedTuple{(:v,:f)}(pode_igs)

    @test parameters(pode) == NullParameters()
    @test invariants(pode) == NullInvariants()
    @test periodicity(pode) == NullPeriodicity()

    @test hasvectorfield(pode) == true
    @test hassolution(pode) == false
    @test hasinitialguess(pode) == true
    @test hasprimary(pode) == false
    @test hassecondary(pode) == false

    @test hasinvariants(pode) == false
    @test hasparameters(pode) == false
    @test hasperiodicity(pode) == false

    @test hashamiltonian(pode) == false
    @test haslagrangian(pode) == false

    funcs = functions(pode)

    @test_nowarn funcs.v(zero(q₀), t₀, q₀, p₀, NullParameters())
    @test_nowarn funcs.f(zero(p₀), t₀, q₀, p₀, NullParameters())

    funcs = functions(pode, NullParameters())

    @test_nowarn funcs.v(zero(q₀), t₀, q₀, p₀)
    @test_nowarn funcs.f(zero(p₀), t₀, q₀, p₀)

    # funcs = functions(pode)
    # @test funcs.v == pode_v == pode.v
    # @test funcs.f == pode_f == pode.f

    # rode = ODE(ode_v, t₀, [x₀])
    # code = convert(ODE, pode)
    # v₁ = zero(x₀)
    # v₂ = zero(x₀)
    # rode.v(rode.t₀, rode.q₀[begin], v₁)
    # code.v(code.t₀, code.q₀[begin], v₂)
    # @test v₁ == v₂

    # rode = SODE((v_sode_1, v_sode_2), t₀, [x₀])
    # code = convert(SODE, pode)
    # v₁ = zero(x₀)
    # v₂ = zero(x₀)
    # rode.v[1](rode.t₀, rode.q₀[begin], v₁)
    # code.v[1](code.t₀, code.q₀[begin], v₂)
    # @test v₁ == v₂
    # rode.v[2](rode.t₀, rode.q₀[begin], v₁)
    # code.v[2](code.t₀, code.q₀[begin], v₂)
    # @test v₁ == v₂


    # Test for periodicity
    podep = PODE(pode_eqs...; periodicity=([0.0,],[2π]))

    @test periodicity(podep) == (Float64[0],Float64[2π])
    @test getperiodicity(podep) == BitArray([true])
    @test hasperiodicity(podep) == true

    podep = PODE(pode_eqs...; periodicity=([-Inf],[+Inf]))

    @test periodicity(podep) == NullPeriodicity()
    @test ismissing(getperiodicity(podep))
    @test hasperiodicity(podep) == false

end


@testset "$(rpad("Implicit Ordinary Differential Equations (IODE)",80))" begin

    iode  = IODE(iode_eqs..., iode_igs..., NullInvariants(), NullParameters(), NullPeriodicity())

    iode1 = IODE(iode_eqs...; v̄ = iode_v)
    iode2 = IODE(iode_eqs...; v̄ = iode_v, invariants=NullInvariants())
    iode3 = IODE(iode_eqs...; v̄ = iode_v, parameters=NullParameters())
    iode4 = IODE(iode_eqs...; v̄ = iode_v, periodicity=NullPeriodicity())
    iode5 = IODE(iode_eqs...; v̄ = iode_v, f̄ = iode_f)
    iode6 = IODE(iode_eqs...; v̄ = iode_v, f̄ = iode_f, invariants=NullInvariants())
    iode7 = IODE(iode_eqs...; v̄ = iode_v, f̄ = iode_f, parameters=NullParameters())
    iode8 = IODE(iode_eqs...; v̄ = iode_v, f̄ = iode_f, periodicity=NullPeriodicity())

    @test iode == iode1
    @test iode == iode2
    @test iode == iode3
    @test iode == iode4
    @test iode == iode5
    @test iode == iode6
    @test iode == iode7
    @test iode == iode8

    @test hash(iode) == hash(iode1)
    @test hash(iode) == hash(iode2)
    @test hash(iode) == hash(iode3)
    @test hash(iode) == hash(iode4)
    @test hash(iode) == hash(iode5)
    @test hash(iode) == hash(iode6)
    @test hash(iode) == hash(iode7)
    @test hash(iode) == hash(iode8)

    @test functions(iode) == NamedTuple{(:ϑ,:f,:g)}(iode_eqs)
    @test solutions(iode) == NamedTuple()
    @test initialguess(iode) == NamedTuple{(:v,:f)}(iode_igs)

    @test parameters(iode) == NullParameters()
    @test invariants(iode) == NullInvariants()
    @test periodicity(iode) == NullPeriodicity()

    @test hasvectorfield(iode) == true
    @test hassolution(iode) == false
    @test hasinitialguess(iode) == true
    @test hasprimary(iode) == false
    @test hassecondary(iode) == false

    @test hasinvariants(iode) == false
    @test hasparameters(iode) == false
    @test hasperiodicity(iode) == false

    @test hashamiltonian(iode) == false
    @test haslagrangian(iode) == false

    funcs = functions(iode)

    @test_nowarn funcs.ϑ(zero(p₀), t₀, q₀, p₀, NullParameters())
    @test_nowarn funcs.f(zero(p₀), t₀, q₀, p₀, NullParameters())
    @test_nowarn funcs.g(zero(p₀), t₀, q₀, p₀, λ₀, NullParameters())

    funcs = functions(iode, NullParameters())

    @test_nowarn funcs.ϑ(zero(p₀), t₀, q₀, p₀)
    @test_nowarn funcs.f(zero(p₀), t₀, q₀, p₀)
    @test_nowarn funcs.g(zero(p₀), t₀, q₀, p₀, λ₀)

    # @test periodicity(iode) == zero(q₀)
    # @test initial_conditions(iode) == (t₀, [q₀], [p₀], [λ₀])

    # funcs = functions(iode)
    # @test funcs.ϑ == iode_ϑ == iode.ϑ
    # @test funcs.f == iode_f == iode.f
    # @test funcs.g == iode_g == iode.g
    # @test funcs.v̄ == iode_v == iode.v̄
    # @test funcs.f̄ == iode_f == iode.f̄

    # @test iode == similar(iode, t₀, q₀, p₀, λ₀)
    # @test iode == similar(iode, t₀, q₀, p₀)
    # @test iode == similar(iode, q₀, p₀)

    # Test for periodicity
    iodep = IODE(iode_eqs...; periodicity=([0.0,],[2π]))

    @test periodicity(iodep) == (Float64[0],Float64[2π])
    @test getperiodicity(iodep) == BitArray([true])
    @test hasperiodicity(iodep) == true

    iodep = IODE(iode_eqs...; periodicity=([-Inf],[+Inf]))

    @test periodicity(iodep) == NullPeriodicity()
    @test ismissing(getperiodicity(iodep))
    @test hasperiodicity(iodep) == false

end


@testset "$(rpad("Hamiltonian Ordinary Differential Equations (HODE)",80))" begin

    hode  = HODE(pode_eqs..., hode_igs..., hode_h, NullInvariants(), NullParameters(), NullPeriodicity())

    hode1 = HODE(hode_eqs...)
    hode2 = HODE(hode_eqs...; invariants=NullInvariants())
    hode3 = HODE(hode_eqs...; parameters=NullParameters())
    hode4 = HODE(hode_eqs...; periodicity=NullPeriodicity())

    @test hode == hode1
    @test hode == hode2
    @test hode == hode3
    @test hode == hode4

    @test hash(hode) == hash(hode1)
    @test hash(hode) == hash(hode2)
    @test hash(hode) == hash(hode3)
    @test hash(hode) == hash(hode4)

    @test functions(hode) == NamedTuple{(:v,:f,:h)}(hode_eqs)
    @test solutions(hode) == NamedTuple()
    @test initialguess(hode) == NamedTuple{(:v,:f)}(hode_igs)

    @test parameters(hode) == NullParameters()
    @test invariants(hode) == NullInvariants()
    @test periodicity(hode) == NullPeriodicity()

    @test hasvectorfield(hode) == true
    @test hassolution(hode) == false
    @test hasinitialguess(hode) == true
    @test hasprimary(hode) == false
    @test hassecondary(hode) == false

    @test hasinvariants(hode) == false
    @test hasparameters(hode) == false
    @test hasperiodicity(hode) == false

    @test hashamiltonian(hode) == true
    @test haslagrangian(hode) == false

    funcs = functions(hode)

    @test_nowarn funcs.v(zero(q₀), t₀, q₀, p₀, NullParameters())
    @test_nowarn funcs.f(zero(p₀), t₀, q₀, p₀, NullParameters())
    @test_nowarn funcs.h(t₀, q₀, p₀, NullParameters())
    # @test_nowarn funcs.poisson(t₀, q₀, p₀, zeros(2,2), NullParameters())

    funcs = functions(hode, NullParameters())

    @test_nowarn funcs.v(zero(q₀), t₀, q₀, p₀)
    @test_nowarn funcs.f(zero(p₀), t₀, q₀, p₀)
    @test_nowarn funcs.h(t₀, q₀, p₀)
    # @test_nowarn funcs.poisson(t₀, q₀, p₀, zeros(2,2))

    # funcs = functions(hode)
    # @test funcs.v == pode_v == hode.v
    # @test funcs.f == pode_f == hode.f
    # @test funcs.h == pode_h == hode.hamiltonian

    # @test hode == similar(hode, t₀, q₀, p₀)
    # @test hode == similar(hode, q₀, p₀)

    # rode = ODE(ode_v, t₀, [x₀])
    # code = convert(ODE, hode)
    # v₁ = zero(x₀)
    # v₂ = zero(x₀)
    # rode.v(rode.t₀, rode.q₀[begin], v₁)
    # code.v(code.t₀, code.q₀[begin], v₂)
    # @test v₁ == v₂

    # pode = convert(PODE, hode)
    # v₁ = zero(q₀)
    # v₂ = zero(q₀)
    # f₁ = zero(p₀)
    # f₂ = zero(p₀)
    # hode.v(hode.t₀, hode.q₀[begin], hode.p₀[begin], v₁)
    # pode.v(pode.t₀, pode.q₀[begin], pode.p₀[begin], v₂)
    # hode.f(hode.t₀, hode.q₀[begin], hode.p₀[begin], f₁)
    # pode.f(pode.t₀, pode.q₀[begin], pode.p₀[begin], f₂)
    # @test v₁ == v₂
    # @test f₁ == f₂

    # rode = SODE((v_sode_1, v_sode_2), t₀, [x₀])
    # code = convert(SODE, hode)
    # v₁ = zero(x₀)
    # v₂ = zero(x₀)
    # rode.v[1](rode.t₀, rode.q₀[begin], v₁)
    # code.v[1](code.t₀, code.q₀[begin], v₂)
    # @test v₁ == v₂
    # rode.v[2](rode.t₀, rode.q₀[begin], v₁)
    # code.v[2](code.t₀, code.q₀[begin], v₂)
    # @test v₁ == v₂

    # Test for periodicity
    hodep = HODE(hode_eqs...; periodicity=([0.0,],[2π]))

    @test periodicity(hodep) == (Float64[0],Float64[2π])
    @test getperiodicity(hodep) == BitArray([true])
    @test hasperiodicity(hodep) == true

    hodep = HODE(hode_eqs...; periodicity=([-Inf],[+Inf]))

    @test periodicity(hodep) == NullPeriodicity()
    @test ismissing(getperiodicity(hodep))
    @test hasperiodicity(hodep) == false

end


@testset "$(rpad("Lagrangian Ordinary Differential Equations (LODE)",80))" begin

    lode = LODE(iode_eqs..., lode_ω, lode_igs..., lode_l, NullInvariants(), NullParameters(), NullPeriodicity())

    lode1 = LODE(lode_eqs...; v̄ = lode_v, f̄ = lode_f)
    lode2 = LODE(lode_eqs...; v̄ = lode_v, f̄ = lode_f, invariants=NullInvariants())
    lode3 = LODE(lode_eqs...; v̄ = lode_v, f̄ = lode_f, parameters=NullParameters())
    lode4 = LODE(lode_eqs...; v̄ = lode_v, f̄ = lode_f, periodicity=NullPeriodicity())

    @test lode == lode1
    @test lode == lode2
    @test lode == lode3
    @test lode == lode4

    @test hash(lode) == hash(lode1)
    @test hash(lode) == hash(lode2)
    @test hash(lode) == hash(lode3)
    @test hash(lode) == hash(lode4)

    @test functions(lode) == NamedTuple{(:ϑ,:f,:g,:ω,:l)}(lode_eqs)
    @test solutions(lode) == NamedTuple()
    @test initialguess(lode) == NamedTuple{(:v,:f)}(lode_igs)

    @test parameters(lode) == NullParameters()
    @test invariants(lode) == NullInvariants()
    @test periodicity(lode) == NullPeriodicity()

    @test hasvectorfield(lode) == true
    @test hassolution(lode) == false
    @test hasinitialguess(lode) == true
    @test hasprimary(lode) == false
    @test hassecondary(lode) == false

    @test hasinvariants(lode) == false
    @test hasparameters(lode) == false
    @test hasperiodicity(lode) == false

    @test hashamiltonian(lode) == false
    @test haslagrangian(lode) == true

    funcs = functions(lode)

    @test_nowarn funcs.ϑ(zero(p₀), t₀, q₀, p₀, NullParameters())
    @test_nowarn funcs.f(zero(p₀), t₀, q₀, p₀, NullParameters())
    @test_nowarn funcs.g(zero(p₀), t₀, q₀, p₀, λ₀, NullParameters())
    @test_nowarn funcs.ω(zeros(2,2), t₀, q₀, p₀, NullParameters())
    @test_nowarn funcs.l(t₀, q₀, p₀, NullParameters())

    funcs = functions(lode, NullParameters())

    @test_nowarn funcs.ϑ(zero(p₀), t₀, q₀, p₀)
    @test_nowarn funcs.f(zero(p₀), t₀, q₀, p₀)
    @test_nowarn funcs.g(zero(p₀), t₀, q₀, p₀, λ₀)
    @test_nowarn funcs.ω(zeros(2,2), t₀, q₀, p₀)
    @test_nowarn funcs.l(t₀, q₀, p₀)

    # funcs = functions(lode)
    # @test funcs.ϑ == iode_ϑ == lode.ϑ
    # @test funcs.f == iode_f == lode.f
    # @test funcs.g == iode_g == lode.g
    # @test funcs.ω == lode_ω == lode.ω
    # @test funcs.v̄ == iode_v == lode.v̄
    # @test funcs.f̄ == iode_f == lode.f̄
    # @test funcs.l == lode_l == lode.lagrangian

    # @test lode == similar(lode, t₀, q₀, p₀, λ₀)
    # @test lode == similar(lode, t₀, q₀, p₀)
    # @test lode == similar(lode, q₀, p₀)

    # iode = convert(IODE, lode)
    # v₁ = zero(q₀)
    # v₂ = zero(q₀)
    # ϑ₁ = zero(q₀)
    # ϑ₂ = zero(q₀)
    # f₁ = zero(p₀)
    # f₂ = zero(p₀)
    # g₁ = zero(p₀)
    # g₂ = zero(p₀)
    # lode.v̄(lode.t₀, lode.q₀[begin], lode.p₀[begin], v₁)
    # iode.v̄(iode.t₀, iode.q₀[begin], iode.p₀[begin], v₂)
    # lode.ϑ(lode.t₀, lode.q₀[begin], v₁, ϑ₁)
    # iode.ϑ(iode.t₀, iode.q₀[begin], v₂, ϑ₂)
    # lode.f(lode.t₀, lode.q₀[begin], v₁, f₁)
    # iode.f(iode.t₀, iode.q₀[begin], v₂, f₂)
    # lode.g(lode.t₀, lode.q₀[begin], v₁, g₁)
    # iode.g(iode.t₀, iode.q₀[begin], v₂, g₂)
    # @test v₁ == v₂
    # @test ϑ₁ == ϑ₂
    # @test f₁ == f₂
    # @test g₁ == g₂

    # Test for periodicity
    lodep = LODE(lode_eqs...; periodicity=([0.0,],[2π]))

    @test periodicity(lodep) == (Float64[0],Float64[2π])
    @test getperiodicity(lodep) == BitArray([true])
    @test hasperiodicity(lodep) == true

    lodep = LODE(lode_eqs...; periodicity=([-Inf],[+Inf]))

    @test periodicity(lodep) == NullPeriodicity()
    @test ismissing(getperiodicity(lodep))
    @test hasperiodicity(lodep) == false

end
