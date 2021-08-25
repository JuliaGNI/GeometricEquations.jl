
using GeometricEquations
using GeometricEquations: symplectic_matrix
using Test

include("functions.jl")
include("initial_conditions.jl")


@testset "$(rpad("Ordinary Differential Equations (ODE)",80))" begin

    ode  = ODE(ode_v, t₀, [x₀], NullInvariants(), NullParameters(), nothing)

    ode1 = ODE(ode_v, t₀, [x₀])
    ode2 = ODE(ode_v, [x₀])
    ode3 = ODE(ode_v, t₀, x₀)
    ode4 = ODE(ode_v, x₀)

    @test axes(ode) == axes(x₀)
    @test ndims(ode) == 2
    @test nsamples(ode) == 1

    @test periodicity(ode) == zero(x₀)
    @test initial_conditions(ode) == (t₀, [x₀])

    @test hasinvariants(ode) == false
    @test hasparameters(ode) == false
    @test hasperiodicity(ode) == false

    @test get_functions(ode) == NamedTuple{(:v,)}((ode_v,))

    @test ode == similar(ode, t₀, x₀)
    @test ode == similar(ode, x₀)

    @test ode == ode1
    @test ode == ode2
    @test ode == ode3
    @test ode == ode4

    @test hash(ode) == hash(ode1)
    @test hash(ode) == hash(ode2)
    @test hash(ode) == hash(ode3)
    @test hash(ode) == hash(ode4)

end


@testset "$(rpad("Split Ordinary Differential Equations (SODE)",80))" begin

    v_sode = (v_sode_1, v_sode_2)
    q_sode = (q_sode_1, q_sode_2)

    sode  = SODE(v_sode, nothing, t₀, [x₀], NullInvariants(), NullParameters(), nothing)
    sode1 = SODE(v_sode, t₀, [x₀])
    sode2 = SODE(v_sode, [x₀])
    sode3 = SODE(v_sode, t₀, x₀)
    sode4 = SODE(v_sode, x₀)

    @test ndims(sode) == 2
    @test nsamples(sode) == 1
    @test periodicity(sode) == zero(x₀)

    @test sode == similar(sode, t₀, x₀)
    @test sode == similar(sode, x₀)

    @test sode == sode1
    @test sode == sode2
    @test sode == sode3
    @test sode == sode4

    @test hash(sode) == hash(sode1)
    @test hash(sode) == hash(sode2)
    @test hash(sode) == hash(sode3)
    @test hash(sode) == hash(sode4)


    sode  = SODE(v_sode, q_sode, t₀, [x₀])
    sode1 = SODE(v_sode, q_sode, t₀, [x₀])
    sode2 = SODE(v_sode, q_sode, [x₀])
    sode3 = SODE(v_sode, q_sode, t₀, x₀)
    sode4 = SODE(v_sode, q_sode, x₀)

    @test sode == similar(sode, t₀, x₀)
    @test sode == similar(sode, x₀)

    @test sode == sode1
    @test sode == sode2
    @test sode == sode3
    @test sode == sode4

    @test hash(sode) == hash(sode1)
    @test hash(sode) == hash(sode2)
    @test hash(sode) == hash(sode3)
    @test hash(sode) == hash(sode4)

end


@testset "$(rpad("Partitioned Ordinary Differential Equations (PODE)",80))" begin

    pode_eqs = (pode_v, pode_f)

    pode  = PODE(pode_eqs..., t₀, [q₀], [p₀], NullInvariants(), NullParameters(), nothing)

    pode1 = PODE(pode_eqs..., t₀, [q₀], [p₀])
    pode2 = PODE(pode_eqs..., [q₀], [p₀])
    pode3 = PODE(pode_eqs..., t₀, q₀, p₀)
    pode4 = PODE(pode_eqs..., q₀, p₀)

    @test axes(pode) == axes(q₀)
    @test ndims(pode) == 1
    @test nsamples(pode) == 1

    @test periodicity(pode) == zero(q₀)
    @test initial_conditions(pode) == (t₀, [q₀], [p₀])

    @test hasinvariants(pode) == false
    @test hasparameters(pode) == false
    @test hasperiodicity(pode) == false

    functions = get_functions(pode)
    @test functions.v == pode_v == pode.v
    @test functions.f == pode_f == pode.f

    @test pode == similar(pode, t₀, q₀, p₀)
    @test pode == similar(pode, q₀, p₀)

    @test pode == pode1
    @test pode == pode2
    @test pode == pode3
    @test pode == pode4

    @test hash(pode) == hash(pode1)
    @test hash(pode) == hash(pode2)
    @test hash(pode) == hash(pode3)
    @test hash(pode) == hash(pode4)

    rode = ODE(ode_v, t₀, [x₀])
    code = convert(ODE, pode)
    v₁ = zero(x₀)
    v₂ = zero(x₀)
    rode.v(rode.t₀, rode.q₀[begin], v₁)
    code.v(code.t₀, code.q₀[begin], v₂)
    @test v₁ == v₂

    rode = SODE((v_sode_1, v_sode_2), t₀, [x₀])
    code = convert(SODE, pode)
    v₁ = zero(x₀)
    v₂ = zero(x₀)
    rode.v[1](rode.t₀, rode.q₀[begin], v₁)
    code.v[1](code.t₀, code.q₀[begin], v₂)
    @test v₁ == v₂
    rode.v[2](rode.t₀, rode.q₀[begin], v₁)
    code.v[2](code.t₀, code.q₀[begin], v₂)
    @test v₁ == v₂

end


@testset "$(rpad("Implicit Ordinary Differential Equations (IODE)",80))" begin

    iode_eqs = (iode_ϑ, iode_f, iode_g)

    iode  = IODE(iode_eqs..., iode_v, iode_f, t₀, [q₀], [p₀], [λ₀], NullInvariants(), NullParameters(), nothing)

    iode1 = IODE(iode_eqs..., t₀, [q₀], [p₀], [λ₀]; v̄=iode_v)
    iode2 = IODE(iode_eqs..., t₀, q₀, p₀, λ₀; v̄=iode_v)
    iode3 = IODE(iode_eqs..., t₀, q₀, p₀; v̄=iode_v)
    iode4 = IODE(iode_eqs..., q₀, p₀; v̄=iode_v)

    @test axes(iode) == axes(q₀)
    @test ndims(iode) == 1
    @test nsamples(iode) == 1

    @test periodicity(iode) == zero(q₀)
    @test initial_conditions(iode) == (t₀, [q₀], [p₀], [λ₀])

    @test hasinvariants(iode) == false
    @test hasparameters(iode) == false
    @test hasperiodicity(iode) == false

    functions = get_functions(iode)
    @test functions.ϑ == iode_ϑ == iode.ϑ
    @test functions.f == iode_f == iode.f
    @test functions.g == iode_g == iode.g
    @test functions.v̄ == iode_v == iode.v̄
    @test functions.f̄ == iode_f == iode.f̄

    @test iode == similar(iode, t₀, q₀, p₀, λ₀)
    @test iode == similar(iode, t₀, q₀, p₀)
    @test iode == similar(iode, q₀, p₀)
    
    @test iode == iode1
    @test iode == iode2
    @test iode == iode3
    @test iode == iode4

    @test hash(iode) == hash(iode1)
    @test hash(iode) == hash(iode2)
    @test hash(iode) == hash(iode3)
    @test hash(iode) == hash(iode4)

end


@testset "$(rpad("Hamiltonian Ordinary Differential Equations (HODE)",80))" begin

    hode_eqs = (pode_v, pode_f, pode_h)

    hode  = HODE(pode_v, pode_f, symplectic_matrix, t₀, [q₀], [p₀], pode_h, NullInvariants(), NullParameters(), nothing)

    hode1 = HODE(hode_eqs..., t₀, [q₀], [p₀])
    hode2 = HODE(hode_eqs..., [q₀], [p₀])
    hode3 = HODE(hode_eqs..., t₀, q₀, p₀)
    hode4 = HODE(hode_eqs..., q₀, p₀)

    @test axes(hode) == axes(q₀)
    @test ndims(hode) == 1
    @test nsamples(hode) == 1

    @test periodicity(hode) == zero(q₀)
    @test initial_conditions(hode) == (t₀, [q₀], [p₀])

    @test hasinvariants(hode) == false
    @test hasparameters(hode) == false
    @test hasperiodicity(hode) == false

    functions = get_functions(hode)
    @test functions.v == pode_v == hode.v
    @test functions.f == pode_f == hode.f
    @test functions.h == pode_h == hode.hamiltonian

    @test hode == similar(hode, t₀, q₀, p₀)
    @test hode == similar(hode, q₀, p₀)

    @test hode == hode1
    @test hode == hode2
    @test hode == hode3
    @test hode == hode4

    @test hash(hode) == hash(hode1)
    @test hash(hode) == hash(hode2)
    @test hash(hode) == hash(hode3)
    @test hash(hode) == hash(hode4)

    rode = ODE(ode_v, t₀, [x₀])
    code = convert(ODE, hode)
    v₁ = zero(x₀)
    v₂ = zero(x₀)
    rode.v(rode.t₀, rode.q₀[begin], v₁)
    code.v(code.t₀, code.q₀[begin], v₂)
    @test v₁ == v₂

    pode = convert(PODE, hode)
    v₁ = zero(q₀)
    v₂ = zero(q₀)
    f₁ = zero(p₀)
    f₂ = zero(p₀)
    hode.v(hode.t₀, hode.q₀[begin], hode.p₀[begin], v₁)
    pode.v(pode.t₀, pode.q₀[begin], pode.p₀[begin], v₂)
    hode.f(hode.t₀, hode.q₀[begin], hode.p₀[begin], f₁)
    pode.f(pode.t₀, pode.q₀[begin], pode.p₀[begin], f₂)
    @test v₁ == v₂
    @test f₁ == f₂

    rode = SODE((v_sode_1, v_sode_2), t₀, [x₀])
    code = convert(SODE, hode)
    v₁ = zero(x₀)
    v₂ = zero(x₀)
    rode.v[1](rode.t₀, rode.q₀[begin], v₁)
    code.v[1](code.t₀, code.q₀[begin], v₂)
    @test v₁ == v₂
    rode.v[2](rode.t₀, rode.q₀[begin], v₁)
    code.v[2](code.t₀, code.q₀[begin], v₂)
    @test v₁ == v₂

end


@testset "$(rpad("Variational Ordinary Differential Equations (LODE)",80))" begin

    lode_eqs = (iode_ϑ, iode_f, iode_g, lode_l, lode_ω)

    lode = LODE(iode_ϑ, iode_f, iode_g, lode_ω, iode_v, iode_f, t₀, [q₀], [p₀], [λ₀], lode_l, NullInvariants(), NullParameters(), nothing)

    lode1 = LODE(lode_eqs..., t₀, [q₀], [p₀], [λ₀]; v̄=iode_v)
    lode2 = LODE(lode_eqs..., t₀, q₀, p₀, λ₀; v̄=iode_v)
    lode3 = LODE(lode_eqs..., t₀, q₀, p₀; v̄=iode_v)
    lode4 = LODE(lode_eqs..., q₀, p₀; v̄=iode_v)

    @test axes(lode) == axes(q₀)
    @test ndims(lode) == 1
    @test nsamples(lode) == 1

    @test periodicity(lode) == zero(q₀)
    @test initial_conditions(lode) == (t₀, [q₀], [p₀], [λ₀])

    @test hasinvariants(lode) == false
    @test hasparameters(lode) == false
    @test hasperiodicity(lode) == false

    functions = get_functions(lode)
    @test functions.ϑ == iode_ϑ == lode.ϑ
    @test functions.f == iode_f == lode.f
    @test functions.g == iode_g == lode.g
    @test functions.ω == lode_ω == lode.ω
    @test functions.v̄ == iode_v == lode.v̄
    @test functions.f̄ == iode_f == lode.f̄
    @test functions.l == lode_l == lode.lagrangian

    @test lode == similar(lode, t₀, q₀, p₀, λ₀)
    @test lode == similar(lode, t₀, q₀, p₀)
    @test lode == similar(lode, q₀, p₀)

    @test lode == lode1
    @test lode == lode2
    @test lode == lode3
    @test lode == lode4

    @test hash(lode) == hash(lode1)
    @test hash(lode) == hash(lode2)
    @test hash(lode) == hash(lode3)
    @test hash(lode) == hash(lode4)

    iode = convert(IODE, lode)
    v₁ = zero(q₀)
    v₂ = zero(q₀)
    ϑ₁ = zero(q₀)
    ϑ₂ = zero(q₀)
    f₁ = zero(p₀)
    f₂ = zero(p₀)
    g₁ = zero(p₀)
    g₂ = zero(p₀)
    lode.v̄(lode.t₀, lode.q₀[begin], lode.p₀[begin], v₁)
    iode.v̄(iode.t₀, iode.q₀[begin], iode.p₀[begin], v₂)
    lode.ϑ(lode.t₀, lode.q₀[begin], v₁, ϑ₁)
    iode.ϑ(iode.t₀, iode.q₀[begin], v₂, ϑ₂)
    lode.f(lode.t₀, lode.q₀[begin], v₁, f₁)
    iode.f(iode.t₀, iode.q₀[begin], v₂, f₂)
    lode.g(lode.t₀, lode.q₀[begin], v₁, g₁)
    iode.g(iode.t₀, iode.q₀[begin], v₂, g₂)
    @test v₁ == v₂
    @test ϑ₁ == ϑ₂
    @test f₁ == f₂
    @test g₁ == g₂

end
