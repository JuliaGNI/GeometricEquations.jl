
using GeometricEquations
using GeometricBase
using Test

import GeometricEquations: check_noise, ntimesteps

include("functions.jl")
include("initial_conditions.jl")

@testset "$(rpad("Wiener process",80))" begin
    w = WienerProcess(3)

    @test w isa WienerProcess{3}
    @test w isa AbstractStochasticProcess
    @test noisedims(w) == 3

    # The dimension is a type parameter, so it is known to the compiler and costs nothing.
    @test @inferred(noisedims(w)) == 3
    @test Base.return_types(noisedims, (typeof(w),)) == [Int]

    # Singleton semantics: equal processes compare and hash equal.
    @test WienerProcess(3) == w
    @test hash(WienerProcess(3)) == hash(w)
    @test WienerProcess(2) != w
    @test length(Set([WienerProcess(3), w])) == 1

    @test WienerProcess{3}() == w

    @test sprint(show, w) == "Wiener process of dimension 3"

    # A dimension must be a positive Int, whichever constructor is used.
    @test_throws AssertionError WienerProcess(0)
    @test_throws AssertionError WienerProcess(-1)
    @test_throws AssertionError WienerProcess{0}()
    @test_throws AssertionError WienerProcess{:foo}()
    @test_throws AssertionError WienerProcess{1.5}()
end

@testset "$(rpad("Grid process",80))" begin
    ΔW = [1.0 2.0 3.0; 4.0 5.0 6.0]
    ΔZ = [7.0 8.0 9.0; 10.0 11.0 12.0]

    p = GridProcess(ΔW, ΔZ)

    @test p isa AbstractStochasticProcess
    @test p.ΔW == ΔW
    @test p.ΔZ == ΔZ
    @test noisedims(p) == 2
    @test ntime(p) == 3

    @test @inferred(noisedims(p)) == 2
    @test @inferred(ntime(p)) == 3

    @test sprint(show, p) ==
          "Prescribed noise increments of dimension 2 for 3 time steps"

    # The one-argument constructor zeroes the second family of increments.
    p1 = GridProcess(ΔW)

    @test p1.ΔW == ΔW
    @test p1.ΔZ == zero(ΔW)
    @test noisedims(p1) == 2
    @test ntime(p1) == 3

    # ΔW and ΔZ must describe the same grid.
    @test_throws AssertionError GridProcess(ΔW, ΔZ[:, 1:2])
    @test_throws AssertionError GridProcess(ΔW, ΔZ[1:1, :])

    # == and hash agree, so that equal processes are one element of a Set and equal problems
    # built from them hash alike.
    q = GridProcess(copy(ΔW), copy(ΔZ))

    @test p == q
    @test hash(p) == hash(q)
    @test length(Set([p, q])) == 1

    r = GridProcess(ΔW, zero(ΔZ))

    @test p != r
    @test length(Set([p, r])) == 2

    # A GridProcess whose two arrays are of different array types: `zero` of a view or of an
    # adjoint is a plain Matrix, so the one-argument constructor produces exactly this pair.
    v = view(ΔW, :, 1:2)
    pv = GridProcess(v)

    @test pv.ΔW === v
    @test pv.ΔZ == zeros(2, 2)
    @test noisedims(pv) == 2
    @test ntime(pv) == 2

    a = ΔW'
    pa = GridProcess(a)

    @test pa.ΔW === a
    @test pa.ΔZ == zeros(3, 2)
    @test noisedims(pa) == 3
    @test ntime(pa) == 2

    # Mixed element types promote rather than failing to dispatch.
    pm = GridProcess(ΔW, zeros(Int, 2, 3))

    @test eltype(pm.ΔW) == Float64
    @test eltype(pm.ΔZ) == Float64
    @test pm == GridProcess(ΔW)

    pi = GridProcess([1 2; 3 4], [5 6; 7 8])

    @test eltype(pi.ΔW) == Int
    @test ntime(pi) == 2
end

@testset "$(rpad("Noise accessors on equations and problems",80))" begin
    W = WienerProcess(1)

    sde = SDE(sde_v, sde_B, W)
    psde = PSDE(psde_v, psde_f, psde_B, psde_G, W)
    spsde = SPSDE(spsde_v, spsde_f1, spsde_f2, spsde_B, spsde_G1, spsde_G2, W)

    # noise() is answered by the equation, noisedims() forwards to the process.
    for equ in (sde, psde, spsde)
        @test noise(equ) == W
        @test noisedims(equ) == 1
    end

    Δt = 0.1
    nt = 10
    tspan = (0.0, Δt * nt)

    sdeprob = SDEProblem(sde_v, sde_B, W, tspan, Δt, x₀)
    psdeprob = PSDEProblem(psde_v, psde_f, psde_B, psde_G, W, tspan, Δt, q₀, p₀)
    spsdeprob = SPSDEProblem(spsde_v, spsde_f1, spsde_f2, spsde_B, spsde_G1, spsde_G2,
        W, tspan, Δt, q₀, p₀)

    # Both accessors reach through the problem wrapper for all three equation types.
    for prob in (sdeprob, psdeprob, spsdeprob)
        @test noise(prob) == W
        @test noisedims(prob) == 1
        @test @inferred(noisedims(prob)) == 1
        @test Base.return_types(noisedims, (typeof(prob),)) == [Int]
    end

    # A two-dimensional noise is reported as such, and reaches the ensembles too.
    W2 = WienerProcess(2)
    prob2 = SDEProblem(sde_v, sde_B, W2, tspan, Δt, x₀)

    @test noisedims(prob2) == 2

    ens = EnsembleProblem(SDE(sde_v, sde_B, W2), tspan, Δt, ode_ics_tpl)

    @test ens isa SDEEnsemble
    @test noise(ens) == W2
    @test noisedims(ens) == 2

    # A GridProcess answers the same accessors.
    grid = GridProcess(zeros(1, nt))
    gridprob = SDEProblem(sde_v, sde_B, grid, tspan, Δt, x₀)

    @test noise(gridprob) == grid
    @test noisedims(gridprob) == 1
    @test ntime(gridprob) == nt
end

@testset "$(rpad("Noise length check",80))" begin
    Δt = 0.1
    nt = 10
    tspan = (0.0, Δt * nt)

    # The check is against ntime(problem), which is what an integrator steps.
    @test ntimesteps(tspan, Δt) == nt

    # Exactly enough increments, and more than enough, are both fine.
    @test_nowarn SDEProblem(sde_v, sde_B, GridProcess(zeros(1, nt)), tspan, Δt, x₀)
    @test_nowarn SDEProblem(sde_v, sde_B, GridProcess(zeros(1, nt + 1)), tspan, Δt, x₀)

    # Too few is rejected when the problem is built, not when an integrator runs off the end.
    @test_throws ArgumentError SDEProblem(
        sde_v, sde_B, GridProcess(zeros(1, nt - 1)), tspan, Δt, x₀)
    @test_throws ArgumentError SDEProblem(
        sde_v, sde_B, GridProcess(zeros(1, 1)), tspan, Δt, x₀)

    # The same check applies to the partitioned and split equations.
    @test_throws ArgumentError PSDEProblem(psde_v, psde_f, psde_B, psde_G,
        GridProcess(zeros(1, nt - 1)), tspan, Δt, q₀, p₀)
    @test_throws ArgumentError SPSDEProblem(
        spsde_v, spsde_f1, spsde_f2, spsde_B, spsde_G1, spsde_G2,
        GridProcess(zeros(1, nt - 1)), tspan, Δt, q₀, p₀)

    # A WienerProcess draws its increments as it goes and is never short.
    @test_nowarn SDEProblem(sde_v, sde_B, WienerProcess(1), (0.0, 1.0e6), Δt, x₀)

    # Neither is a process that says nothing about its length.
    @test_nowarn SDEProblem(sde_v, sde_B, TestNoise(), tspan, Δt, x₀)
    @test check_noise(TestNoise(), tspan, Δt)

    # Deterministic equations are unaffected by the check.
    @test_nowarn ODEProblem(ode_v, tspan, Δt, x₀)
end
