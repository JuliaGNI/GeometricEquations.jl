@doc raw"""
`SPDAE`: Split Partitioned Differential Algebraic Equation *EXPERIMENTAL*

Defines a split differential algebraic initial value problem, that is
a canonical Hamiltonian system of equations subject to Dirac constraints,
```math
\begin{aligned}
\dot{q} (t) &= v_1(t, q(t), p(t)) + v_2(t, q(t), p(t), \lambda(t)) + v_3(t, q(t), p(t), \lambda(t), \gamma(t)) , & q(t_{0}) &= q_{0} , \\
\dot{p} (t) &= f_1(t, q(t), p(t)) + f_2(t, q(t), p(t), \lambda(t)) + f_3(t, q(t), p(t), \lambda(t), \gamma(t)) , & p(t_{0}) &= p_{0} , \\
0 &= \phi (t, q(t), p(t)) , \\
0 &= \psi (t, q(t), p(t), \dot{q}(t), \dot{p}(t)) ,
\end{aligned}
```
with vector fields ``v_i`` and ``f_i`` for ``i = 1 ... 3``,
primary constraint ``\phi(q,p)=0`` and secondary constraint ``\psi(q,p,\lambda)=0``,
initial conditions ``(q_{0}, p_{0})``, the dynamical variables ``(q,p)``
taking values in ``\mathbb{R}^{d} \times \mathbb{R}^{d}`` and
the algebraic variables ``(\lambda, \gamma)`` taking values in
``\mathbb{R}^{n} \times \mathbb{R}^{d}``.

### Parameters

* `vType <: Callable`: type of `v`
* `fType <: Callable`: type of `f`
* `ϕType <: Callable`: type of `ϕ`
* `ψType <: Callable`: type of `ψ`
* `invType <: OptionalInvariants`: invariants type
* `parType <: OptionalParameters`: parameters type
* `perType <: OptionalPeriodicity`: periodicity type

### Fields

* `v`: tuple of functions computing the vector fields ``v_i``, ``i = 1 ... 3``
* `f`: tuple of functions computing the vector fields ``f_i``, ``i = 1 ... 3``
* `ϕ`: primary constraints
* `ψ`: secondary constraints
* `invariants`: functions for the computation of invariants, either a `NamedTuple` containing the equation's invariants or `NullInvariants`
* `parameters`: type constraints for parameters, either a `NamedTuple` containing the equation's parameters or `NullParameters`
* `periodicity`: determines the periodicity of the state vector `q` for cutting periodic solutions, either a `AbstractArray` or `NullPeriodicity`

### Constructors

```julia
SPDAE(v, f, ϕ, ψ, t₀, q₀::StateVector, p₀::StateVector, λ₀::StateVector; kwargs...)
SPDAE(v, f, ϕ, ψ, q₀::StateVector, p₀::StateVector, λ₀::StateVector; kwargs...)
SPDAE(v, f, ϕ, ψ, t₀, q₀::State, p₀::State, λ₀::State=zero(q₀); kwargs...)
SPDAE(v, f, ϕ, ψ, q₀::State, p₀::State, λ₀::State=zero(q₀); kwargs...)
```

### Keyword arguments:

* `invariants = nothing`
* `parameters = nothing`
* `periodicity = nothing`

"""
struct SPDAE{dType <: Number, tType <: Real, arrayType <: AbstractArray{dType},
             vType <: Tuple, fType <: Tuple, ϕType <: Function, ψType <: OptionalFunction,
             invType <: OptionalInvariants,
             parType <: OptionalParameters,
             perType <: OptionalPeriodicity} <: AbstractEquationPDAE{invType,parType,perType,ψType}

    d::Int
    m::Int

    v::vType
    f::fType
    ϕ::ϕType
    ψ::ψType

    t₀::tType
    q₀::Vector{arrayType}
    p₀::Vector{arrayType}
    λ₀::Vector{arrayType}
    μ₀::Vector{arrayType}

    invariants::invType
    parameters::parType
    periodicity::perType

    function SPDAE(v::vType, f::fType, ϕ::ϕType, ψ::ψType,
                   t₀::tType, q₀::Vector{arrayType}, p₀::Vector{arrayType}, λ₀::Vector{arrayType}, μ₀::Vector{arrayType},
                   invariants::invType, parameters::parType, periodicity::perType) where {
                        dType <: Number, tType <: Real, arrayType <: AbstractArray{dType},
                        vType <: Tuple, fType <: Tuple, ϕType <: Function, ψType <: Function,
                        invType <: OptionalInvariants,
                        parType <: OptionalParameters,
                        perType <: OptionalPeriodicity}

        d = length(q₀[begin])
        m = length(λ₀[begin])

        @assert 2d ≥ m

        @assert length(q₀) == length(p₀)
        @assert length(λ₀) == length(μ₀)

        @assert all(length(q) == d for q in q₀)
        @assert all(length(p) == d for p in p₀)
        @assert all(length(λ) == m for λ in λ₀)
        @assert all(length(μ) == m for μ in μ₀)

        @assert all([ndims(q) == ndims(p) == ndims(λ) == ndims(μ) for (q,p,λ,μ) in zip(q₀,p₀,λ₀,μ₀)])

        new{dType, tType, arrayType, vType, fType, ϕType, ψType, invType, parType, perType}(d, m, v, f, ϕ, ψ, t₀, q₀, p₀, λ₀, μ₀, invariants, parameters, periodicity)
    end
end

_SPDAE(v, f, ϕ, ψ, t₀, q₀, p₀, λ₀, μ₀; invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity()) = SPDAE(v, f, ϕ, ψ, t₀, q₀, p₀, λ₀, μ₀, invariants, parameters, periodicity)

SPDAE(v, f, ϕ, ψ, t₀, q₀::StateVector, p₀::StateVector, λ₀::StateVector, μ₀::StateVector=zero(λ₀); kwargs...) = _SPDAE(v, f, ϕ, ψ, t₀, q₀, p₀, λ₀, μ₀; kwargs...)
SPDAE(v, f, ϕ, ψ, q₀::StateVector, p₀::StateVector, λ₀::StateVector, μ₀::StateVector=zero(λ₀); kwargs...) = SPDAE(v, f, ϕ, ψ, 0.0, q₀, p₀, λ₀, μ₀; kwargs...)
SPDAE(v, f, ϕ, ψ, t₀, q₀::State, p₀::State, λ₀::State, μ₀::State=zero(λ₀); kwargs...) = SPDAE(v, f, ϕ, ψ, t₀, [q₀], [p₀], [λ₀], [μ₀]; kwargs...)
SPDAE(v, f, ϕ, ψ, q₀::State, p₀::State, λ₀::State, μ₀::State=zero(λ₀); kwargs...) = SPDAE(v, f, ϕ, ψ, 0.0, q₀, p₀, λ₀, μ₀; kwargs...)

Base.hash(dae::SPDAE, h::UInt) = hash(dae.d, hash(dae.m,
                    hash(dae.v, hash(dae.f, hash(dae.ϕ, hash(dae.ψ,
                    hash(dae.t₀, hash(dae.q₀, hash(dae.p₀, 
                    hash(dae.invariants, hash(dae.parameters, hash(dae.periodicity, h))))))))))))

Base.:(==)(dae1::SPDAE, dae2::SPDAE) = (
                                dae1.d == dae2.d
                             && dae1.m == dae2.m
                             && dae1.v == dae2.v
                             && dae1.f == dae2.f
                             && dae1.ϕ == dae2.ϕ
                             && dae1.ψ == dae2.ψ
                             && dae1.t₀ == dae2.t₀
                             && dae1.q₀ == dae2.q₀
                             && dae1.p₀ == dae2.p₀
                             && dae1.λ₀ == dae2.λ₀
                             && dae1.μ₀ == dae2.μ₀
                             && dae1.invariants  == dae2.invariants
                             && dae1.parameters  == dae2.parameters
                             && dae1.periodicity == dae2.periodicity)

function Base.similar(equ::SPDAE, t₀::Real, q₀::StateVector, p₀::StateVector, λ₀::StateVector, μ₀::StateVector=initial_multiplier(q₀, equ.μ₀); parameters=equ.parameters)
    @assert all([length(q) == equ.d for q in q₀])
    @assert all([length(p) == equ.d for p in p₀])
    @assert all([length(λ) == equ.m for λ in λ₀])
    @assert all([length(μ) == equ.m for μ in μ₀])
    _SPDAE(equ.v, equ.f, equ.ϕ, equ.ψ, t₀, q₀, p₀, λ₀, μ₀;
          invariants=equ.invariants, parameters=parameters, periodicity=equ.periodicity)
end

Base.similar(equ::SPDAE, q₀, p₀, λ₀=initial_multiplier(q₀, equ.λ₀), μ₀=initial_multiplier(q₀, equ.μ₀); kwargs...) = similar(equ, equ.t₀, q₀, p₀, λ₀, μ₀; kwargs...)
Base.similar(equ::SPDAE, t₀::Real, q₀::State, p₀::State, λ₀::State=initial_multiplier(q₀, equ.λ₀), μ₀::State=initial_multiplier(λ₀, equ.μ₀); kwargs...) = similar(equ, equ.t₀, [q₀], [p₀], [λ₀], [μ₀]; kwargs...)

@inline Base.axes(equation::SPDAE) = axes(equation.q₀[begin])
@inline Base.ndims(equation::SPDAE) = equation.d
@inline GeometricBase.nsamples(equation::SPDAE) = length(eachindex(equation.q₀))
@inline GeometricBase.nconstraints(equation::SPDAE) = equation.m

@inline GeometricBase.periodicity(equation::SPDAE) = hasperiodicity(equation) ? equation.periodicity : zero(equation.q₀[begin])
@inline initial_conditions(equation::SPDAE) = (equation.t₀, equation.q₀, equation.p₀, equation.λ₀, equation.μ₀)

function functions(equation::SPDAE{DT,TT,AT,VT,FT,ϕT,ψT,IT,PT}) where {DT, TT, AT, VT, FT, ϕT, ψT, IT, PT <: NullParameters}
    NamedTuple{(:v, :f, :ϕ, :ψ)}((equation.v, equation.f, equation.ϕ, equation.ψ))
end

function functions(equation::SPDAE{DT,TT,AT,VT,FT,ϕT,ψT,IT,PT}) where {DT, TT, AT, VT, FT, ϕT, ψT, IT, PT <: NamedTuple}
    vₚ = (v, t, q, p)   -> equation.v(v, t, q, p, equation.parameters)
    fₚ = (f, t, q, p)   -> equation.f(f, t, q, p, equation.parameters)
    ϕₚ = (ϕ, t, q, p)   -> equation.ϕ(ϕ, t, q, p, equation.parameters)
    ψₚ = (ψ, t, q, p, v, f) -> equation.ψ(ψ, t, q, p, v, f, equation.parameters)

    names = (:v, :f, :ϕ, :ψ)
    equs  = (vₚ, fₚ, ϕₚ, ψₚ)

    NamedTuple{names}(equs)
end


"""
"""
const SPDAEProblem = GeometricProblem{SPDAE}

function GeometricBase.periodicity(prob::SPDAEProblem)
    (q = periodicity(equation(prob)), p = NullPeriodicity(), λ = NullPeriodicity())
end
