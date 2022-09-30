@doc raw"""
`IDAE`: Implicit Differential Algebraic Equation

Defines an implicit differential algebraic initial value problem
```math
\begin{aligned}
\dot{q} (t) &= v(t) + u(t, q(t), p(t), \lambda(t)) , & q(t_{0}) &= q_{0} , \\
\dot{p} (t) &= f(t, q(t), v(t)) + g(t, q(t), p(t), \lambda(t)) , & p(t_{0}) &= p_{0} , \\
p(t) &= p(t, q(t), v(t)) , && \\
0 &= \phi (t, q(t), p(t)) , & \lambda(t_{0}) &= \lambda_{0} ,
\end{aligned}
```
with force field ``f``, the momentum defined by ``p``, projection ``g`` and ``r``,
algebraic constraint ``\phi=0``,
conditions ``(q_{0}, p_{0})`` and ``\lambda_{0}``, the dynamical variables
``(q,p)`` taking values in ``\mathbb{R}^{d} \times \mathbb{R}^{d}`` and
the algebraic variable ``\lambda`` taking values in ``\mathbb{R}^{m}``.

### Parameters

* `DT <: Number`: data type
* `TT <: Real`: time step type
* `AT <: AbstractArray{DT}`: array type
* `ϑType <: Function`: type of `ϑ`
* `fType <: Function`: type of `f`
* `uType <: Function`: type of `u`
* `gType <: Function`: type of `g`
* `ϕType <: Function`: type of `ϕ`
* `ūType <: Function`: type of `ū`
* `ḡType <: Function`: type of `ḡ`
* `ψType <: Function`: type of `ψ`
* `v̄Type <: Function`: type of `v̄`
* `f̄Type <: Function`: type of `f̄`
* `invType <: OptionalNamedTuple`: invariants type
* `parType <: OptionalNamedTuple`: parameters type
* `perType <: OptionalArray{AT}`: periodicity type

### Fields

* `d`: dimension of dynamical variables ``q`` and ``p`` as well as the vector fields ``f`` and ``p``
* `m`: dimension of algebraic variable ``\lambda`` and the constraint ``\phi``
* `ϑ`: function determining the momentum
* `f`: function computing the vector field ``f``
* `u`: function computing the projection for ``q``
* `g`: function computing the projection for ``p``
* `ϕ`: algebraic constraints
* `ū`: function computing the secondary projection field ``\bar{u}`` (optional)
* `ḡ`: function computing the secondary projection field ``\bar{g}`` (optional)
* `ψ`: secondary constraints (optional)
* `v̄`: function computing an initial guess for the velocity field ``v`` (optional)
* `f̄`: function computing an initial guess for the force field ``f`` (optional)
* `t₀`: initial time (optional)
* `q₀`: initial condition for dynamical variable ``q``
* `p₀`: initial condition for dynamical variable ``p``
* `λ₀`: initial condition for algebraic variable `λ`
* `μ₀`: initial condition for algebraic variable `μ` (optional)
* `invariants`: either a `NamedTuple` containing the equation's invariants or `nothing`
* `parameters`: either a `NamedTuple` containing the equation's parameters or `nothing`
* `periodicity`: determines the periodicity of the state vector `q` for cutting periodic solutions

### Constructors

```julia
IDAE(ϑ, f, u, g, ϕ, ū, ḡ, ψ, v̄, f̄, t₀, q₀, p₀, λ₀, μ₀, invariants, parameters, periodicity)

IDAE(ϑ, f, u, g, ϕ, t₀, q₀::StateVector, p₀::StateVector, λ₀::StateVector, μ₀::StateVector=zero(λ₀); kwargs...)
IDAE(ϑ, f, u, g, ϕ, q₀::StateVector, p₀::StateVector, λ₀::StateVector, μ₀::StateVector=zero(λ₀); kwargs...)
IDAE(ϑ, f, u, g, ϕ, t₀, q₀::State, p₀::State, λ₀::State, μ₀::State=zero(λ₀); kwargs...)
IDAE(ϑ, f, u, g, ϕ, q₀::State, p₀::State, λ₀::State, μ₀::State=zero(λ₀); kwargs...)

IDAE(ϑ, f, u, g, ϕ, ū, ḡ, ψ, t₀, q₀::StateVector, p₀::StateVector, λ₀::StateVector, μ₀::StateVector=zero(λ₀); kwargs...)
IDAE(ϑ, f, u, g, ϕ, ū, ḡ, ψ, q₀::StateVector, p₀::StateVector, λ₀::StateVector, μ₀::StateVector=zero(λ₀); kwargs...)
IDAE(ϑ, f, u, g, ϕ, ū, ḡ, ψ, t₀, q₀::State, p₀::State, λ₀::State, μ₀::State=zero(λ₀); kwargs...)
IDAE(ϑ, f, u, g, ϕ, ū, ḡ, ψ, q₀::State, p₀::State, λ₀::State, μ₀::State=zero(λ₀); kwargs...)
```

### Keyword arguments

* `v̄ = (t,q,v) -> nothing`
* `f̄ = f`
* `invariants = nothing`
* `parameters = nothing`
* `periodicity = nothing`


"""
struct IDAE{ϑType <: Callable, fType <: Callable,
            uType <: Callable, gType <: Callable, ϕType <: Callable,
            ūType <: OptionalCallable, ḡType <: OptionalCallable, ψType <: OptionalCallable,
            v̄Type <: Callable, f̄Type <: Callable,
            invType <: OptionalInvariants,
            parType <: OptionalParameters,
            perType <: OptionalPeriodicity} <: AbstractEquationPDAE{invType,parType,perType,ψType}

    ϑ::ϑType
    f::fType
    u::uType
    g::gType
    ϕ::ϕType
    ū::ūType
    ḡ::ḡType
    ψ::ψType
    v̄::v̄Type
    f̄::f̄Type

    invariants::invType
    parameters::parType
    periodicity::perType

    function IDAE(ϑ, f, u, g, ϕ, ū, ḡ, ψ, v̄, f̄, invariants, parameters, periodicity)
        @assert !isempty(methods(ϑ))
        @assert !isempty(methods(f))
        @assert !isempty(methods(u))
        @assert !isempty(methods(g))
        @assert !isempty(methods(ϕ))
        @assert !isempty(methods(ū)) || ū === nothing
        @assert !isempty(methods(ḡ)) || ḡ === nothing
        @assert !isempty(methods(ψ)) || ψ === nothing
        @assert !isempty(methods(v̄))
        @assert !isempty(methods(f̄))

        new{typeof(ϑ), typeof(f),
            typeof(u), typeof(g), typeof(ϕ),
            typeof(ū), typeof(ḡ), typeof(ψ),
            typeof(v̄), typeof(f̄),
            typeof(invariants), typeof(parameters), typeof(periodicity)}(
                ϑ, f, u, g, ϕ, ū, ḡ, ψ, v̄, f̄, invariants, parameters, periodicity)
    end
end

_idae_default_v̄(t,q,v,params) = nothing

IDAE(ϑ, f, u, g, ϕ, ū, ḡ, ψ; v̄=_idae_default_v̄, f̄=f, invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity()) = IDAE(ϑ, f, u, g, ϕ, ū, ḡ, ψ, v̄, f̄, invariants, parameters, periodicity)
IDAE(ϑ, f, u, g, ϕ; kwargs...) = IDAE(ϑ, f, u, g, ϕ, nothing, nothing, nothing; kwargs...)

GeometricBase.invariants(equation::IDAE) = equation.invariants
GeometricBase.parameters(equation::IDAE) = equation.parameters
GeometricBase.periodicity(equation::IDAE) = equation.periodicity

hasvectorfield(::IDAE) = true

function check_initial_conditions(equ::IDAE, ics::NamedTuple)
    haskey(ics, :q) || return false
    haskey(ics, :p) || return false
    haskey(ics, :λ) || return false
    if hassecondary(equ)
        haskey(ics, :μ) || return false
        eltype(ics.λ) == eltype(ics.μ) || return false
        typeof(ics.λ) == typeof(ics.μ) || return false
        axes(ics.λ) == axes(ics.μ) || return false
    end
    return true
end

function check_methods(equ::IDAE, tspan, ics::NamedTuple, params)
    applicable(equ.ϑ, zero(ics.p), tspan[begin], ics.q, vectorfield(ics.q), params) || return false
    applicable(equ.f, zero(ics.p), tspan[begin], ics.q, vectorfield(ics.q), params) || return false
    applicable(equ.u, zero(ics.q), tspan[begin], ics.q, ics.p, ics.λ, params) || return false
    applicable(equ.g, zero(ics.p), tspan[begin], ics.q, ics.p, ics.λ, params) || return false
    applicable(equ.ϕ, zero(ics.λ), tspan[begin], ics.q, ics.p, params) || return false
    applicable(equ.v̄, zero(ics.q), tspan[begin], ics.q, params) || return false
    applicable(equ.f̄, zero(ics.p), tspan[begin], ics.q, vectorfield(ics.q), params) || return false
    equ.ū === nothing || applicable(equ.ū, zero(ics.q), tspan[begin], ics.q, ics.p, ics.λ, params) || return false
    equ.ḡ === nothing || applicable(equ.ḡ, zero(ics.p), tspan[begin], ics.q, ics.p, ics.λ, params) || return false
    # equ.ψ === nothing || applicable(equ.ψ, zero(ics.λ), tspan[begin], ics.q, ics.p, vectorfield(ics.q), vectorfield(ics.p), params) || return false
    # TODO add missing methods
    return true
end

function GeometricBase.datatype(equ::IDAE, ics::NamedTuple)
    @assert check_initial_conditions(equ, ics)
    return eltype(ics.q)
end

function GeometricBase.arrtype(equ::IDAE, ics::NamedTuple)
    @assert check_initial_conditions(equ, ics)
    return typeof(ics.q)
end

_get_ϑ(equ::IDAE, params) = (ϑ, t, q, v)       -> equ.ϑ(ϑ, t, q, v, params)
_get_f(equ::IDAE, params) = (f, t, q, v)       -> equ.f(f, t, q, v, params)
_get_u(equ::IDAE, params) = (u, t, q, p, λ)    -> equ.u(u, t, q, p, λ, params)
_get_g(equ::IDAE, params) = (g, t, q, p, λ)    -> equ.g(g, t, q, p, λ, params)
_get_ϕ(equ::IDAE, params) = (ϕ, t, q, p)       -> equ.ϕ(ϕ, t, q, p, params)
_get_ū(equ::IDAE, params) = (u, t, q, p, λ)    -> equ.ū(u, t, q, p, λ, params)
_get_ḡ(equ::IDAE, params) = (g, t, q, p, λ)    -> equ.ḡ(g, t, q, p, λ, params)
_get_ψ(equ::IDAE, params) = (ψ, t, q, p, v, f) -> equ.ψ(ψ, t, q, p, v, f, params)
_get_v̄(equ::IDAE, params) = (v, t, q)          -> equ.v̄(v, t, q, params)
_get_f̄(equ::IDAE, params) = (f, t, q, v)       -> equ.f̄(f, t, q, v, params)
_get_invariant(::IDAE, inv, params) = (t, q, v) -> inv(t, q, v, params)

function _functions(equ::IDAE)
    if hassecondary(equ)
        (ϑ = equ.ϑ, f = equ.f, u = equ.u, g = equ.g, ϕ = equ.ϕ, ū = equ.ū, ḡ = equ.ḡ, ψ = equ.ψ, v̄ = equ.v̄, f̄ = equ.f̄)
    else
        (ϑ = equ.ϑ, f = equ.f, u = equ.u, g = equ.g, ϕ = equ.ϕ, v̄ = equ.v̄, f̄ = equ.f̄)
    end
end

function _functions(equ::IDAE, params::OptionalParameters)
    if hassecondary(equ)
        (ϑ = _get_ϑ(equ, params), f = _get_f(equ, params), u = _get_u(equ, params), g = _get_g(equ, params), ϕ = _get_ϕ(equ, params), ū = _get_ū(equ, params), ḡ = _get_ḡ(equ, params), ψ = _get_ψ(equ, params), v̄ = _get_v̄(equ, params), f̄ = _get_f̄(equ, params))
    else
        (ϑ = _get_ϑ(equ, params), f = _get_f(equ, params), u = _get_u(equ, params), g = _get_g(equ, params), ϕ = _get_ϕ(equ, params), v̄ = _get_v̄(equ, params), f̄ = _get_f̄(equ, params))
    end
end
