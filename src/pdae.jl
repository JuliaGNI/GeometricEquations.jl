@doc raw"""
`PDAE`: Partitioned Differential Algebraic Equation

Defines a partitioned differential algebraic initial value problem
```math
\begin{aligned}
\dot{q} (t) &= v(t, q(t), p(t)) + u(t, q(t), p(t), \lambda(t)) , & q(t_{0}) &= q_{0} , \\
\dot{p} (t) &= f(t, q(t), p(t)) + g(t, q(t), p(t), \lambda(t)) , & p(t_{0}) &= p_{0} , \\
0 &= \phi (t, q(t), p(t), \lambda(t)) , & \lambda(t_{0}) &= \lambda_{0} ,
\end{aligned}
```
with vector fields ``v`` and ``f``, projection ``u`` and ``g``,
algebraic constraint ``\phi=0``,
conditions ``(q_{0}, p_{0})`` and ``\lambda_{0}``, the dynamical variables
``(q,p)`` taking values in ``\mathbb{R}^{d} \times \mathbb{R}^{d}`` and
the algebraic variable ``\lambda`` taking values in ``\mathbb{R}^{m}``.

### Parameters

* `DT <: Number`: data type
* `TT <: Real`: time step type
* `AT <: AbstractArray{DT}`: array type
* `vType <: Function`: type of `v`
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
* `v`: function computing the vector field ``v``
* `f`: function computing the vector field ``f``
* `u`: function computing the projection for ``q``
* `g`: function computing the projection for ``p``
* `ϕ`: algebraic constraints
* `ū`: function computing the secondary projection field ``\bar{u}`` (optional)
* `ḡ`: function computing the secondary projection field ``\bar{g}`` (optional)
* `ψ`: secondary constraints (optional)
* `v̄`: function computing an initial guess for the velocity field ``v`` (optional)
* `f̄`: function computing an initial guess for the force field ``f`` (optional)
* `t₀`: initial time
* `q₀`: initial condition for dynamical variable ``q``
* `p₀`: initial condition for dynamical variable ``p``
* `λ₀`: initial condition for algebraic variable ``\lambda``
* `μ₀`: initial condition for algebraic variable ``μ`` (optional)
* `invariants`: either a `NamedTuple` containing the equation's invariants or `nothing`
* `parameters`: either a `NamedTuple` containing the equation's parameters or `nothing`
* `periodicity`: determines the periodicity of the state vector `q` for cutting periodic solutions

### Constructors

```julia
PDAE(v, f, u, g, ϕ, ū, ḡ, ψ, v̄, f̄, t₀, q₀, p₀, λ₀, invariants, parameters, periodicity)

PDAE(v, f, u, g, ϕ, t₀, q₀::StateVector, p₀::StateVector, λ₀::StateVector, μ₀::StateVector=zero(λ₀); kwargs...)
PDAE(v, f, u, g, ϕ, q₀::StateVector, p₀::StateVector, λ₀::StateVector, μ₀::StateVector=zero(λ₀); kwargs...)
PDAE(v, f, u, g, ϕ, t₀, q₀::State, p₀::State, λ₀::State, μ₀::State=zero(λ₀); kwargs...)
PDAE(v, f, u, g, ϕ, q₀::State, p₀::State, λ₀::State, μ₀::State=zero(λ₀); kwargs...)

PDAE(v, f, u, g, ϕ, ū, ḡ, ψ, t₀, q₀::StateVector, p₀::StateVector, λ₀::StateVector, μ₀::StateVector=zero(λ₀); kwargs...)
PDAE(v, f, u, g, ϕ, ū, ḡ, ψ, q₀::StateVector, p₀::StateVector, λ₀::StateVector, μ₀::StateVector=zero(λ₀); kwargs...)
PDAE(v, f, u, g, ϕ, ū, ḡ, ψ, t₀, q₀::State, p₀::State, λ₀::State, μ₀::State=zero(λ₀); kwargs...)
PDAE(v, f, u, g, ϕ, ū, ḡ, ψ, q₀::State, p₀::State, λ₀::State, μ₀::State=zero(λ₀); kwargs...)
```

"""
struct PDAE{vType <: Callable, fType <: Callable,
            uType <: Callable, gType <: Callable, ϕType <: Callable,
            ūType <: OptionalCallable, ḡType <: OptionalCallable, ψType <: OptionalCallable,
            v̄Type <: Callable, f̄Type <: Callable,
            invType <: OptionalInvariants,
            parType <: OptionalParameters,
            perType <: OptionalPeriodicity} <: AbstractEquationPDAE{invType,parType,perType,ψType}

    v::vType
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

    function PDAE(v, f, u, g, ϕ, ū, ḡ, ψ, v̄, f̄, invariants, parameters, periodicity)
        @assert !isempty(methods(v))
        @assert !isempty(methods(f))
        @assert !isempty(methods(u))
        @assert !isempty(methods(g))
        @assert !isempty(methods(ϕ))
        @assert !isempty(methods(ū)) || ū === nothing
        @assert !isempty(methods(ḡ)) || ḡ === nothing
        @assert !isempty(methods(ψ)) || ψ === nothing
        @assert !isempty(methods(v̄))
        @assert !isempty(methods(f̄))

        new{typeof(v), typeof(f),
            typeof(u), typeof(g), typeof(ϕ),
            typeof(ū), typeof(ḡ), typeof(ψ),
            typeof(v̄), typeof(f̄),
            typeof(invariants), typeof(parameters), typeof(periodicity)}(
                v, f, u, g, ϕ, ū, ḡ, ψ, v̄, f̄, invariants, parameters, periodicity)
    end
end

PDAE(v, f, u, g, ϕ, ū, ḡ, ψ; v̄=v, f̄=f, invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity()) = PDAE(v, f, u, g, ϕ, ū, ḡ, ψ, v̄, f̄, invariants, parameters, periodicity)
PDAE(v, f, u, g, ϕ; kwargs...) = PDAE(v, f, u, g, ϕ, nothing, nothing, nothing; kwargs...)

GeometricBase.invariants(equation::PDAE) = equation.invariants
GeometricBase.parameters(equation::PDAE) = equation.parameters
GeometricBase.periodicity(equation::PDAE) = equation.periodicity

hasvectorfield(::PDAE) = true

function check_initial_conditions(equ::PDAE, ics::NamedTuple)
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

function check_methods(equ::PDAE, tspan, ics::NamedTuple, params)
    applicable(equ.v, tspan[begin], ics.q, ics.p, zero(ics.q), params) || return false
    applicable(equ.f, tspan[begin], ics.q, ics.p, zero(ics.p), params) || return false
    applicable(equ.ϕ, tspan[begin], ics.q, ics.p, zero(ics.λ), params) || return false
    # TODO add missing methods
    return true
end

function datatype(equ::PDAE, ics::NamedTuple)
    @assert check_initial_conditions(equ, ics)
    return eltype(ics.q)
end

function arrtype(equ::PDAE, ics::NamedTuple)
    @assert check_initial_conditions(equ, ics)
    return typeof(ics.q)
end

_get_v(equ::PDAE, params) = (t,q,p,v)     -> equ.v(t, q, p, v, params)
_get_f(equ::PDAE, params) = (t,q,p,f)     -> equ.f(t, q, p, f, params)
_get_u(equ::PDAE, params) = (t,q,p,λ,u)   -> equ.u(t, q, p, λ, u, params)
_get_g(equ::PDAE, params) = (t,q,p,λ,g)   -> equ.g(t, q, p, λ, g, params)
_get_ϕ(equ::PDAE, params) = (t,q,p,ϕ)     -> equ.ϕ(t, q, p, ϕ, params)
_get_ū(equ::PDAE, params) = (t,q,p,λ,u)   -> equ.ū(t, q, p, λ, u, params)
_get_ḡ(equ::PDAE, params) = (t,q,p,λ,g)   -> equ.ḡ(t, q, p, λ, g, params)
_get_ψ(equ::PDAE, params) = (t,q,p,v,f,ψ) -> equ.ψ(t, q, p, v, f, ψ, params)
_get_v̄(equ::PDAE, params) = (t,q,p,v)     -> equ.v̄(t, q, p, v, params)
_get_f̄(equ::PDAE, params) = (t,q,p,f)     -> equ.f̄(t, q, p, f, params)
_get_invariant(::PDAE, inv, params) = (t,q,p) -> inv(t, q, p, params)

function _functions(equ::PDAE)
    if hassecondary(equ)
        (v = equ.v, f = equ.f, u = equ.u, g = equ.g, ϕ = equ.ϕ, ū = equ.ū, ḡ = equ.ḡ, ψ = equ.ψ, v̄ = equ.v̄, f̄ = equ.f̄)
    else
        (v = equ.v, f = equ.f, u = equ.u, g = equ.g, ϕ = equ.ϕ, v̄ = equ.v̄, f̄ = equ.f̄)
    end
end

function _functions(equ::PDAE, params::OptionalParameters)
    if hassecondary(equ)
        (v = _get_v(equ, params), f = _get_f(equ, params), u = _get_u(equ, params), g = _get_g(equ, params), ϕ = _get_ϕ(equ, params), ū = _get_ū(equ, params), ḡ = _get_ḡ(equ, params), ψ = _get_ψ(equ, params), v̄ = _get_v̄(equ, params), f̄ = _get_f̄(equ, params))
    else
        (v = _get_v(equ, params), f = _get_f(equ, params), u = _get_u(equ, params), g = _get_g(equ, params), ϕ = _get_ϕ(equ, params), v̄ = _get_v̄(equ, params), f̄ = _get_f̄(equ, params))
    end
end
