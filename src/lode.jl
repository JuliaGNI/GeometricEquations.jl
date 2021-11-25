@doc raw"""
`LODE`: Lagrangian Ordinary Differential Equation

Defines an implicit initial value problem
```math
\begin{aligned}
\dot{q} (t) &= v(t) , &
q(t_{0}) &= q_{0} , \\
\dot{p} (t) &= f(t, q(t), v(t)) , &
p(t_{0}) &= p_{0} , \\
p(t) &= ϑ(t, q(t), v(t))
\end{aligned}
```
with momentum ``p`` and force field ``f``, given by
```math
\begin{aligned}
p &= \frac{\partial L}{\partial v} , &
f &= \frac{\partial L}{\partial q} ,
\end{aligned}
```
initial conditions ``(q_{0}, p_{0})`` and the solution ``(q,p)`` taking values
in ``\mathbb{R}^{d} \times \mathbb{R}^{d}``.
This is a special case of a differential algebraic equation with dynamical
variables ``(q,p)`` and algebraic variable ``v``.

### Parameters

* `DT <: Number`: data type
* `TT <: Real`: time step type
* `AT <: AbstractArray{DT}`: array type
* `ϑType <: Function`: type of `ϑ`
* `fType <: Function`: type of `f`
* `gType <: Function`: type of `g`
* `ωType <: Function`: type of `ω`
* `v̄Type <: Function`: type of `v̄`
* `f̄Type <: Function`: type of `f̄`
* `lagType <: Function`: Lagrangian type
* `invType <: OptionalNamedTuple`: invariants type
* `parType <: OptionalNamedTuple`: parameters type
* `perType <: OptionalArray{AT}`: periodicity type

### Fields

* `d`: dimension of dynamical variables ``q`` and ``p`` as well as the vector fields ``f`` and ``p``
* `ϑ`: function determining the momentum
* `f`: function computing the vector field
* `g`: function determining the projection, given by ``\nabla \vartheta (q) \cdot \lambda``
* `ω`: function computing the symplectic matrix
* `v̄`: function computing an initial guess for the velocity field ``v`` (optional)
* `f̄`: function computing an initial guess for the force field ``f`` (optional)
* `t₀`: initial time (optional)
* `q₀`: initial condition for `q`
* `p₀`: initial condition for `p`
* `λ₀`: initial condition for `λ` (optional)
* `lagrangian`: function computing the Lagrangian ``L``
* `invariants`: either a `NamedTuple` containing the equation's invariants or `nothing`
* `parameters`: either a `NamedTuple` containing the equation's parameters or `nothing`
* `periodicity`: determines the periodicity of the state vector `q` for cutting periodic solutions

The functions `ϑ` and `f` must have the interface
```julia
    function ϑ(t, q, v, p)
        p[1] = ...
        p[2] = ...
        ...
    end
```
and
```julia
    function f(t, q, v, f)
        f[1] = ...
        f[2] = ...
        ...
    end
```
where `t` is the current time, `q` is the current solution vector, `v` is the
current velocity and `f` and `p` are the vectors which hold the result of
evaluating the functions ``f`` and ``ϑ`` on `t`, `q` and `v`.
The funtions `g` and `v` are specified by
```julia
    function g(t, q, λ, g)
        g[1] = ...
        g[2] = ...
        ...
    end
```
and
```julia
    function v(t, q, p, v)
        v[1] = ...
        v[2] = ...
        ...
    end
```

### Constructors

```julia
LODE(ϑ, f, ω, v̄, f̄, t₀, q₀, p₀, λ₀, lagrangian, invariants, parameters, periodicity)

LODE(ϑ, f, l, ω, t₀, q₀::StateVector, p₀::StateVector, λ₀::StateVector=zero(q₀); kwargs...)
LODE(ϑ, f, l, ω, q₀::StateVector, p₀::StateVector, λ₀::StateVector=zero(q₀); kwargs...)
LODE(ϑ, f, l, ω, t₀, q₀::State, p₀::State, λ₀::StateVector=zero(q₀); kwargs...)
LODE(ϑ, f, l, ω, q₀::State, p₀::State, λ₀::StateVector=zero(q₀); kwargs...)
```

### Keyword arguments

* `v̄ = (t,q,v) -> nothing`
* `f̄ = f`
* `invariants = nothing`
* `parameters = nothing`
* `periodicity = nothing`

"""
struct LODE{ϑType <: Callable, fType <: Callable, gType <: Callable,
            ωType <: Callable, v̄Type <: Callable, f̄Type <: Callable,
            lagType <: Callable,
            invType <: OptionalInvariants,
            parType <: OptionalParameters,
            perType <: OptionalPeriodicity} <: AbstractEquationPODE{invType,parType,perType}

    ϑ::ϑType
    f::fType
    g::gType
    ω::ωType

    v̄::v̄Type
    f̄::f̄Type

    lagrangian::lagType
    invariants::invType
    parameters::parType
    periodicity::perType

    function LODE(ϑ, f, g, ω, v̄, f̄, lagrangian, invariants, parameters, periodicity)
        @assert !isempty(methods(ϑ))
        @assert !isempty(methods(f))
        @assert !isempty(methods(g))
        @assert !isempty(methods(ω))
        @assert !isempty(methods(v̄))
        @assert !isempty(methods(f̄))
        @assert !isempty(methods(lagrangian))

        new{typeof(ϑ), typeof(f), typeof(g), typeof(ω), typeof(v̄), typeof(f̄),
            typeof(lagrangian), typeof(invariants), typeof(parameters), typeof(periodicity)}(
                ϑ, f, g, ω, v̄, f̄, lagrangian, invariants, parameters, periodicity)
    end
end

_lode_default_v̄(t,q,v,params) = nothing

LODE(ϑ, f, g, ω, l; invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity(), v̄=_lode_default_v̄, f̄=f) = LODE(ϑ, f, g, ω, v̄, f̄, l, invariants, parameters, periodicity)

GeometricBase.invariants(equation::LODE) = equation.invariants
GeometricBase.parameters(equation::LODE) = equation.parameters
GeometricBase.periodicity(equation::LODE) = equation.periodicity

hasvectorfield(::LODE) = true
haslagrangian(::LODE) = true

function check_initial_conditions(::LODE, ics::NamedTuple)
    haskey(ics, :q) || return false
    haskey(ics, :v) || haskey(ics, :p) || return false
    haskey(ics, :λ) || return false
    eltype(ics.q) == eltype(ics.p) == eltype(ics.λ) || return false
    typeof(ics.q) == typeof(ics.p) == typeof(ics.λ) || return false
    axes(ics.q) == axes(ics.p) == axes(ics.λ) || return false
    return true
end

function check_methods(equ::LODE, tspan, ics, params)
    applicable(equ.ϑ, tspan[begin], ics.q, zero(ics.q), zero(ics.p), params) || return false
    applicable(equ.f, tspan[begin], ics.q, zero(ics.q), zero(ics.p), params) || return false
    applicable(equ.g, tspan[begin], ics.q, zero(ics.q), zero(ics.p), params) || return false
    applicable(equ.lagrangian, tspan[begin], ics.q, zero(ics.q), params) || return false
    return true
end

function datatype(equ::LODE, ics::NamedTuple)
    @assert check_initial_conditions(equ, ics)
    return eltype(ics.q)
end

function arrtype(equ::LODE, ics::NamedTuple)
    @assert check_initial_conditions(equ, ics)
    typeof(ics.q)
end

_get_ϑ(equ::LODE, params) = (t,q,v,ϑ) -> equ.ϑ(t, q, v, ϑ, params)
_get_f(equ::LODE, params) = (t,q,v,f) -> equ.f(t, q, v, f, params)
_get_g(equ::LODE, params) = (t,q,v,g) -> equ.g(t, q, v, g, params)
_get_ω(equ::LODE, params) = (t,q,v,ω) -> equ.ω(t, q, v, ω, params)
_get_v̄(equ::LODE, params) = (t,q,v)   -> equ.v̄(t, q, v, params)
_get_f̄(equ::LODE, params) = (t,q,v,f) -> equ.f̄(t, q, v, f, params)
_get_l(equ::LODE, params) = (t,q,v)   -> equ.lagrangian(t, q, v, params)

_functions(equ::LODE) = (ϑ = equ.ϑ, f = equ.f, g = equ.g, ω = equ.ω, v̄ = equ.v̄, f̄ = equ.f̄, l = equ.lagrangian)
_functions(equ::LODE, params::OptionalParameters) = (ϑ = _get_ϑ(equ, params),
                                                     f = _get_f(equ, params),
                                                     g = _get_g(equ, params),
                                                     ω = _get_ω(equ, params),
                                                     v̄ = _get_v̄(equ, params),
                                                     f̄ = _get_f̄(equ, params),
                                                     l = _get_l(equ, params))
