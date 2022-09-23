@doc raw"""
`IODE`: Implicit Ordinary Differential Equation

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
with force field ``f``, the momentum defined by ``p``, initial conditions ``(q_{0}, p_{0})`` and the solution
``(q,p)`` taking values in ``\mathbb{R}^{d} \times \mathbb{R}^{d}``.
This is a special case of a differential algebraic equation with dynamical
variables ``(q,p)`` and algebraic variable ``v``.

### Parameters

* `DT <: Number`: data type
* `TT <: Real`: time step type
* `AT <: AbstractArray{DT}`: array type
* `ϑType <: Function`: type of `ϑ`
* `fType <: Function`: type of `f`
* `gType <: Function`: type of `g`
* `v̄Type <: Function`: type of `v̄`
* `f̄Type <: Function`: type of `f̄`
* `hType <: OptionalFunction`: type of `h`
* `pType <: Union{NamedTuple,Nothing}`: parameters type

### Fields

* `d`: dimension of dynamical variables ``q`` and ``p`` as well as the vector fields ``f`` and ``p``
* `ϑ`: function determining the momentum
* `f`: function computing the vector field
* `g`: function determining the projection, given by ``\nabla \vartheta (q) \cdot \lambda``
* `v̄`: function computing an initial guess for the velocity field ``v`` (optional)
* `f̄`: function computing an initial guess for the force field ``f`` (optional)
* `h`: function computing the Hamiltonian (optional)
* `t₀`: initial time (optional)
* `q₀`: initial condition for `q`
* `p₀`: initial condition for `p`
* `parameters`: either a `NamedTuple` containing the equations parameters or `nothing`
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
In addition, the functions `g`, `v̄` and `f̄` are specified by
```julia
    function g(t, q, λ, g)
        g[1] = ...
        g[2] = ...
        ...
    end

    function v̄(t, q, v)
        v[1] = ...
        v[2] = ...
        ...
    end

    function f̄(t, q, v, f)
        f[1] = ...
        f[2] = ...
        ...
    end
```
The function `g` is used in projection methods that enforce ``p = ϑ(q)``.
The functions `v̄` and `f̄` are used for initial guesses in nonlinear implicit solvers.

### Constructors

```julia
IODE(ϑ, f, v̄, f̄, t₀, q₀, p₀, λ₀, invariants, parameters, periodicity)

IODE(ϑ, f, t₀, q₀::StateVector, p₀::StateVector, λ₀::StateVector=zero(q₀); kwargs...)
IODE(ϑ, f, q₀::StateVector, p₀::StateVector, λ₀::StateVector=zero(q₀); kwargs...)
IODE(ϑ, f, t₀, q₀::State, p₀::State, λ₀::StateVector=zero(q₀); kwargs...)
IODE(ϑ, f, q₀::State, p₀::State, λ₀::StateVector=zero(q₀); kwargs...)
```

### Keyword arguments

* `v̄ = (t,q,v) -> nothing`
* `f̄ = f`
* `invariants = nothing`
* `parameters = nothing`
* `periodicity = nothing`

"""
struct IODE{ϑType <: Callable, fType <: Callable, gType <: Callable,
            v̄Type <: Callable, f̄Type <: Callable,
            invType <: OptionalInvariants,
            parType <: OptionalParameters,
            perType <: OptionalPeriodicity} <: AbstractEquationPODE{invType,parType,perType}

    ϑ::ϑType
    f::fType
    g::gType
    v̄::v̄Type
    f̄::f̄Type

    invariants::invType
    parameters::parType
    periodicity::perType

    function IODE(ϑ, f, g, v̄, f̄, invariants, parameters, periodicity)
        new{typeof(ϑ), typeof(f), typeof(g), typeof(v̄), typeof(f̄),
            typeof(invariants), typeof(parameters), typeof(periodicity)}(
                ϑ, f, g, v̄, f̄, invariants, parameters, periodicity)
    end
end

_iode_default_v̄(t,q,v,params) = nothing

IODE(ϑ, f, g; invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity(), v̄=_iode_default_v̄, f̄=f) = IODE(ϑ, f, g, v̄, f̄, invariants, parameters, periodicity)

GeometricBase.invariants(equation::IODE) = equation.invariants
GeometricBase.parameters(equation::IODE) = equation.parameters
GeometricBase.periodicity(equation::IODE) = equation.periodicity

hasvectorfield(::IODE) = true

function check_initial_conditions(::IODE, ics::NamedTuple)
    haskey(ics, :q) || return false
    haskey(ics, :p) || return false
    haskey(ics, :λ) || return false
    eltype(ics.q) == eltype(ics.p) == eltype(ics.λ) || return false
    typeof(ics.q) == typeof(ics.p) == typeof(ics.λ) || return false
    axes(ics.q) == axes(ics.p) == axes(ics.λ) || return false
    return true
end

function check_methods(equ::IODE, tspan, ics::NamedTuple, params)
    applicable(equ.ϑ, tspan[begin], ics.q, zero(ics.q), zero(ics.p), params) || return false
    applicable(equ.f, tspan[begin], ics.q, zero(ics.q), zero(ics.p), params) || return false
    applicable(equ.g, tspan[begin], ics.q, zero(ics.q), zero(ics.p), params) || return false
    return true
end

function GeometricBase.datatype(equ::IODE, ics::NamedTuple)
    @assert check_initial_conditions(equ, ics)
    return eltype(ics.q)
end

function GeometricBase.arrtype(equ::IODE, ics::NamedTuple)
    @assert check_initial_conditions(equ, ics)
    typeof(ics.q)
end

_get_ϑ(equ::IODE, params) = (t,q,v,ϑ) -> equ.ϑ(t, q, v, ϑ, params)
_get_f(equ::IODE, params) = (t,q,v,f) -> equ.f(t, q, v, f, params)
_get_g(equ::IODE, params) = (t,q,v,g) -> equ.g(t, q, v, g, params)
_get_v̄(equ::IODE, params) = (t,q,v) -> equ.v̄(t, q, v, params)
_get_f̄(equ::IODE, params) = (t,q,v,f) -> equ.f̄(t, q, v, f, params)
_get_invariant(::IODE, inv, params) = (t,q,v) -> inv(t, q, v, params)

_functions(equ::IODE) = (ϑ = equ.ϑ, f = equ.f, g = equ.g, v̄ = equ.v̄, f̄ = equ.f̄)
_functions(equ::IODE, params::OptionalParameters) = (ϑ = _get_ϑ(equ, params),
                                                     f = _get_f(equ, params),
                                                     g = _get_g(equ, params),
                                                     v̄ = _get_v̄(equ, params),
                                                     f̄ = _get_f̄(equ, params))
