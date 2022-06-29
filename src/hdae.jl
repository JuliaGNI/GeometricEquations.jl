@doc raw"""
`HDAE`: Hamiltonian Differential Algebraic Equation

Defines a Hamiltonian differential algebraic initial value problem, that is
a canonical Hamiltonian system of equations subject to Dirac constraints,
```math
\begin{aligned}
\dot{q} (t) &= v(t, q(t), p(t)) + u(t, q(t), p(t), \lambda(t)) + \bar{g}(t, q(t), p(t), \lambda(t), \gamma(t)) , & q(t_{0}) &= q_{0} , \\
\dot{p} (t) &= f(t, q(t), p(t)) + g(t, q(t), p(t), \lambda(t)) + \bar{f}(t, q(t), p(t), \lambda(t), \gamma(t)) , & p(t_{0}) &= p_{0} , \\
0 &= \phi (t, q(t), p(t)) , \\
0 &= \psi (t, q(t), p(t), \dot{q}(t), \dot{p}(t)) ,
\end{aligned}
```
with vector fields ``v``, ``u``, ``\bar{u}`` and ``f``, ``g``, ``\bar{g}``,
primary constraint ``\phi(q,p)=0`` and secondary constraint ``\psi(q,p,\lambda)=0``,
initial conditions ``(q_{0}, p_{0})``, the dynamical variables ``(q,p)``
taking values in ``\mathbb{R}^{d} \times \mathbb{R}^{d}`` and
the algebraic variables ``(\lambda, \gamma)`` taking values in
``\mathbb{R}^{m} \times \mathbb{R}^{m}``.

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
* `PType <: Function`: type of `P`
* `v̄Type <: Function`: type of `v̄`
* `f̄Type <: Function`: type of `f̄`
* `hamType <: Function`: Hamiltonian type
* `invType <: OptionalNamedTuple`: invariants type
* `parType <: OptionalNamedTuple`: parameters type
* `perType <: OptionalArray{AT}`: periodicity type

### Fields

* `d`: dimension of dynamical variables ``q`` and ``p`` as well as the vector fields ``v`` and ``f``
* `m`: dimension of algebraic variables ``\lambda`` and ``\gamma`` and the constraints ``\phi`` and ``\psi``
* `v`: function computing the Hamiltonian vector field ``v``
* `f`: function computing the Hamiltonian vector field ``f``
* `u`: function computing the projection for ``q``
* `g`: function computing the primary projection field ``g``
* `ϕ`: primary constraints
* `ū`: function computing the secondary projection field ``\bar{u}`` (optional)
* `ḡ`: function computing the secondary projection field ``\bar{g}`` (optional)
* `ψ`: secondary constraints (optional)
* `P`: function computing the Poisson matrix ``P``
* `v̄`: function computing an initial guess for the velocity field ``v``` (optional, defaults to `v`)
* `f̄`: function computing an initial guess for the force field ``f`` (optional, defaults to `f`)
* `t₀`: initial time (optional)
* `q₀`: initial condition for dynamical variable ``q``
* `p₀`: initial condition for dynamical variable ``p``
* `λ₀`: initial condition for algebraic variable ``λ``
* `μ₀`: initial condition for algebraic variable ``μ`` (optional)
* `hamiltonian`: function computing the Hamiltonian ``H``
* `invariants`: either a `NamedTuple` containing the equation's invariants or `nothing`
* `parameters`: either a `NamedTuple` containing the equation's parameters or `nothing`
* `periodicity`: determines the periodicity of the state vector `q` for cutting periodic solutions

### Constructors

```julia
HDAE(v, f, u, g, ϕ, ū, ḡ, ψ, v̄, f̄, P, t₀, q₀, p₀, λ₀, hamiltonian, invariants, parameters, periodicity)

HDAE(v, f, u, g, ϕ, h, t₀, q₀::StateVector, p₀::StateVector, λ₀::StateVector; kwargs...)
HDAE(v, f, u, g, ϕ, h, q₀::StateVector, p₀::StateVector, λ₀::StateVector; kwargs...)
HDAE(v, f, u, g, ϕ, h, t₀, q₀::State, p₀::State, λ₀::State; kwargs...)
HDAE(v, f, u, g, ϕ, h, q₀::State, p₀::State, λ₀::State; kwargs...)

HDAE(v, f, u, g, ϕ, ū, ḡ, ψ, h, t₀, q₀::StateVector, p₀::StateVector, λ₀::StateVector; kwargs...)
HDAE(v, f, u, g, ϕ, ū, ḡ, ψ, h, q₀::StateVector, p₀::StateVector, λ₀::StateVector; kwargs...)
HDAE(v, f, u, g, ϕ, ū, ḡ, ψ, h, t₀, q₀::State, p₀::State, λ₀::State; kwargs...)
HDAE(v, f, u, g, ϕ, ū, ḡ, ψ, h, q₀::State, p₀::State, λ₀::State; kwargs...)
```

"""
struct HDAE{vType <: Callable, fType <: Callable,
            uType <: Callable, gType <: Callable, ϕType <: Callable,
            ūType <: OptionalCallable, ḡType <: OptionalCallable, ψType <: OptionalCallable,
            v̄Type <: Callable, f̄Type <: Callable,
            poiType <: Callable,
            hamType <: Callable,
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

    poisson::poiType
    hamiltonian::hamType
    invariants::invType
    parameters::parType
    periodicity::perType

    function HDAE(v, f, u, g, ϕ, ū, ḡ, ψ, v̄, f̄, poisson, hamiltonian, invariants, parameters, periodicity)
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
        @assert !isempty(methods(poisson))
        @assert !isempty(methods(hamiltonian))

        new{typeof(v), typeof(f),
            typeof(u), typeof(g), typeof(ϕ),
            typeof(ū), typeof(ḡ), typeof(ψ),
            typeof(v̄), typeof(f̄),
            typeof(poisson), typeof(hamiltonian), typeof(invariants), typeof(parameters), typeof(periodicity)}(
                v, f, u, g, ϕ, ū, ḡ, ψ, v̄, f̄, poisson, hamiltonian, invariants, parameters, periodicity)
    end
end

HDAE(v, f, u, g, ϕ, ū, ḡ, ψ, poisson, hamiltonian; v̄=v, f̄=f, invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity()) = HDAE(v, f, u, g, ϕ, ū, ḡ, ψ, v̄, f̄, poisson, hamiltonian, invariants, parameters, periodicity)
HDAE(v, f, u, g, ϕ, poisson, hamiltonian; kwargs...) = HDAE(v, f, u, g, ϕ, nothing, nothing, nothing, poisson, hamiltonian; kwargs...)

GeometricBase.invariants(equation::HDAE) = equation.invariants
GeometricBase.parameters(equation::HDAE) = equation.parameters
GeometricBase.periodicity(equation::HDAE) = equation.periodicity

hasvectorfield(::HDAE) = true
hashamiltonian(::HDAE) = true

function check_initial_conditions(equ::HDAE, ics::NamedTuple)
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

function check_methods(equ::HDAE, tspan, ics::NamedTuple, params)
    applicable(equ.v, tspan[begin], ics.q, ics.p, zero(ics.q), params) || return false
    applicable(equ.f, tspan[begin], ics.q, ics.p, zero(ics.p), params) || return false
    applicable(equ.ϕ, tspan[begin], ics.q, ics.p, zero(ics.λ), params) || return false
    # TODO add missing methods
    return true
end

function datatype(equ::HDAE, ics::NamedTuple)
    @assert check_initial_conditions(equ, ics)
    return eltype(ics.q)
end

function arrtype(equ::HDAE, ics::NamedTuple)
    @assert check_initial_conditions(equ, ics)
    return typeof(ics.q)
end

_get_v(equ::HDAE, params) = (t,q,p,v)     -> equ.v(t, q, p, v, params)
_get_f(equ::HDAE, params) = (t,q,p,f)     -> equ.f(t, q, p, f, params)
_get_u(equ::HDAE, params) = (t,q,p,λ,u)   -> equ.u(t, q, p, λ, u, params)
_get_g(equ::HDAE, params) = (t,q,p,λ,g)   -> equ.g(t, q, p, λ, g, params)
_get_ϕ(equ::HDAE, params) = (t,q,p,ϕ)     -> equ.ϕ(t, q, p, ϕ, params)
_get_ū(equ::HDAE, params) = (t,q,p,λ,u)   -> equ.ū(t, q, p, λ, u, params)
_get_ḡ(equ::HDAE, params) = (t,q,p,λ,g)   -> equ.ḡ(t, q, p, λ, g, params)
_get_ψ(equ::HDAE, params) = (t,q,p,v,f,ψ) -> equ.ψ(t, q, p, v, f, ψ, params)
_get_v̄(equ::HDAE, params) = (t,q,p,v)     -> equ.v̄(t, q, p, v, params)
_get_f̄(equ::HDAE, params) = (t,q,p,f)     -> equ.f̄(t, q, p, f, params)
_get_h(equ::HDAE, params) = (t,q,p) -> equ.hamiltonian(t, q, p, params)
_get_poisson(equ::HDAE, params) = (t,q,p,ω) -> equ.poisson(t, q, p, ω, params)
_get_invariant(::HDAE, inv, params) = (t,q,p) -> inv(t, q, p, params)

function _functions(equ::HDAE)
    if hassecondary(equ)
        (v = equ.v, f = equ.f, u = equ.u, g = equ.g, ϕ = equ.ϕ, ū = equ.ū, ḡ = equ.ḡ, ψ = equ.ψ, v̄ = equ.v̄, f̄ = equ.f̄, poisson = equ.poisson, h = equ.hamiltonian)
    else
        (v = equ.v, f = equ.f, u = equ.u, g = equ.g, ϕ = equ.ϕ, v̄ = equ.v̄, f̄ = equ.f̄, poisson = equ.poisson, h = equ.hamiltonian)
    end
end

function _functions(equ::HDAE, params::OptionalParameters)
    if hassecondary(equ)
        (v = _get_v(equ, params), f = _get_f(equ, params), u = _get_u(equ, params), g = _get_g(equ, params), ϕ = _get_ϕ(equ, params), ū = _get_ū(equ, params), ḡ = _get_ḡ(equ, params), ψ = _get_ψ(equ, params), v̄ = _get_v̄(equ, params), f̄ = _get_f̄(equ, params), poisson = _get_poisson(equ, params), h = _get_h(equ, params))
    else
        (v = _get_v(equ, params), f = _get_f(equ, params), u = _get_u(equ, params), g = _get_g(equ, params), ϕ = _get_ϕ(equ, params), v̄ = _get_v̄(equ, params), f̄ = _get_f̄(equ, params), poisson = _get_poisson(equ, params), h = _get_h(equ, params))
    end
end
