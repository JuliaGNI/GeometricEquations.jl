@doc raw"""
`DAE`: Differential Algebraic Equation

Defines a differential algebraic initial value problem
```math
\begin{aligned}
\dot{q} (t) &= v(t, q(t)) + u(t, q(t), \lambda(t)) , & q(t_{0}) &= q_{0} , \\
0 &= \phi (t, q(t), \lambda(t)) , & \lambda(t_{0}) &= \lambda_{0} ,
\end{aligned}
```
with vector field ``v``, projection ``u``, algebraic constraint ``\phi=0``,
initial conditions ``q_{0}`` and ``\lambda_{0}``, the dynamical variable ``q``
taking values in ``\mathbb{R}^{d}`` and the algebraic variable ``\lambda``
taking values in ``\mathbb{R}^{m}``.

### Parameters

* `DT <: Number`: data type
* `TT <: Real`: time step type
* `AT <: AbstractArray{DT}`: array type
* `vType <: Function`: type of `v`
* `uType <: Function`: type of `u`
* `ūType <: OptionalFunction`: type of `ū`
* `ϕType <: Function`: type of `ϕ`
* `ψType <: OptionalFunction`: type of `ψ`
* `v̄Type <: Function`: type of `v̄`
* `invType <: OptionalNamedTuple`: invariants type
* `parType <: OptionalNamedTuple`: parameters type
* `perType <: OptionalArray{AT}`: periodicity type

### Fields

* `d`: dimension of dynamical variable ``q`` and the vector field ``v``
* `m`: dimension of algebraic variable ``\lambda`` and the constraint ``\phi``
* `v`: function computing the vector field
* `u`: function computing the projection
* `ū`: function computing the secondary projection field ``\bar{u}`` (optional)
* `ϕ`: algebraic constraint
* `ψ`: secondary constraints (optional)
* `v̄`: function computing an initial guess for the velocity field ``v`` (defaults to `v`)
* `t₀`: initial time
* `q₀`: initial condition for dynamical variable ``q``
* `λ₀`: initial condition for algebraic variable ``\lambda``
* `μ₀`: initial condition for algebraic variable ``μ`` (optional)
* `invariants`: either a `NamedTuple` containing the equation's invariants or `nothing`
* `parameters`: either a `NamedTuple` containing the equation's parameters or `nothing`
* `periodicity`: determines the periodicity of the state vector `q` for cutting periodic solutions

The function `v`, providing the vector field, takes three arguments,
`v(t, q, v)`, the functions `u` and `ϕ`, providing the projection and the
algebraic constraint take four arguments, `u(t, q, λ, u)` and `ϕ(t, q, λ, ϕ)`,
where `t` is the current time, `q` and `λ` are the current solution vectors,
and `v`, `u` and `ϕ` are the vectors which hold the result of evaluating the
vector field ``v``, the projection ``u`` and the algebraic constraint ``\phi``
on `t`, `q` and `λ`.

### Constructors

```julia
DAE(v, u, ū, ϕ, ψ, v̄, t₀, q₀, λ₀, invariants, parameters, periodicity)

DAE(v, u, ϕ, t₀, q₀::StateVector, λ₀::StateVector, μ₀::StateVector=zero(λ₀); kwargs...)
DAE(v, u, ϕ, q₀::StateVector, λ₀::StateVector, μ₀::StateVector=zero(λ₀); kwargs...)
DAE(v, u, ϕ, t₀, q₀::State, λ₀::State, μ₀::State=zero(λ₀); kwargs...)
DAE(v, u, ϕ, q₀::State, λ₀::State, μ₀::State=zero(λ₀); kwargs...)

DAE(v, u, ū, ϕ, ψ, t₀, q₀::StateVector, λ₀::StateVector, μ₀::StateVector=zero(λ₀); kwargs...)
DAE(v, u, ū, ϕ, ψ, q₀::StateVector, λ₀::StateVector, μ₀::StateVector=zero(λ₀); kwargs...)
DAE(v, u, ū, ϕ, ψ, t₀, q₀::State, λ₀::State, μ₀::State=zero(λ₀); kwargs...)
DAE(v, u, ū, ϕ, ψ, q₀::State, λ₀::State, μ₀::State=zero(λ₀); kwargs...)
```

### Example

```julia
    function v(t, q, v)
        v[1] = q[1]
        v[2] = q[2]
    end

    function u(t, q, λ, u)
        u[1] = +λ[1]
        u[2] = -λ[1]
    end

    function ϕ(t, q, λ, ϕ)
        ϕ[1] = q[2] - q[1]
    end

    t₀ = 0.
    q₀ = [1., 1.]
    λ₀ = [0.]

    dae = DAE(v, u, ϕ, t₀, q₀, λ₀)

```

"""
struct DAE{vType <: Callable,
           uType <: Callable, ϕType <: Callable,
           ūType <: OptionalCallable, ψType <: OptionalCallable,
           v̄Type <: Callable,
           invType <: OptionalInvariants,
           parType <: OptionalParameters,
           perType <: OptionalPeriodicity} <: AbstractEquationDAE

    v::vType
    u::uType
    ϕ::ϕType
    ū::ūType
    ψ::ψType
    v̄::v̄Type

    invariants::invType
    parameters::parType
    periodicity::perType

    function DAE(v, u, ϕ, ū, ψ, v̄, invariants, parameters, periodicity)
        @assert !isempty(methods(v))
        @assert !isempty(methods(u))
        @assert !isempty(methods(ϕ))
        @assert !isempty(methods(ū)) || ū === nothing
        @assert !isempty(methods(ψ)) || ψ === nothing
        @assert !isempty(methods(v̄))

        new{typeof(v), typeof(u), typeof(ϕ), typeof(ū), typeof(ψ), typeof(v̄),
            typeof(invariants), typeof(parameters), typeof(periodicity)}(
                v, u, ϕ, ū, ψ, v̄, invariants, parameters, periodicity)
    end
end



DAE(v, u, ϕ, ū, ψ; v̄=v, invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity()) = DAE(v, u, ϕ, ū, ψ, v̄, invariants, parameters, periodicity)
DAE(v, u, ϕ; kwargs...) = DAE(v, u, ϕ, nothing, nothing; kwargs...)

GeometricBase.invariants(equation::DAE) = equation.invariants
GeometricBase.parameters(equation::DAE) = equation.parameters
GeometricBase.periodicity(equation::DAE) = equation.periodicity

const DAEsecType{ΨT,VT,UT,ΦT,ŪT,V̄T,invT,parT,perT} = DAE{VT,UT,ΦT,ŪT,ΨT,V̄T,invT,parT,perT} # type alias for dispatch on secondary constraint type parameter
const DAEinvType{invT,VT,UT,ΦT,ŪT,ΨT,V̄T,parT,perT} = DAE{VT,UT,ΦT,ŪT,ΨT,V̄T,invT,parT,perT} # type alias for dispatch on invariants type parameter
const DAEparType{parT,VT,UT,ΦT,ŪT,ΨT,V̄T,invT,perT} = DAE{VT,UT,ΦT,ŪT,ΨT,V̄T,invT,parT,perT} # type alias for dispatch on parameters type parameter
const DAEperType{perT,VT,UT,ΦT,ŪT,ΨT,V̄T,invT,parT} = DAE{VT,UT,ΦT,ŪT,ΨT,V̄T,invT,parT,perT} # type alias for dispatch on periodicity type parameter

hasvectorfield(::DAE) = true

hassecondary(::DAEsecType{<:Nothing}) = false
hassecondary(::DAEsecType{<:Callable}) = true

hasinvariants(::DAEinvType{<:NullInvariants}) = false
hasinvariants(::DAEinvType{<:NamedTuple}) = true

hasparameters(::DAEparType{<:NullParameters}) = false
hasparameters(::DAEparType{<:NamedTuple}) = true

hasperiodicity(::DAEperType{<:NullPeriodicity}) = false
hasperiodicity(::DAEperType{<:AbstractArray}) = true

function check_initial_conditions(equ::DAE, ics::NamedTuple)
    haskey(ics, :q) || return false
    haskey(ics, :λ) || return false
    if hassecondary(equ)
        haskey(ics, :μ) || return false
        eltype(ics.λ) == eltype(ics.μ) || return false
        typeof(ics.λ) == typeof(ics.μ) || return false
        axes(ics.λ) == axes(ics.μ) || return false
    end
    return true
end

function check_methods(equ::DAE, tspan, ics::NamedTuple, params)
    applicable(equ.v, tspan[begin], ics.q, zero(ics.q), params) || return false
    applicable(equ.u, tspan[begin], ics.q, ics.λ, zero(ics.q), params) || return false
    applicable(equ.ϕ, tspan[begin], ics.q, zero(ics.λ), params) || return false
    return true
end

function datatype(equ::DAE, ics::NamedTuple)
    @assert check_initial_conditions(equ, ics)
    return eltype(ics.q)
end

function arrtype(equ::DAE, ics::NamedTuple)
    @assert check_initial_conditions(equ, ics)
    return typeof(ics.q)
end

_get_v(equ::DAE, params) = (t,q,v)   -> equ.v(t, q, v, params)
_get_u(equ::DAE, params) = (t,q,λ,u) -> equ.u(t, q, λ, u, params)
_get_ϕ(equ::DAE, params) = (t,q,ϕ)   -> equ.ϕ(t, q, ϕ, params)
_get_ū(equ::DAE, params) = (t,q,λ,u) -> equ.ū(t, q, λ, u, params)
_get_ψ(equ::DAE, params) = (t,q,v,ψ) -> equ.ψ(t, q, v, ϕ, params)
_get_v̄(equ::DAE, params) = (t,q,v)   -> equ.v̄(t, q, v, params)

_functions(equ::DAEsecType{<:Nothing}) = (v = equ.v, u = equ.u, ϕ = equ.ϕ, v̄ = equ.v̄)
_functions(equ::DAEsecType{<:Callable}) = (v = equ.v, u = equ.u, ϕ = equ.ϕ, ū = equ.ū, ψ = equ.ψ, v̄ = equ.v̄)
_functions(equ::DAEsecType{<:Nothing}, params::OptionalParameters) = (v = _get_v(equ, params), u = _get_u(equ, params), ϕ = _get_ϕ(equ, params), v̄ = _get_v̄(equ, params))
_functions(equ::DAEsecType{<:Callable}, params::OptionalParameters) = (v = _get_v(equ, params), u = _get_u(equ, params), ϕ = _get_ϕ(equ, params), ū = _get_ū(equ, params), ψ = _get_ψ(equ, params), v̄ = _get_v̄(equ, params))
