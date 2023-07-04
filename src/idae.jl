
const idae_equations = raw"""
An implicit differential algebraic initial value problem takes the form
```math
\begin{aligned}
\dot{q} (t) &= v(t) + u(t, q(t), v(t), p(t), \lambda(t)) , & q(t_{0}) &= q_{0} , \\
\dot{p} (t) &= f(t, q(t), v(t)) + g(t, q(t), v(t), p(t), \lambda(t)) , & p(t_{0}) &= p_{0} , \\
p(t) &= ϑ(t, q(t), v(t)) , && \\
0 &= \phi (t, q(t), v(t), p(t)) , & \lambda(t_{0}) &= \lambda_{0} ,
\end{aligned}
```
with force field ``f``, the momentum defined by ``ϑ``, projections ``u`` and ``g``,
algebraic constraint ``\phi(t,q,v,p)=0``, the dynamical variables ``(q,p)`` taking
values in ``\mathbb{R}^{d} \times \mathbb{R}^{d}``, the algebraic variable ``\lambda``
taking values in ``\mathbb{R}^{m}``, and initial conditions ``(q_{0}, p_{0})`` and
``\lambda_{0}``.

Some integrators also enforce the secondary constraint ``\psi``, that is the time
derivative of the algebraic constraint ``\phi``.
In this case, the system of equations is modified as follows
```math
\begin{aligned}
\dot{q} (t) &= v(t) + u(t, q(t), v(t), p(t), \lambda(t)) + \bar{u} (t, q(t), v(t), p(t), \gamma(t)) , & q(t_{0}) &= q_{0} , \\
\dot{p} (t) &= f(t, q(t), v(t)) + g(t, q(t), v(t), p(t), \lambda(t)) + \bar{g} (t, q(t), v(t), p(t), \gamma(t)) , & p(t_{0}) &= p_{0} , \\
p(t) &= ϑ(t, q(t), v(t)) , && \\
0 &= \phi (t, q(t), v(t), p(t)) , & \lambda(t_{0}) &= \lambda_{0} , \\
0 &= \psi (t, q(t), v(t), p(t), \dot{q} (t), \dot{p} (t)) , & \gamma(t_{0}) &= \gamma_{0} ,
\end{aligned}
```
with the second algebraic variable ``\gamma`` also taking values in ``\mathbb{R}^{m}``.
"""

const idae_constructors = raw"""
The function `ϑ` computes the momentum, `f` computes the force field, `u` and `g` compute
the projections, and `ϕ` provides the algebraic constraint.
The functions `ψ`, `ū` and `ḡ` are optional and provide the secondary constraint, that is
the time derivative of the algebraic constraint, and the corresponding projection.
"""

const idae_functions = raw"""
The functions `ϑ` and `f` must have the interface
```julia
function ϑ(p, t, q, v)
    p[1] = ...
    p[2] = ...
    ...
end

function f(f, t, q, v)
    f[1] = ...
    f[2] = ...
    ...
end
```
where `t` is the current time, `q` is the current solution vector, `v` is the
current velocity and `f` and `p` are the vectors which hold the result of
evaluating the functions ``f`` and ``ϑ`` on `t`, `q` and `v`.
The funtions `g`, `v̄` and `f̄` are specified by
```julia
function u(u, t, q, v, p, λ)
    u[1] = ...
    u[2] = ...
    ...
end

function g(g, t, q, v, p, λ)
    g[1] = ...
    g[2] = ...
    ...
end

function v̄(v, t, q, p)
    v[1] = ...
    v[2] = ...
    ...
end

function f̄(f, t, q, v)
    f[1] = ...
    f[2] = ...
    ...
end
```

Some integrators also enforce the secondary constraint ``\psi`` and require
the following additional functions
```
function ū(u, t, q, v, p, γ, params)
    u[1] = ...
    u[2] = ...
    ...
end

function ḡ(g, t, q, v, p, γ, params)
    g[1] = ...
    g[2] = ...
    ...
end

function ψ(ψ, t, q, v, p, q̇, ṗ, params)
    ψ[1] = ...
end
```
"""


@doc """
`IDAE`: Implicit Differential Algebraic Equation

$(idae_equations)

### Parameters

* `ϑType <: Callable`: type of `ϑ`
* `fType <: Callable`: type of `f`
* `uType <: Callable`: type of `u`
* `gType <: Callable`: type of `g`
* `ϕType <: Callable`: type of `ϕ`
* `ūType <: Callable`: type of `ū`
* `ḡType <: Callable`: type of `ḡ`
* `ψType <: Callable`: type of `ψ`
* `v̄Type <: Callable`: type of `v̄`
* `f̄Type <: Callable`: type of `f̄`
* `invType <: OptionalInvariants`: invariants type
* `parType <: OptionalParameters`: parameters type
* `perType <: OptionalPeriodicity`: periodicity type

### Fields

* `ϑ`: function determining the momentum
* `f`: function computing the vector field ``f``
* `u`: function computing the projection for ``q``
* `g`: function computing the projection for ``p``
* `ϕ`: algebraic constraints
* `ū`: function computing the secondary projection field ``\\bar{u}`` (*optional*)
* `ḡ`: function computing the secondary projection field ``\\bar{g}`` (*optional*)
* `ψ`: secondary constraints (*optional*)
* `v̄`: function computing an initial guess for the velocity field ``v`` (*optional*)
* `f̄`: function computing an initial guess for the force field ``f`` (*optional*)
* `invariants`: functions for the computation of invariants, either a `NamedTuple` containing the equation's invariants or `NullInvariants`
* `parameters`: type constraints for parameters, either a `NamedTuple` containing the equation's parameters or `NullParameters`
* `periodicity`: determines the periodicity of the state vector `q` for cutting periodic solutions, either a `AbstractArray` or `NullPeriodicity`

### Constructors

```julia
IDAE(ϑ, f, u, g, ϕ, ū, ḡ, ψ, v̄, f̄, invariants, parameters, periodicity)
IDAE(ϑ, f, u, g, ϕ, ū, ḡ, ψ; v̄ = _idae_default_v̄, f̄ = f, invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
IDAE(ϑ, f, u, g, ϕ; v̄ = _idae_default_v̄, f̄ = f, invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
```

where 

```julia
_idae_default_v̄(v, t, q, params) = nothing
```

$(idae_constructors)

### Function Definitions

$(idae_functions)

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

_idae_default_v̄(t, q, v, p, params) = nothing

IDAE(ϑ, f, u, g, ϕ, ū, ḡ, ψ, v̄, f̄; invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity()) = IDAE(ϑ, f, u, g, ϕ, ū, ḡ, ψ, v̄, f̄, invariants, parameters, periodicity)
IDAE(ϑ, f, u, g, ϕ, ū, ḡ, ψ; v̄=_idae_default_v̄, f̄=f, kwargs...) = IDAE(ϑ, f, u, g, ϕ, ū, ḡ, ψ, v̄, f̄; kwargs...)
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
    applicable(equ.u, zero(ics.q), tspan[begin], ics.q, vectorfield(ics.q), ics.p, ics.λ, params) || return false
    applicable(equ.g, zero(ics.p), tspan[begin], ics.q, vectorfield(ics.q), ics.p, ics.λ, params) || return false
    applicable(equ.ϕ, zero(ics.λ), tspan[begin], ics.q, vectorfield(ics.q), ics.p, params) || return false
    applicable(equ.v̄, zero(ics.q), tspan[begin], ics.q, ics.p, params) || return false
    applicable(equ.f̄, zero(ics.p), tspan[begin], ics.q, vectorfield(ics.q), params) || return false
    equ.ū === nothing || applicable(equ.ū, zero(ics.q), tspan[begin], ics.q, vectorfield(ics.q), ics.p, ics.λ, params) || return false
    equ.ḡ === nothing || applicable(equ.ḡ, zero(ics.p), tspan[begin], ics.q, vectorfield(ics.q), ics.p, ics.λ, params) || return false
    equ.ψ === nothing || applicable(equ.ψ, zero(ics.λ), tspan[begin], ics.q, vectorfield(ics.q), ics.p, vectorfield(ics.q), vectorfield(ics.p), params) || return false
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
_get_u(equ::IDAE, params) = (u, t, q, v, p, λ) -> equ.u(u, t, q, v, p, λ, params)
_get_g(equ::IDAE, params) = (g, t, q, v, p, λ) -> equ.g(g, t, q, v, p, λ, params)
_get_ϕ(equ::IDAE, params) = (ϕ, t, q, v, p)    -> equ.ϕ(ϕ, t, q, v, p, params)
_get_ū(equ::IDAE, params) = (u, t, q, v, p, λ)    -> equ.ū(u, t, q, v, p, λ, params)
_get_ḡ(equ::IDAE, params) = (g, t, q, v, p, λ)    -> equ.ḡ(g, t, q, v, p, λ, params)
_get_ψ(equ::IDAE, params) = (ψ, t, q, v, p, q̇, ṗ) -> equ.ψ(ψ, t, q, v, p, q̇, ṗ, params)
_get_v̄(equ::IDAE, params) = (v, t, q, p)       -> equ.v̄(v, t, q, p, params)
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
        (
            ϑ = _get_ϑ(equ, params),
            f = _get_f(equ, params),
            u = _get_u(equ, params),
            g = _get_g(equ, params),
            ϕ = _get_ϕ(equ, params),
            ū = _get_ū(equ, params),
            ḡ = _get_ḡ(equ, params),
            ψ = _get_ψ(equ, params),
            v̄ = _get_v̄(equ, params),
            f̄ = _get_f̄(equ, params)
        )
    else
        (
            ϑ = _get_ϑ(equ, params),
            f = _get_f(equ, params),
            u = _get_u(equ, params),
            g = _get_g(equ, params),
            ϕ = _get_ϕ(equ, params),
            v̄ = _get_v̄(equ, params),
            f̄ = _get_f̄(equ, params)
        )
    end
end


@doc """
`IDAEProblem`: Implicit Differential Algebraic Equation Problem

$(idae_equations)

### Constructors

```julia
IDAEProblem(ϑ, f, u, g, ϕ, ū, ḡ, ψ, tspan, tstep, ics; kwargs...)
IDAEProblem(ϑ, f, u, g, ϕ, ū, ḡ, ψ, tspan, tstep, q₀::State, p₀::State, λ₀::State = zero(q₀); kwargs...)
IDAEProblem(ϑ, f, u, g, ϕ, tspan, tstep, ics; kwargs...)
IDAEProblem(ϑ, f, u, g, ϕ, tspan, tstep, q₀::State, p₀::State, λ₀::State = zero(q₀); kwargs...)
```

$(idae_constructors)

`tspan` is the time interval `(t₀,t₁)` for the problem to be solved in,
`tstep` is the time step to be used in the simulation, and
`ics` is a `NamedTuple` with entries `q` and `p`.
The initial conditions `q₀` and `p₀` can also be prescribed
directly, with `State` an `AbstractArray{<:Number}`.

In addition to the standard keyword arguments for [`GeometricProblem`](@ref GeometricEquations.GeometricProblem) subtypes,
an `IDAEProblem` accepts functions `v̄` and `f̄` for the computation of initial guesses for the vector fields with default
values `v̄ = _idae_default_v̄` and `f̄ = f`.

### Function Definitions

$(idae_functions)

"""
const IDAEProblem = GeometricProblem{IDAE}

function IDAEProblem(ϑ, f, u, g, ϕ, ū, ḡ, ψ, tspan, tstep, ics::NamedTuple;
                     v̄ = _idae_default_v̄, f̄ = f, invariants = NullInvariants(),
                     parameters = NullParameters(), periodicity = NullPeriodicity())
    equ = IDAE(ϑ, f, u, g, ϕ, ū, ḡ, ψ, v̄, f̄, invariants, parameter_types(parameters),
               periodicity)
    GeometricProblem(equ, tspan, tstep, ics, parameters)
end

function IDAEProblem(ϑ, f, u, g, ϕ, ū, ḡ, ψ, tspan, tstep, q₀::State, p₀::State,
                     λ₀::State, μ₀::State = zero(λ₀); kwargs...)
    ics = (q = q₀, p = p₀, λ = λ₀, μ = μ₀)
    IDAEProblem(ϑ, f, u, g, ϕ, ū, ḡ, ψ, tspan, tstep, ics; kwargs...)
end

function IDAEProblem(ϑ, f, u, g, ϕ, tspan, tstep, ics::NamedTuple; kwargs...)
    IDAEProblem(ϑ, f, u, g, ϕ, nothing, nothing, nothing, tspan, tstep, ics; kwargs...)
end

function IDAEProblem(ϑ, f, u, g, ϕ, tspan, tstep, q₀::State, p₀::State, λ₀::State;
                     kwargs...)
    ics = (q = q₀, p = p₀, λ = λ₀)
    IDAEProblem(ϑ, f, u, g, ϕ, tspan, tstep, ics; kwargs...)
end

function GeometricBase.periodicity(prob::IDAEProblem)
    (q = periodicity(equation(prob)), p = NullPeriodicity(), λ = NullPeriodicity(),
     μ = NullPeriodicity())
end


const IDAEEnsemble  = GeometricEnsemble{IDAE}

function IDAEEnsemble(ϑ, f, u, g, ϕ, ū, ḡ, ψ, tspan, tstep, ics::AbstractVector{<:NamedTuple}; v̄ = _idae_default_v̄, f̄ = f,
        invariants = NullInvariants(),
        parameters = NullParameters(),
        periodicity = NullPeriodicity())
    equ = IDAE(ϑ, f, u, g, ϕ, ū, ḡ, ψ, v̄, f̄, invariants, parameter_types(parameters), periodicity)
    GeometricEnsemble(equ, tspan, tstep, ics, parameters)
end

function IDAEEnsemble(ϑ, f, u, g, ϕ, tspan, tstep, ics::AbstractVector{<:NamedTuple}; kwargs...)
    IDAEEnsemble(ϑ, f, u, g, ϕ, nothing, nothing, nothing, tspan, tstep, ics; kwargs...)
end
