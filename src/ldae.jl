
const ldae_equations = raw"""
A special case of an implicit initial value problem is a Lagrangian differential
algebraic equation of the form
```math
\begin{aligned}
\dot{q} (t) &= v(t) + u(t, q(t), v(t), p(t), \lambda(t)) + \bar{u} (t, q(t), v(t), p(t), \mu(t)) , \\
\dot{p} (t) &= f(t, q(t), v(t)) + g(t, q(t), v(t), p(t), \lambda(t)) + \bar{g} (t, q(t), v(t), p(t), \mu(t)) , \\
p(t) &= ϑ(t, q(t), v(t)) , \\
0 &= \phi (t, q(t), v(t), p(t)) , \\
0 &= \psi (t, q(t), v(t), p(t), \dot{q}(t), \dot{p}(t)) ,
\end{aligned}
```
with momentum ``p`` and force field ``f``, given by
```math
\begin{aligned}
p &= \frac{\partial L}{\partial v} (q,v) , &
f &= \frac{\partial L}{\partial q} (q,v) ,
\end{aligned}
```
projection fields ``u``, ``\bar{u}`` and ``g``, ``\bar{g}``.
This is a special case of a differential algebraic equation with dynamical
variables ``(q,p)`` and algebraic variables ``v``, ``\lambda`` and ``\mu``.
"""

const ldae_constructors = raw"""
The function `ϑ` computes the momentum, `f` computes the force field, `u` and `g` compute
the projections, and `ϕ` provides the algebraic constraint.
The functions `ψ`, `ū` and `ḡ` are optional and provide the secondary constraint, that is
the time derivative of the algebraic constraint, and the corresponding projection.
"""

const ldae_functions = raw"""
The functions `ϑ` and `f` must have the interface
```julia
function ϑ(p, t, q, v, params)
    p[1] = ...
    p[2] = ...
    ...
end

function f(f, t, q, v, params)
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
function u(u, t, q, v, p, μ, params)
    u[1] = ...
    u[2] = ...
    ...
end

function g(g, t, q, v, p, μ, params)
    g[1] = ...
    g[2] = ...
    ...
end

function v̄(v, t, q, p, params)
    v[1] = ...
    v[2] = ...
    ...
end

function f̄(f, t, q, v, params)
    f[1] = ...
    f[2] = ...
    ...
end
```
and the functions `ω` and `l`, computing the symplectic matrix and the Lagrangian,
have the following signature
```julia
function ω(f, t, q, v, params)
    ω[1,1] = ...
    ω[1,2] = ...
    ...
end

function l(t, q, v, params)
    return ...
end
```

Some integrators also enforce the secondary constraint ``\psi`` and require
the following additional functions
```
function ū(u, t, q, v, p, μ, params)
    u[1] = ...
    u[2] = ...
    ...
end

function ḡ(g, t, q, v, p, μ, params)
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
`LDAE`: Lagrangian Differential Algebraic Equation

$(ldae_equations)

### Parameters

* `ϑType <: Callable`: type of `ϑ`
* `fType <: Callable`: type of `f`
* `uType <: Callable`: type of `u`
* `gType <: Callable`: type of `g`
* `ϕType <: Callable`: type of `ϕ`
* `ūType <: Callable`: type of `ū`
* `ḡType <: Callable`: type of `ḡ`
* `ψType <: Callable`: type of `ψ`
* `ωType <: Callable`: type of `ω`
* `v̄Type <: Callable`: type of `v̄`
* `f̄Type <: Callable`: type of `f̄`
* `lagType <: Callable`: Lagrangian type
* `invType <: OptionalInvariants`: invariants type
* `parType <: OptionalParameters`: parameters type
* `perType <: OptionalPeriodicity`: periodicity type

### Fields

* `f`: function computing the vector field
* `u`: function computing the projection for ``q``, for a degenerate system given by ``\\lambda``
* `g`: function computing the projection for ``p``, for a degenerate system given by ``\\nabla \\vartheta (q) \\cdot \\lambda``
* `ϕ`: primary constraints, for a degenerate system given by ``p - \\vartheta (t,q)``
* `ū`: function computing the secondary projection field ``\\bar{u}``, for a degenerate system given by ``\\lambda`` (*optional*)
* `ḡ`: function computing the secondary projection field ``\\bar{g}``, for a degenerate system given by ``\\lambda \\cdot \\nabla \\vartheta (t,q)`` (*optional*)
* `ψ`: secondary constraints, for a degenerate system given by ``\\dot{p} - \\dot{q} \\cdot \\nabla \\vartheta (t,q)`` (*optional*)
* `ω`: function computing the symplectic matrix
* `v̄`: function computing an initial guess for the velocity field ``v`` (*optional*)
* `f̄`: function computing an initial guess for the force field ``f`` (*optional*)
* `lagrangian`: function computing the Lagrangian ``L``
* `invariants`: functions for the computation of invariants, either a `NamedTuple` containing the equation's invariants or `NullInvariants`
* `parameters`: type constraints for parameters, either a `NamedTuple` containing the equation's parameters or `NullParameters`
* `periodicity`: determines the periodicity of the state vector `q` for cutting periodic solutions, either a `AbstractArray` or `NullPeriodicity`


### Constructors

```julia
LDAE(ϑ, f, u, g, ϕ, ū, ḡ, ψ, ω, v̄, f̄, lagrangian, invariants, parameters, periodicity)
LDAE(ϑ, f, u, g, ϕ, ū, ḡ, ψ, ω, lagrangian; v̄ = _ldae_default_v̄, f̄ = f, invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
LDAE(ϑ, f, u, g, ϕ, ω, lagrangian; v̄ = _ldae_default_v̄, f̄ = f, invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
```

where 

```julia
_ldae_default_v̄(v, t, q, params) = nothing
```

$(ldae_constructors)

### Function Definitions

$(ldae_functions)

"""
struct LDAE{ϑType <: Callable, fType <: Callable,
            uType <: Callable, gType <: Callable, ϕType <: Callable,
            ūType <: OptionalCallable, ḡType <: OptionalCallable, ψType <: OptionalCallable,
            ωType <: Callable, v̄Type <: Callable, f̄Type <: Callable,
            lagType <: Callable,
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
    ω::ωType

    v̄::v̄Type
    f̄::f̄Type

    lagrangian::lagType
    invariants::invType
    parameters::parType
    periodicity::perType

    function LDAE(ϑ, f, u, g, ϕ, ū, ḡ, ψ, ω, v̄, f̄, lagrangian, invariants, parameters, periodicity)
        @assert !isempty(methods(ϑ))
        @assert !isempty(methods(f))
        @assert !isempty(methods(u))
        @assert !isempty(methods(g))
        @assert !isempty(methods(ϕ))
        @assert !isempty(methods(ū)) || ū === nothing
        @assert !isempty(methods(ḡ)) || ḡ === nothing
        @assert !isempty(methods(ψ)) || ψ === nothing
        @assert !isempty(methods(ω))
        @assert !isempty(methods(v̄))
        @assert !isempty(methods(f̄))
        @assert !isempty(methods(lagrangian))

        new{typeof(ϑ), typeof(f),
            typeof(u), typeof(g), typeof(ϕ),
            typeof(ū), typeof(ḡ), typeof(ψ),
            typeof(ω), typeof(v̄), typeof(f̄),
            typeof(lagrangian), typeof(invariants), typeof(parameters), typeof(periodicity)}(
                ϑ, f, u, g, ϕ, ū, ḡ, ψ, ω, v̄, f̄, lagrangian, invariants, parameters, periodicity)
    end
end

_ldae_default_v̄(t, q, v, p, params) = nothing

LDAE(ϑ, f, u, g, ϕ, ū, ḡ, ψ, ω, v̄, f̄, lagrangian; invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity()) = LDAE(ϑ, f, u, g, ϕ, ū, ḡ, ψ, ω, v̄, f̄, lagrangian, invariants, parameters, periodicity)
LDAE(ϑ, f, u, g, ϕ, ū, ḡ, ψ, ω, lagrangian; v̄=_ldae_default_v̄, f̄=f, kwargs...) = LDAE(ϑ, f, u, g, ϕ, ū, ḡ, ψ, ω, v̄, f̄, lagrangian; kwargs...)
LDAE(ϑ, f, u, g, ϕ, ω, lagrangian; kwargs...) = LDAE(ϑ, f, u, g, ϕ, nothing, nothing, nothing, ω, lagrangian; kwargs...)

GeometricBase.invariants(equation::LDAE) = equation.invariants
GeometricBase.parameters(equation::LDAE) = equation.parameters
GeometricBase.periodicity(equation::LDAE) = equation.periodicity

hasvectorfield(::LDAE) = true
haslagrangian(::LDAE) = true

function Base.show(io::IO, equation::LDAE)
    print(io, "Lagrangian Differential Algebraic Equation (LDAE)", "\n")
    print(io, "\n")
    print(io, " with vector fields")
    print(io, "\n")
    print(io, "   ϑ = ", equation.ϑ, "\n")
    print(io, "   f = ", equation.f, "\n")
    print(io, "   u = ", equation.u, "\n")
    print(io, "   g = ", equation.g, "\n")
    print(io, "   ū = ", equation.ū, "\n")
    print(io, "   ḡ = ", equation.ḡ, "\n")
    print(io, "\n")
    print(io, " and constraints")
    print(io, "\n")
    print(io, "   ϕ = ", equation.ϕ, "\n")
    print(io, "   ψ = ", equation.ψ, "\n")
    print(io, "\n")
    print(io, " Lagrangian: L = ", equation.lagrangian, "\n")
    print(io, "\n")
    print(io, " Invariants: \n")
    print(io, "   ", invariants(equation))
end

function check_initial_conditions(equ::LDAE, ics::NamedTuple)
    haskey(ics, :q) || return false
    haskey(ics, :v) || haskey(ics, :p) || return false
    haskey(ics, :λ) || return false
    if hassecondary(equ)
        haskey(ics, :μ) || return false
        eltype(ics.λ) == eltype(ics.μ) || return false
        typeof(ics.λ) == typeof(ics.μ) || return false
        axes(ics.λ) == axes(ics.μ) || return false
    end
    return true
end

function check_methods(equ::LDAE, tspan, ics::NamedTuple, params)
    applicable(equ.ϑ, zero(ics.p), tspan[begin], ics.q, zero(ics.q), params) || return false
    applicable(equ.f, zero(ics.p), tspan[begin], ics.q, zero(ics.q), params) || return false
    applicable(equ.u, zero(ics.q), tspan[begin], ics.q, vectorfield(ics.q), ics.p, ics.λ, params) || return false
    applicable(equ.g, zero(ics.p), tspan[begin], ics.q, vectorfield(ics.q), ics.p, ics.λ, params) || return false
    applicable(equ.ϕ, zero(ics.λ), tspan[begin], ics.q, vectorfield(ics.q), ics.p, params) || return false
    applicable(equ.v̄, zero(ics.q), tspan[begin], ics.q, ics.p, params) || return false
    applicable(equ.f̄, zero(ics.p), tspan[begin], ics.q, vectorfield(ics.q), params) || return false
    applicable(equ.lagrangian, tspan[begin], ics.q, vectorfield(ics.q), params) || return false
    equ.ū === nothing || applicable(equ.ū, zero(ics.q), tspan[begin], ics.q, vectorfield(ics.q), ics.p, ics.λ, params) || return false
    equ.ḡ === nothing || applicable(equ.ḡ, zero(ics.p), tspan[begin], ics.q, vectorfield(ics.q), ics.p, ics.λ, params) || return false
    equ.ψ === nothing || applicable(equ.ψ, zero(ics.λ), tspan[begin], ics.q, vectorfield(ics.q), ics.p, vectorfield(ics.q), vectorfield(ics.p), params) || return false
    # TODO add missing methods
    return true
end

function GeometricBase.datatype(equ::LDAE, ics::NamedTuple)
    @assert check_initial_conditions(equ, ics)
    return eltype(ics.q)
end

function GeometricBase.arrtype(equ::LDAE, ics::NamedTuple)
    @assert check_initial_conditions(equ, ics)
    return typeof(ics.q)
end

_get_ϑ(equ::LDAE, params) = (ϑ, t, q, v)       -> equ.ϑ(ϑ, t, q, v, params)
_get_f(equ::LDAE, params) = (f, t, q, v)       -> equ.f(f, t, q, v, params)
_get_u(equ::LDAE, params) = (u, t, q, v, p, λ) -> equ.u(u, t, q, v, p, λ, params)
_get_g(equ::LDAE, params) = (g, t, q, v, p, λ) -> equ.g(g, t, q, v, p, λ, params)
_get_ϕ(equ::LDAE, params) = (ϕ, t, q, v, p)    -> equ.ϕ(ϕ, t, q, v, p, params)
_get_ū(equ::LDAE, params) = (u, t, q, v, p, λ)    -> equ.ū(u, t, q, v, p, λ, params)
_get_ḡ(equ::LDAE, params) = (g, t, q, v, p, λ)    -> equ.ḡ(g, t, q, v, p, λ, params)
_get_ψ(equ::LDAE, params) = (ψ, t, q, v, p, q̇, ṗ) -> equ.ψ(ψ, t, q, v, p, q̇, ṗ, params)
_get_v̄(equ::LDAE, params) = (v, t, q, p)       -> equ.v̄(v, t, q, p, params)
_get_f̄(equ::LDAE, params) = (f, t, q, v)       -> equ.f̄(f, t, q, v, params)
_get_ω(equ::LDAE, params) = (ω, t, q, v)       -> equ.ω(ω, t, q, v, params)
_get_l(equ::LDAE, params) = (t, q, v)          -> equ.lagrangian(t, q, v, params)
_get_invariant(::LDAE, inv, params) = (t,q,v) -> inv(t, q, v, params)

function _functions(equ::LDAE)
    if hassecondary(equ)
        (ϑ = equ.ϑ, f = equ.f, u = equ.u, g = equ.g, ϕ = equ.ϕ, ū = equ.ū, ḡ = equ.ḡ, ψ = equ.ψ, v̄ = equ.v̄, f̄ = equ.f̄, ω = equ.ω, l = equ.lagrangian)
    else
        (ϑ = equ.ϑ, f = equ.f, u = equ.u, g = equ.g, ϕ = equ.ϕ, v̄ = equ.v̄, f̄ = equ.f̄, ω = equ.ω, l = equ.lagrangian)
    end
end

function _functions(equ::LDAE, params::OptionalParameters)
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
            f̄ = _get_f̄(equ, params),
            ω = _get_ω(equ, params),
            l = _get_l(equ, params)
        )
    else
        (
            ϑ = _get_ϑ(equ, params),
            f = _get_f(equ, params),
            u = _get_u(equ, params),
            g = _get_g(equ, params),
            ϕ = _get_ϕ(equ, params),
            v̄ = _get_v̄(equ, params),
            f̄ = _get_f̄(equ, params),
            ω = _get_ω(equ, params),
            l = _get_l(equ, params)
        )
    end
end


@doc """
`LDAEProblem`: Lagrangian Differential Algebraic Equation Problem

$(ldae_equations)

The dynamical variables ``(q,p)`` with initial conditions ``(q(t_{0}) = q_{0}, p(t_{0}) = p_{0})``
take values in ``\\mathbb{R}^{d} \\times \\mathbb{R}^{d}``. The algebraic variables ``(λ,μ)``
with initial condition ``(λ(t_{0}) = λ_{0}, μ(t_{0}) = μ_{0})`` take values in ``\\mathbb{R}^{m} \\times \\mathbb{R}^{m}``.

### Constructors

```julia
LDAEProblem(ϑ, f, u, g, ϕ, ū, ḡ, ψ, ω, l, tspan, tstep, ics; kwargs...)
LDAEProblem(ϑ, f, u, g, ϕ, ū, ḡ, ψ, ω, l, tspan, tstep, q₀::State, p₀::State, λ₀::State = zero(q₀), μ₀::State = zero(λ₀); kwargs...)
LDAEProblem(ϑ, f, u, g, ϕ, ω, l, tspan, tstep, ics; kwargs...)
LDAEProblem(ϑ, f, u, g, ϕ, ω, l, tspan, tstep, q₀::State, p₀::State, λ₀::State = zero(q₀); kwargs...)
```

$(ldae_constructors)

`tspan` is the time interval `(t₀,t₁)` for the problem to be solved in,
`tstep` is the time step to be used in the simulation, and
`ics` is a `NamedTuple` with entries `q` and `p`.
The initial conditions `q₀`, `p₀`, `λ₀` and `μ₀` can also be prescribed
directly, with `State` an `AbstractArray{<:Number}`.
For the interfaces of the functions `ϑ`, `f`, `u`, `g`, `ϕ`, `ū`, `ḡ`, `ψ`, `ω` and `l` see [`LDAE`](@ref).

In addition to the standard keyword arguments for [`EquationProblem`](@ref GeometricEquations.EquationProblem) subtypes,
a `LDAEProblem` accepts functions `v̄` and `f̄` for the computation of initial guesses for the vector fields with default
values `v̄ = _ldae_default_v̄` and `f̄ = f`.

### Function Definitions

$(ldae_functions)

With the above function definitions the `LDAEProblem` can be created by

```julia
tspan = (0.0, 1.0)
tstep = 0.1
q₀ = [1., 1.]
p₀ = [1., 0.]
λ₀ = [0.]
μ₀ = [0.]

prob = LDAEProblem(ϑ, f, u, g, ϕ, ω, l, tspan, tstep, q₀, p₀, λ₀)
```
or
```julia
prob = LDAEProblem(ϑ, f, u, g, ϕ, ū, ḡ, ψ, ω, l, tspan, tstep, q₀, p₀, λ₀, μ₀)
```    
"""
const LDAEProblem = EquationProblem{LDAE}

function LDAEProblem(ϑ, f, u, g, ϕ, ū, ḡ, ψ, ω, lagrangian, tspan, tstep, ics::NamedTuple;
                     v̄ = _ldae_default_v̄, f̄ = f, invariants = NullInvariants(),
                     parameters = NullParameters(), periodicity = NullPeriodicity())
    equ = LDAE(ϑ, f, u, g, ϕ, ū, ḡ, ψ, ω, v̄, f̄, lagrangian, invariants,
               parameter_types(parameters), periodicity)
    EquationProblem(equ, tspan, tstep, ics, parameters)
end

function LDAEProblem(ϑ, f, u, g, ϕ, ū, ḡ, ψ, ω, lagrangian, tspan, tstep, q₀::State,
                     p₀::State, λ₀::State, μ₀::State = zero(λ₀); kwargs...)
    ics = (q = q₀, p = p₀, λ = λ₀, μ = μ₀)
    LDAEProblem(ϑ, f, u, g, ϕ, ū, ḡ, ψ, ω, lagrangian, tspan, tstep, ics; kwargs...)
end

function LDAEProblem(ϑ, f, u, g, ϕ, ω, lagrangian, tspan, tstep, ics::NamedTuple; kwargs...)
    LDAEProblem(ϑ, f, u, g, ϕ, nothing, nothing, nothing, ω, lagrangian, tspan, tstep, ics;
                kwargs...)
end

function LDAEProblem(ϑ, f, u, g, ϕ, ω, lagrangian, tspan, tstep, q₀::State, p₀::State,
                     λ₀::State; kwargs...)
    ics = (q = q₀, p = p₀, λ = λ₀)
    LDAEProblem(ϑ, f, u, g, ϕ, ω, lagrangian, tspan, tstep, ics; kwargs...)
end

function GeometricBase.periodicity(prob::LDAEProblem)
    (q = periodicity(equation(prob)), p = NullPeriodicity(), λ = NullPeriodicity(),
     μ = NullPeriodicity())
end


const LDAEEnsemble  = EnsembleProblem{LDAE}

function LDAEEnsemble(ϑ, f, u, g, ϕ, ū, ḡ, ψ, ω, lagrangian, tspan, tstep, ics::AbstractVector{<:NamedTuple}; v̄ = _ldae_default_v̄, f̄ = f,
        invariants = NullInvariants(),
        parameters = NullParameters(),
        periodicity = NullPeriodicity())
    equ = LDAE(ϑ, f, u, g, ϕ, ū, ḡ, ψ, ω, v̄, f̄, lagrangian, invariants, parameter_types(parameters), periodicity)
    EnsembleProblem(equ, tspan, tstep, ics, parameters)
end

function LDAEEnsemble(ϑ, f, u, g, ϕ, ω, lagrangian, tspan, tstep, ics::AbstractVector{<:NamedTuple}; kwargs...)
    LDAEEnsemble(ϑ, f, u, g, ϕ, nothing, nothing, nothing, ω, lagrangian, tspan, tstep, ics; kwargs...)
end
