
const pdae_equations = raw"""
A partitioned differential algebraic equation has the form
```math
\begin{aligned}
\dot{q} (t) &= v(t, q(t), p(t)) + u(t, q(t), p(t), \lambda(t)) , & q(t_{0}) &= q_{0} , \\
\dot{p} (t) &= f(t, q(t), p(t)) + g(t, q(t), p(t), \lambda(t)) , & p(t_{0}) &= p_{0} , \\
0 &= \phi (t, q(t), p(t)) , & \lambda(t_{0}) &= \lambda_{0} ,
\end{aligned}
```
with vector fields ``v`` and ``f``, projection ``u`` and ``g``,
algebraic constraint ``\phi=0``, initial conditions ``(q_{0}, p_{0})``
and ``\lambda_{0}``, the dynamical variables ``(q,p)`` taking values
in ``\mathbb{R}^{d} \times \mathbb{R}^{d}`` and the algebraic variable
``\lambda`` taking values in ``\mathbb{R}^{m}``.

Some integrators also enforce the secondary constraint ``\psi``, that is the time
derivative of the algebraic constraint ``\phi``.
In this case, the system of equations is modified as follows
```math
\begin{aligned}
\dot{q} (t) &= v(t, q(t), p(t)) + u(t, q(t), p(t), \lambda(t)) + \bar{u} (t, q(t), p(t), \gamma(t)) , & q(t_{0}) &= q_{0} , \\
\dot{p} (t) &= f(t, q(t), p(t)) + g(t, q(t), p(t), \lambda(t)) + \bar{g} (t, q(t), p(t), \gamma(t)) , & p(t_{0}) &= p_{0} , \\
0 &= \phi (t, q(t), p(t)) , & \lambda(t_{0}) &= \lambda_{0} , \\
0 &= \psi (t, q(t), p(t), \dot{q} (t), \dot{p} (t)) , & \gamma(t_{0}) &= \gamma_{0} ,
\end{aligned}
```
with the second algebraic variable ``\gamma`` also taking values in ``\mathbb{R}^{m}``.
"""

const pdae_constructors = raw"""
The functions `v` and `f` compute the vector field, `u` and `g` compute the projections,
and `ϕ` provides the algebraic constraint.
The functions `ψ`, `ū` and `ḡ` are optional and provide the secondary constraint and the
corresponding projection.
"""

const pdae_functions = raw"""
The functions `v`, `f`, `u`, `g` and `ϕ` must have the interface
```julia
function v(v, t, q, p, params)
    v[1] = ...
    v[2] = ...
    ...
end

function f(g, t, q, p, params)
    f[1] = ...
    f[2] = ...
    ...
end

function u(u, t, q, p, λ, params)
    u[1] = ...
    u[2] = ...
    ...
end

function g(g, t, q, p, λ, params)
    g[1] = ...
    g[2] = ...
    ...
end

function ϕ(ϕ, t, q, p, params)
    ϕ[1] = ...
end
```
where `t` is the current time, `q`, `p` and `λ` are the current solution vectors,
`v`, `f`, `u` and `g` are the vectors which hold the result of evaluating the
vector fields ``v`` and ``f``, the projections ``u`` and ``g``, and `ϕ` holds the
algebraic constraint ``\phi``, evaluated on `t`, `q`, `p` and `λ`.

Some integrators also enforce the secondary constraint ``\psi`` and require
the following additional functions
```
function ū(u, t, q, p, γ, params)
    u[1] = ...
    u[2] = ...
    ...
end

function ḡ(g, t, q, p, γ, params)
    g[1] = ...
    g[2] = ...
    ...
end

function ψ(ψ, t, q, p, v, f, params)
    ψ[1] = ...
end
```
"""


@doc """
`PDAE`: Partitioned Differential Algebraic Equation

$(pdae_equations)

### Parameters

* `vType <: Callable`: type of `v`
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

* `v`: function computing the vector field ``v``
* `f`: function computing the vector field ``f``
* `u`: function computing the projection for ``q``
* `g`: function computing the projection for ``p``
* `ϕ`: algebraic constraints
* `ū`: function computing the secondary projection field ``\\bar{u}`` (*optional*)
* `ḡ`: function computing the secondary projection field ``\\bar{g}`` (*optional*)
* `ψ`: secondary constraints (*optional*)
* `v̄`: function computing an initial guess for the velocity field ``v`` (*optional*, defaults to `v`)
* `f̄`: function computing an initial guess for the force field ``f`` (*optional*, defaults to `f`)
* `invariants`: functions for the computation of invariants, either a `NamedTuple` containing the equation's invariants or `NullInvariants`
* `parameters`: type constraints for parameters, either a `NamedTuple` containing the equation's parameters or `NullParameters`
* `periodicity`: determines the periodicity of the state vector `q` for cutting periodic solutions, either a `AbstractArray` or `NullPeriodicity`

### Constructors

```julia
PDAE(v, f, u, g, ϕ, ū, ḡ, ψ, v̄, f̄, invariants, parameters, periodicity)
PDAE(v, f, u, g, ϕ, ū, ḡ, ψ; kwargs...)
PDAE(v, f, u, g, ϕ; kwargs...)
```

$(pdae_constructors)

### Function Definitions

The functions are defined by

$(pdae_functions)

The `PDAE` is created by

```julia
equ = PDAE(v, f, u, g, ϕ)
```
or
```julia
equ = PDAE(v, f, u, g, ϕ, ū, ḡ, ψ)
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
    applicable(equ.v, zero(ics.q), tspan[begin], ics.q, ics.p, params) || return false
    applicable(equ.f, zero(ics.p), tspan[begin], ics.q, ics.p, params) || return false
    applicable(equ.u, zero(ics.q), tspan[begin], ics.q, ics.p, ics.λ, params) || return false
    applicable(equ.g, zero(ics.p), tspan[begin], ics.q, ics.p, ics.λ, params) || return false
    applicable(equ.ϕ, zero(ics.λ), tspan[begin], ics.q, ics.p, params) || return false
    # TODO add missing methods
    return true
end

function GeometricBase.datatype(equ::PDAE, ics::NamedTuple)
    @assert check_initial_conditions(equ, ics)
    return eltype(ics.q)
end

function GeometricBase.arrtype(equ::PDAE, ics::NamedTuple)
    @assert check_initial_conditions(equ, ics)
    return typeof(ics.q)
end

_get_v(equ::PDAE, params) = (v, t, q, p)       -> equ.v(v, t, q, p, params)
_get_f(equ::PDAE, params) = (f, t, q, p)       -> equ.f(f, t, q, p, params)
_get_u(equ::PDAE, params) = (u, t, q, p, λ)    -> equ.u(u, t, q, p, λ, params)
_get_g(equ::PDAE, params) = (g, t, q, p, λ)    -> equ.g(g, t, q, p, λ, params)
_get_ϕ(equ::PDAE, params) = (ϕ, t, q, p)       -> equ.ϕ(ϕ, t, q, p, params)
_get_ū(equ::PDAE, params) = (u, t, q, p, λ)    -> equ.ū(u, t, q, p, λ, params)
_get_ḡ(equ::PDAE, params) = (g, t, q, p, λ)    -> equ.ḡ(g, t, q, p, λ, params)
_get_ψ(equ::PDAE, params) = (ψ, t, q, p, v, f) -> equ.ψ(ψ, t, q, p, v, f, params)
_get_v̄(equ::PDAE, params) = (v, t, q, p)       -> equ.v̄(v, t, q, p, params)
_get_f̄(equ::PDAE, params) = (f, t, q, p)       -> equ.f̄(f, t, q, p, params)
_get_invariant(::PDAE, inv, params) = (t, q, p) -> inv(t, q, p, params)

function _functions(equ::PDAE)
    if hassecondary(equ)
        (v = equ.v, f = equ.f, u = equ.u, g = equ.g, ϕ = equ.ϕ, ū = equ.ū, ḡ = equ.ḡ, ψ = equ.ψ, v̄ = equ.v̄, f̄ = equ.f̄)
    else
        (v = equ.v, f = equ.f, u = equ.u, g = equ.g, ϕ = equ.ϕ, v̄ = equ.v̄, f̄ = equ.f̄)
    end
end

function _functions(equ::PDAE, params::OptionalParameters)
    if hassecondary(equ)
        (
            v = _get_v(equ, params),
            f = _get_f(equ, params),
            u = _get_u(equ, params),
            g = _get_g(equ, params),
            ϕ = _get_ϕ(equ, params),
            ū = _get_ū(equ, params),
            ḡ = _get_ḡ(equ, params),
            ψ = _get_ψ(equ, params),
            v̄ = _get_v̄(equ, params),
            f̄ = _get_f̄(equ, params))
    else
        (
            v = _get_v(equ, params),
            f = _get_f(equ, params),
            u = _get_u(equ, params),
            g = _get_g(equ, params),
            ϕ = _get_ϕ(equ, params),
            v̄ = _get_v̄(equ, params),
            f̄ = _get_f̄(equ, params))
    end
end


@doc """
`PDAEProblem`: Partitioned Differential Algebraic Equation Problem

$(pdae_equations)

### Constructors

```julia
PDAEProblem(v, f, u, g, ϕ, ū, ḡ, ψ, tspan, tstep, ics::NamedTuple; kwargs...)
PDAEProblem(v, f, u, g, ϕ, ū, ḡ, ψ, tspan, tstep, q₀::State, p₀::State, λ₀::State; kwargs...)
PDAEProblem(v, f, u, g, ϕ, tspan, tstep, ics::NamedTuple; kwargs...)
PDAEProblem(v, f, u, g, ϕ, tspan, tstep, q₀::State, p₀::State, λ₀::State; kwargs...)
```

$(pdae_constructors)

`tspan` is the time interval `(t₀,t₁)` for the problem to be solved in,
`tstep` is the time step to be used in the simulation, and
`ics` is a `NamedTuple` with entries `q`, `p` and `λ`.
The initial conditions `q₀`, `p₀` and `λ₀` can also be prescribed directly,
with `State` an `AbstractArray{<:Number}`.

In addition to the standard keyword arguments for [`GeometricProblem`](@ref GeometricEquations.GeometricProblem) subtypes,
a `PDAEProblem` accepts functions `v̄` and `f̄` for the computation of initial guesses for the vector fields with default
values `v̄ = v` and `f̄ = f`.
    
### Function Definitions

$(pdae_functions)

With the above function definitions the `PDAEProblem` can be created by

```julia
tspan = (0.0, 1.0)
tstep = 0.1
q₀ = [1., 1.]
p₀ = [1., 0.]
λ₀ = [0.]
γ₀ = [0.]

prob = PDAEProblem(v, f, u, g, ϕ, tspan, tstep, q₀, p₀, λ₀)
```
or
```julia
prob = PDAEProblem(v, f, u, g, ϕ, ū, ḡ, ψ, tspan, tstep, q₀, p₀, λ₀, γ₀)
```
"""
const PDAEProblem = GeometricProblem{PDAE}

function PDAEProblem(v, f, u, g, ϕ, ū, ḡ, ψ, tspan, tstep, ics::NamedTuple; v̄ = v,
                     f̄ = f, invariants = NullInvariants(), parameters = NullParameters(),
                     periodicity = NullPeriodicity())
    equ = PDAE(v, f, u, g, ϕ, ū, ḡ, ψ, v̄, f̄, invariants, parameter_types(parameters),
               periodicity)
    GeometricProblem(equ, tspan, tstep, ics, parameters)
end

function PDAEProblem(v, f, u, g, ϕ, ū, ḡ, ψ, tspan, tstep, q₀::State, p₀::State,
                     λ₀::State, μ₀::State = zero(λ₀); kwargs...)
    ics = (q = q₀, p = p₀, λ = λ₀, μ = μ₀)
    PDAEProblem(v, f, u, g, ϕ, ū, ḡ, ψ, tspan, tstep, ics; kwargs...)
end

function PDAEProblem(v, f, u, g, ϕ, tspan, tstep, ics::NamedTuple; kwargs...)
    PDAEProblem(v, f, u, g, ϕ, nothing, nothing, nothing, tspan, tstep, ics; kwargs...)
end

function PDAEProblem(v, f, u, g, ϕ, tspan, tstep, q₀::State, p₀::State, λ₀::State;
                     kwargs...)
    ics = (q = q₀, p = p₀, λ = λ₀)
    PDAEProblem(v, f, u, g, ϕ, tspan, tstep, ics; kwargs...)
end

function GeometricBase.periodicity(prob::PDAEProblem)
    (q = periodicity(equation(prob)), p = NullPeriodicity(), λ = NullPeriodicity())
end
