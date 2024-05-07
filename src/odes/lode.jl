
const lode_equations = raw"""
A Lagrangian system of equations is a special case of an implicit ordinary differential equations,
that is an implicit initial value problem of the form
```math
\begin{aligned}
\dot{q} (t) &= v(t) , \\
\dot{p} (t) &= f(t, q(t), v(t)) , \\
p(t) &= ϑ(t, q(t), v(t)) ,
\end{aligned}
```
with momentum ``p`` and force field ``f``, given by
```math
\begin{aligned}
p &= \frac{\partial L}{\partial v} , &
f &= \frac{\partial L}{\partial q} .
\end{aligned}
```
This is a special case of an implicit ordinary differential equation, that
is defined by a Lagrangian, as well as a special case of a differential algebraic
equation with dynamical variables ``(q,p)`` and algebraic variable ``v``, that is
determined such that the constraint ``p(t) = ϑ(t, q(t), v(t))`` is satisfied.

Many integrators perform a projection step in order to enforce this constraint. To this end,
the system is extended to
```math
\begin{aligned}
\dot{q} (t) &= v(t) + \lambda(t) , \\
\dot{p} (t) &= f(t, q(t), v(t)) + g(t, q(t), v(t), \lambda(t)) , \\
p(t) &= ϑ(t, q(t), v(t)) ,
\end{aligned}
```
where the vector field defining the projection step is usually given as
```math
\begin{aligned}
g(t, q(t), v(t), λ(t)) &= λ(t) \cdot \nabla ϑ(t, q(t), v(t)) .
\end{aligned}
```
"""

const lode_functions = raw"""
The functions `ϑ` and `f` must have the interface
```julia
function ϑ(p, t, q, v)
    p[1] = ...
    p[2] = ...
    ...
end
```
and
```julia
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
function g(g, t, q, λ)
    g[1] = ...
    g[2] = ...
    ...
end
```
and
```julia
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
The function `g` is used in projection methods that enforce ``p = ϑ(q)``.
The functions `v̄` and `f̄` are used for initial guesses in nonlinear implicit solvers.
Finally, the functions `ω` and `l`, computing the symplectic matrix and the Lagrangian,
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
"""


@doc """
`LODE`: Lagrangian Ordinary Differential Equation

$(lode_equations)

### Parameters

* `ϑType <: Function`: type of `ϑ`
* `fType <: Function`: type of `f`
* `gType <: Function`: type of `g`
* `ωType <: Function`: type of `ω`
* `v̄Type <: Function`: type of `v̄`
* `f̄Type <: Function`: type of `f̄`
* `lagType <: Function`: Lagrangian type
* `invType <: OptionalInvariants`: invariants type
* `parType <: OptionalParameters`: parameters type
* `perType <: OptionalPeriodicity`: periodicity type

### Fields

* `ϑ`: function determining the momentum
* `f`: function computing the vector field
* `g`: function determining the projection, given by ``\\nabla \\vartheta (q) \\cdot \\lambda``
* `ω`: function computing the symplectic matrix
* `v̄`: function computing an initial guess for the velocity field ``v`` (optional)
* `f̄`: function computing an initial guess for the force field ``f`` (optional)
* `lagrangian`: function computing the Lagrangian ``L``
* `invariants`: functions for the computation of invariants, either a `NamedTuple` containing the equation's invariants or `NullInvariants`
* `parameters`: type constraints for parameters, either a `NamedTuple` containing the equation's parameters or `NullParameters`
* `periodicity`: determines the periodicity of the state vector `q` for cutting periodic solutions, either a `AbstractArray` or `NullPeriodicity`

### Constructors

```julia
LODE(ϑ, f, g, ω, l, v̄, f̄, invariants, parameters, periodicity)
LODE(ϑ, f, g, ω, l; v̄ = _lode_default_v̄, f̄ = f, invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
```

where 

```julia
_lode_default_v̄(v, t, q, params) = nothing
```

### Function Definitions

$(lode_functions)

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

_lode_default_g(g, t, q, v, λ, params) = nothing
_lode_default_v̄(v, t, q, p, params) = nothing

LODE(ϑ, f, g, ω, l; invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity(), v̄=_lode_default_v̄, f̄=f) = LODE(ϑ, f, g, ω, v̄, f̄, l, invariants, parameters, periodicity)

GeometricBase.invariants(equation::LODE) = equation.invariants
GeometricBase.parameters(equation::LODE) = equation.parameters
GeometricBase.periodicity(equation::LODE) = equation.periodicity

hasvectorfield(::LODE) = true
haslagrangian(::LODE) = true

function Base.show(io::IO, equation::LODE)
    print(io, "Lagrangian Ordinary Differential Equation (LODE)", "\n")
    print(io, "\n")
    print(io, " with vector fields", "\n")
    print(io, "   ϑ = ", equation.ϑ, "\n")
    print(io, "   f = ", equation.f, "\n")
    print(io, "   g = ", equation.g, "\n")
    print(io, "\n")
    print(io, " Lagrangian: L = ", equation.lagrangian, "\n")
    print(io, "\n")
    print(io, " Invariants: \n")
    print(io, "   ", invariants(equation))
end

function initialstate(::LODE, q₀::InitialState, p₀::InitialState, λ₀::InitialAlgebraic = zeroalgebraic(q₀))
    (
        q = _statevariable(q₀),
        p = _statevariable(p₀),
        λ = _algebraicvariable(λ₀),
    )
end

function initialstate(::LODE, q₀::InitialStateVector, p₀::InitialStateVector, λ₀::InitialAlgebraicVector = zeroalgebraic(q₀))
    [(
        q = _statevariable(q),
        p = _statevariable(p),
        λ = _algebraicvariable(λ)
    ) for (q,p,λ) in zip(q₀,p₀,λ₀)]
end

function check_initial_conditions(::LODE, ics::NamedTuple)
    haskey(ics, :q) || return false
    haskey(ics, :v) || haskey(ics, :p) || return false
    haskey(ics, :λ) || return false
    eltype(ics.q) == eltype(ics.p) == eltype(ics.λ) || return false
    axes(ics.q) == axes(ics.p) == axes(ics.λ) || return false
    typeof(ics.q) <: StateVariable || return false
    typeof(ics.p) <: StateVariable || return false
    typeof(ics.λ) <: AlgebraicVariable || return false
    return true
end

function check_methods(equ::LODE, tspan, ics, params)
    applicable(equ.ϑ, zero(ics.p), tspan[begin], ics.q, vectorfield(ics.q), params) || return false
    applicable(equ.f, vectorfield(ics.p), tspan[begin], ics.q, vectorfield(ics.q), params) || return false
    applicable(equ.g, vectorfield(ics.p), tspan[begin], ics.q, vectorfield(ics.q), ics.λ, params) || return false
    # TODO add missing methods (namely ω)
    applicable(equ.v̄, vectorfield(ics.q), tspan[begin], ics.q, ics.p, params) || return false
    applicable(equ.f̄, vectorfield(ics.p), tspan[begin], ics.q, vectorfield(ics.q), params) || return false
    applicable(equ.lagrangian, tspan[begin], ics.q, vectorfield(ics.q), params) || return false
    return true
end

function GeometricBase.datatype(equ::LODE, ics::NamedTuple)
    @assert check_initial_conditions(equ, ics)
    return eltype(ics.q)
end

function GeometricBase.arrtype(equ::LODE, ics::NamedTuple)
    @assert check_initial_conditions(equ, ics)
    typeof(ics.q)
end

_get_ϑ(equ::LODE, params) = (ϑ, t, q, v) -> equ.ϑ(ϑ, t, q, v, params)
_get_f(equ::LODE, params) = (f, t, q, v) -> equ.f(f, t, q, v, params)
_get_g(equ::LODE, params) = (g, t, q, v, λ) -> equ.g(g, t, q, v, λ, params)
_get_ω(equ::LODE, params) = (ω, t, q, v) -> equ.ω(ω, t, q, v, params)
_get_v̄(equ::LODE, params) = (v, t, q, p) -> equ.v̄(v, t, q, p, params)
_get_f̄(equ::LODE, params) = (f, t, q, v) -> equ.f̄(f, t, q, v, params)
_get_l(equ::LODE, params) = (t, q, v)    -> equ.lagrangian(t, q, v, params)
_get_invariant(::LODE, inv, params) = (t, q, v) -> inv(t, q, v, params)

_functions(equ::LODE) = (ϑ = equ.ϑ, f = equ.f, g = equ.g, ω = equ.ω, v̄ = equ.v̄, f̄ = equ.f̄, l = equ.lagrangian)

_functions(equ::LODE, params::OptionalParameters) = (
        ϑ = _get_ϑ(equ, params),
        f = _get_f(equ, params),
        g = _get_g(equ, params),
        ω = _get_ω(equ, params),
        v̄ = _get_v̄(equ, params),
        f̄ = _get_f̄(equ, params),
        l = _get_l(equ, params)
    )


@doc """
`LODEProblem`: Lagrangian Ordinary Differential Equation Problem

$(lode_equations)

The dynamical variables ``(q,p)`` with initial conditions ``(q(t_{0}) = q_{0}, p(t_{0}) = p_{0})``
take values in ``T^{*} Q \\simeq \\mathbb{R}^{d} \\times \\mathbb{R}^{d}``. The algebraic variable ``λ``
with initial condition ``λ(t_{0}) = λ_{0}`` takes values in ``\\mathbb{R}^{m}``.

### Constructors

```julia
LODEProblem(ϑ, f, ω, l, tspan, tstep, ics; kwargs...)
LODEProblem(ϑ, f, ω, l, tspan, tstep, q₀::StateVariable, p₀::StateVariable, λ₀::AlgebraicVariable; kwargs...)
LODEProblem(ϑ, f, ω, l, tspan, tstep, q₀::AbstractArray, p₀::AbstractArray, λ₀::AbstractArray = zero(q₀); kwargs...)
LODEProblem(ϑ, f, g, ω, l, tspan, tstep, ics; kwargs...)
LODEProblem(ϑ, f, g, ω, l, tspan, tstep, q₀::StateVariable, p₀::StateVariable, λ₀::AlgebraicVariable; kwargs...)
LODEProblem(ϑ, f, g, ω, l, tspan, tstep, q₀::AbstractArray, p₀::AbstractArray, λ₀::AbstractArray = zero(q₀); kwargs...)
```
where `ϑ`, `f` and `g` are the functions computing the momentum and the vector fields, respectively,
`ω` determines the symplectic matrix, and `l` returns the Lagrangian,
`tspan` is the time interval `(t₀,t₁)` for the problem to be solved in,
`tstep` is the time step to be used in the simulation, and
`ics` is a `NamedTuple` with entries `q`, `p` and `λ`.
The initial conditions `q₀`, `p₀` and `λ₀` can also be prescribed
directly, with `StateVariable` an `AbstractArray{<:Number}`, where `λ₀` can also be omitted.
For the interfaces of the functions `ϑ`, `f`, `g`, `ω` and `l` see [`LODE`](@ref).

In addition to the standard keyword arguments for [`EquationProblem`](@ref GeometricEquations.EquationProblem) subtypes,
a `LODEProblem` accepts functions `v̄` and `f̄` for the computation of initial guesses for the vector fields with default
values `v̄ = _lode_default_v̄` and `f̄ = f`.

"""
const LODEProblem = EquationProblem{LODE}

function LODEProblem(ϑ, f, g, ω, l, tspan::Tuple, tstep::Real, ics...;
        invariants = NullInvariants(),
        parameters = NullParameters(),
        periodicity = NullPeriodicity(),
        v̄ = _lode_default_v̄, f̄ = f)
    equ = LODE(ϑ, f, g, ω, v̄, f̄, l, invariants, parameter_types(parameters), periodicity)
    EquationProblem(equ, tspan, tstep, initialstate(equ, ics...), parameters)
end

function LODEProblem(ϑ, f, ω, l, args...; kwargs...)
    LODEProblem(ϑ, f, _lode_default_g, ω, l, args...; kwargs...)
end


function GeometricBase.periodicity(prob::LODEProblem)
    (q = periodicity(equation(prob)), p = NullPeriodicity(), λ = NullPeriodicity())
end

@inline GeometricBase.nconstraints(prob::LODEProblem) = ndims(prob)


@doc """
`LODEEnsemble`: Lagrangian Ordinary Differential Equation Ensemble

$(lode_equations)

The dynamical variables ``(q,p)`` take values in ``T^{*} Q \\simeq \\mathbb{R}^{d} \\times \\mathbb{R}^{d}``.
The algebraic variable ``λ`` takes values in ``\\mathbb{R}^{m}``.

### Constructors

```julia
LODEProblem(ϑ, f, ω, l, tspan, tstep, ics; kwargs...)
LODEProblem(ϑ, f, ω, l, tspan, tstep, q₀::StateVariable, p₀::StateVariable, λ₀::AlgebraicVariable; kwargs...)
LODEProblem(ϑ, f, ω, l, tspan, tstep, q₀::AbstractArray, p₀::AbstractArray, λ₀::AbstractArray = zero(q₀); kwargs...)
LODEProblem(ϑ, f, g, ω, l, tspan, tstep, ics; kwargs...)
LODEProblem(ϑ, f, g, ω, l, tspan, tstep, q₀::StateVariable, p₀::StateVariable, λ₀::AlgebraicVariable; kwargs...)
LODEProblem(ϑ, f, g, ω, l, tspan, tstep, q₀::AbstractArray, p₀::AbstractArray, λ₀::AbstractArray = zero(q₀); kwargs...)
```
where `ϑ`, `f` and `g` are the functions computing the momentum and the vector fields, respectively,
`ω` determines the symplectic matrix, and `l` returns the Lagrangian,
`tspan` is the time interval `(t₀,t₁)` for the problem to be solved in,
`tstep` is the time step to be used in the simulation, and
`ics` is a `NamedTuple` with entries `q`, `p` and `λ`.
The initial conditions `q₀`, `p₀` and `λ₀` can also be prescribed
directly, with `StateVariable` an `AbstractArray{<:Number}`, where `λ₀` can also be omitted.
For the interfaces of the functions `ϑ`, `f`, `g`, `ω` and `l` see [`LODE`](@ref).

In addition to the standard keyword arguments for [`EquationProblem`](@ref GeometricEquations.EquationProblem) subtypes,
a `LODEProblem` accepts functions `v̄` and `f̄` for the computation of initial guesses for the vector fields with default
values `v̄ = _lode_default_v̄` and `f̄ = f`.

"""
const LODEEnsemble = EnsembleProblem{LODE}

function LODEEnsemble(ϑ, f, g, ω, l, tspan::Tuple, tstep::Real, ics...;
        invariants = NullInvariants(),
        parameters = NullParameters(),
        periodicity = NullPeriodicity(),
        v̄ = _lode_default_v̄,
        f̄ = f
    )
    equ = LODE(ϑ, f, g, ω, v̄, f̄, l, invariants, parameter_types(parameters), periodicity)
    EnsembleProblem(equ, tspan, tstep, initialstate(equ, ics...), parameters)
end

function LODEEnsemble(ϑ, f, ω, l, args...; kwargs...)
    LODEEnsemble(ϑ, f, _lode_default_g, ω, l, args...; kwargs...)
end
