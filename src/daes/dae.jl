
const dae_equations = raw"""
Defines a differential algebraic initial value problem
```math
\begin{aligned}
\dot{q} (t) &= v(t, q(t)) + u(t, q(t), \lambda(t)) , \\
0 &= \phi (t, q(t)) ,
\end{aligned}
```
with vector field ``v``, projection ``u``, algebraic constraint ``\phi=0``.

Some integrators also enforce the secondary constraint ``\psi``, that is the time
derivative of the algebraic constraint ``\phi``.
In this case, the system of equations is modified as follows
```math
\begin{aligned}
\dot{q} (t) &= v(t, q(t)) + u(t, q(t), \lambda(t)) + \bar{u} (t, q(t), \dot{q} (t), \dot{p} (t), \mu(t)) , \\
0 &= \phi (t, q(t)) , \\
0 &= \psi (t, q(t), \dot{q} (t)) .
\end{aligned}
```
"""

const dae_constructors = raw"""
The functions `v` and `u` compute the vector field and the projection, respectively,
`ϕ` provides the algebraic constraint.
The functions `ψ` and `ū` are optional and provide the secondary constraint, that is the time
derivative of the algebraic constraint, and the corresponding projection.
"""

const dae_functions = raw"""
The functions `v`, `u` and `ϕ` must have the interface
```julia
function v(v, t, q, params)
    v[1] = ...
    v[2] = ...
    ...
end

function u(u, t, q, λ, params)
    u[1] = ...
    u[2] = ...
    ...
end

function ϕ(ϕ, t, q, params)
    ϕ[1] = ...
end
```
where `t` is the current time, `q` and `λ` are the current solution vectors,
and `v`, `u` and `ϕ` are the vectors which hold the result of evaluating the
vector field ``v``, the projection ``u`` and the algebraic constraint ``\phi``
on `t`, `q` and `λ`.

Some integrators also enforce the secondary constraint ``\psi`` and require
the following additional functions
```
function ū(u, t, q, μ, params)
    u[1] = ...
    u[2] = ...
    ...
end

function ψ(ψ, t, q, v, params)
    ψ[1] = ...
end
```
"""


@doc """
`DAE`: Differential Algebraic Equation

$(dae_equations)

### Parameters

* `vType <: Callable`: type of `v`
* `uType <: Callable`: type of `u`
* `ϕType <: Callable`: type of `ϕ`
* `ūType <: OptionalCallable`: type of `ū`
* `ψType <: OptionalCallable`: type of `ψ`
* `v̄Type <: Callable`: type of `v̄`
* `invType <: OptionalInvariants`: invariants type
* `parType <: OptionalParameters`: parameters type
* `perType <: OptionalPeriodicity`: periodicity type

### Fields

* `v`: function computing the vector field `v(v, t, q, params)`
* `u`: function computing the projection `u(u, t, q, λ, params)`
* `ϕ`: algebraic constraint `ϕ(ϕ, t, q, params)`
* `ū`: function computing the secondary projection field `ū(ū, t, q, λ, params)` (*optional*)
* `ψ`: secondary constraint `ψ(ψ, t, q, v, params)` (*optional*)
* `v̄`: function computing an initial guess for the velocity field ``v`` (defaults to `v`)
* `invariants`: functions for the computation of invariants, either a `NamedTuple` containing the equation's invariants or `NullInvariants`
* `parameters`: type constraints for parameters, either a `NamedTuple` containing the equation's parameters or `NullParameters`
* `periodicity`: determines the periodicity of the state vector `q` for cutting periodic solutions, either a `AbstractArray` or `NullPeriodicity`

### Constructors

```julia
DAE(v, u, ϕ, ū, ψ, v̄, invariants, parameters, periodicity)
DAE(v, u, ϕ, ū, ψ; kwargs...)
DAE(v, u, ϕ; kwargs...)
```

$(dae_constructors)

### Function Definitions

The functions are defined by

$(dae_functions)

The `DAE` is created by

```julia
equ = DAE(v, u, ϕ)
```
or
```julia
equ = DAE(v, u, ϕ, ū, ψ)
```
"""
struct DAE{vType <: Callable,
           uType <: Callable, ϕType <: Callable,
           ūType <: OptionalCallable, ψType <: OptionalCallable,
           v̄Type <: OptionalCallable,
           invType <: OptionalInvariants,
           parType <: OptionalParameters,
           perType <: OptionalPeriodicity} <: AbstractEquationDAE{invType,parType,perType,ψType}

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

DAE(v, u, ϕ, ū, ψ, v̄; invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity()) = DAE(v, u, ϕ, ū, ψ, v̄, invariants, parameters, periodicity)
DAE(v, u, ϕ, ū, ψ; v̄=v, kwargs...) = DAE(v, u, ϕ, ū, ψ, v̄; kwargs...)
DAE(v, u, ϕ; kwargs...) = DAE(v, u, ϕ, nothing, nothing; kwargs...)

GeometricBase.invariants(equation::DAE) = equation.invariants
GeometricBase.parameters(equation::DAE) = equation.parameters
GeometricBase.periodicity(equation::DAE) = equation.periodicity

hasvectorfield(::DAE) = true
hasinitialguess(::DAE{vType, uType, ϕType, ūType, ψType, <:Callable}) where {vType, uType, ϕType, ūType, ψType} = true

function Base.show(io::IO, equation::DAE)
    print(io, "Differential Algebraic Equation (DAE)", "\n")
    print(io, "\n")
    print(io, " with vector fields")
    print(io, "\n")
    print(io, "   v = ", equation.v, "\n")
    print(io, "   u = ", equation.u, "\n")
    print(io, "   ū = ", equation.ū, "\n")
    print(io, "\n")
    print(io, " and constraints")
    print(io, "\n")
    print(io, "   ϕ = ", equation.ϕ, "\n")
    print(io, "   ψ = ", equation.ψ, "\n")
    print(io, "\n")
    print(io, " Invariants: \n")
    print(io, "   ", invariants(equation))
end

function initialstate(::DAE, t::InitialTime, ics::NamedTuple, params::OptionalParameters)
    ics = haskey(ics, :μ) ? ics : merge(ics, (μ = zeroalgebraic(ics.λ),))

    (
        q = _statevariable(ics.q),
        λ = _algebraicvariable(ics.λ),
        μ = _algebraicvariable(ics.μ),
    )
end

function initialstate(equ::DAE, q₀::InitialState, λ₀::InitialAlgebraic, μ₀::InitialAlgebraic = zeroalgebraic(λ₀))
    initialstate(equ, (q = q₀, λ = λ₀, μ = μ₀))
end

function initialstate(equ::DAE, q₀::InitialStateVector, λ₀::InitialAlgebraicVector, μ₀::InitialAlgebraicVector = zeroalgebraic(λ₀))
    [initialstate(equ, q, λ, μ) for (q,λ,μ) in zip(q₀,λ₀,μ₀)]
end

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

function check_methods(equ::DAE, timespan, ics::NamedTuple, params)
    applicable(equ.v, zero(ics.q), timespan[begin], ics.q, params) || return false
    applicable(equ.u, zero(ics.q), timespan[begin], ics.q, ics.λ, params) || return false
    applicable(equ.ϕ, zero(ics.λ), timespan[begin], ics.q, params) || return false
    equ.ū === nothing || applicable(equ.u, zero(ics.q), timespan[begin], ics.q, ics.λ, params) || return false
    equ.ψ === nothing || applicable(equ.ψ, zero(ics.λ), timespan[begin], ics.q, vectorfield(ics.q), params) || return false
    equ.v̄ === nothing || applicable(equ.v̄, zero(ics.q), timespan[begin], ics.q, params) || return false
    return true
end

function GeometricBase.datatype(equ::DAE, ics::NamedTuple)
    @assert check_initial_conditions(equ, ics)
    return eltype(ics.q)
end

function GeometricBase.arrtype(equ::DAE, ics::NamedTuple)
    @assert check_initial_conditions(equ, ics)
    return typeof(ics.q)
end

_get_v(equ::DAE, params) = (v, t, q)    -> equ.v(v, t, q, params)
_get_u(equ::DAE, params) = (u, t, q, λ) -> equ.u(u, t, q, λ, params)
_get_ϕ(equ::DAE, params) = (ϕ, t, q)    -> equ.ϕ(ϕ, t, q, params)
_get_ū(equ::DAE, params) = (u, t, q, λ) -> equ.ū(u, t, q, λ, params)
_get_ψ(equ::DAE, params) = (ψ, t, q, v) -> equ.ψ(ψ, t, q, v, params)
_get_v̄(equ::DAE, params) = (v, t, q)    -> equ.v̄(v, t, q, params)
_get_invariant(::DAE, inv, params) = (t,q) -> inv(t, q, params)

function _functions(equ::DAE)
    if hassecondary(equ)
        (v = equ.v, u = equ.u, ϕ = equ.ϕ, ū = equ.ū, ψ = equ.ψ)
    else
        (v = equ.v, u = equ.u, ϕ = equ.ϕ)
    end
end

function _functions(equ::DAE, params::OptionalParameters)
    if hassecondary(equ)
        (
            v = _get_v(equ, params),
            u = _get_u(equ, params),
            ϕ = _get_ϕ(equ, params),
            ū = _get_ū(equ, params),
            ψ = _get_ψ(equ, params),
            v̄ = _get_v̄(equ, params)
        )
    else
        (
            v = _get_v(equ, params),
            u = _get_u(equ, params),
            ϕ = _get_ϕ(equ, params),
            v̄ = _get_v̄(equ, params)
        )
    end
end

_initialguess(equ::DAE) = (v = equ.v̄,)
_initialguess(equ::DAE, params::OptionalParameters) = (v = _get_v̄(equ, params),)


@doc """
`DAEProblem`: Differential Algebraic Equation Problem

$(dae_equations)

The dynamical variable ``q`` with initial conditions ``q(t_{0}) = q_{0}``
takes values in ``\\mathbb{R}^{d}``. The algebraic variables ``(λ,μ)``
with initial condition ``(λ(t_{0}) = λ_{0}, μ(t_{0}) = μ_{0})`` take
values in ``\\mathbb{R}^{m} \\times \\mathbb{R}^{m}``.

### Constructors

```julia
DAEProblem(v, u, ϕ, ū, ψ, timespan, timestep, ics::NamedTuple; kwargs...)
DAEProblem(v, u, ϕ, ū, ψ, timespan, timestep, q₀::StateVariable, λ₀::StateVariable, μ₀::StateVariable = zero(λ₀); kwargs...)
DAEProblem(v, u, ϕ, timespan, timestep, ics::NamedTuple; kwargs...)
DAEProblem(v, u, ϕ, timespan, timestep, q₀::StateVariable, λ₀::StateVariable; kwargs...)
```

$(dae_constructors)

`timespan` is the time interval `(t₀,t₁)` for the problem to be solved in,
`timestep` is the time step to be used in the simulation, and
`ics` is a `NamedTuple` with entries `q`, `λ` and `μ`.
The initial conditions `q₀`, `λ₀` and `μ₀` can also be prescribed directly,
with `StateVariable` an `AbstractArray{<:Number}`.
For the interfaces of the functions `v`, `u`, `ϕ`, `ū`, `ψ` see [`DAE`](@ref).

In addition to the standard keyword arguments for [`EquationProblem`](@ref GeometricEquations.EquationProblem) subtypes,
a `DAEProblem` accepts a function `v̄` for the computation of an initial guess for the vector field with default value `v̄ = v`.

### Function Definitions

$(dae_functions)

With the above function definitions the `DAEProblem` can be created by

```julia
timespan = (0.0, 1.0)
timestep = 0.1
q₀ = [1., 1.]
λ₀ = [0.]
μ₀ = [0.]

prob = DAEProblem(v, u, ϕ, timespan, timestep, q₀, λ₀)
```
or
```julia
prob = DAEProblem(v, u, ϕ, ū, ψ, timespan, timestep, q₀, λ₀, μ₀)
```
"""
const DAEProblem = EquationProblem{DAE}

function DAEProblem(v, u, ϕ, ū, ψ, timespan::Tuple, timestep::Real, ics...; v̄ = v,
                    invariants = NullInvariants(),
                    parameters = NullParameters(),
                    periodicity = NullPeriodicity())
    equ = DAE(v, u, ϕ, ū, ψ, v̄, invariants, parameter_types(parameters), periodicity)
    EquationProblem(equ, timespan, timestep, initialstate(equ, ics...), parameters)
end

function DAEProblem(v, u, ϕ, args...; kwargs...)
    DAEProblem(v, u, ϕ, nothing, nothing, args...; kwargs...)
end

function GeometricBase.periodicity(prob::DAEProblem)
    (q = periodicity(equation(prob)), λ = NullPeriodicity(), μ = NullPeriodicity())
end

@inline GeometricBase.nconstraints(prob::DAEProblem) = length(initial_conditions(prob).λ)


const DAEEnsemble = EnsembleProblem{DAE}

function DAEEnsemble(v, u, ϕ, ū, ψ, timespan::Tuple, timestep::Real, ics...; v̄ = v,
        invariants = NullInvariants(),
        parameters = NullParameters(),
        periodicity = NullPeriodicity())
    equ = DAE(v, u, ϕ, ū, ψ, v̄)
    EnsembleProblem(equ, timespan, timestep, initialstate(equ, ics...), parameters)
end

function DAEEnsemble(v, u, ϕ, args...; kwargs...)
    DAEEnsemble(v, u, ϕ, nothing, nothing, args...; kwargs...)
end
