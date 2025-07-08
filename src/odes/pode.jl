
const pode_equations = raw"""
A partitioned ordinary differential equation is an initial value problem of the form
```math
\begin{aligned}
\dot{q} (t) &= v(t, q(t), p(t)) , \\
\dot{p} (t) &= f(t, q(t), p(t)) ,
\end{aligned}
```
with vector fields ``v`` and ``f``.
"""

const pode_functions = raw"""
The functions `v` and `f` must have the interface
```julia
function v(v, t, q, p, params)
    v[1] = ...
    v[2] = ...
    ...
end

function f(f, t, q, p, params)
    f[1] = ...
    f[2] = ...
    ...
end
```
where `t` is the current time, `q` and `p` are the current solution vectors,
`v` and `f` are the vectors which hold the result of evaluating the vector
fields ``v`` and ``f`` on `t`, `q` and `p`, and params is a `NamedTuple` of
additional parameters.
"""


@doc """
`PODE`: Partitioned Ordinary Differential Equation

$(pode_equations)

### Parameters

* `vType <: Callable`: type of `v`
* `fType <: Callable`: type of `f`
* `invType <: OptionalInvariants`: invariants type
* `parType <: OptionalParameters`: parameters type
* `perType <: OptionalPeriodicity`: periodicity type

### Fields

* `v`: function computing the vector field ``v``
* `f`: function computing the vector field ``f``
* `invariants`: functions for the computation of invariants, either a `NamedTuple` containing the equation's invariants or `NullInvariants`
* `parameters`: type constraints for parameters, either a `NamedTuple` containing the equation's parameters or `NullParameters`
* `periodicity`: determines the periodicity of the state vector `q` for cutting periodic solutions, either a `AbstractArray` or `NullPeriodicity`

### Constructors

```julia
PODE(v, f, invariants, parameters, periodicity)
PODE(v, f; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
```

### Function Definitions

$(pode_functions)

"""
struct PODE{vType <: Callable, fType <: Callable,
            v̄Type <: OptionalCallable,
            f̄Type <: OptionalCallable,
            invType <: OptionalInvariants,
            parType <: OptionalParameters,
            perType <: OptionalPeriodicity} <: AbstractEquationPODE{invType,parType,perType}

    v::vType
    f::fType
    v̄::v̄Type
    f̄::f̄Type

    invariants::invType
    parameters::parType
    periodicity::perType

    function PODE(v, f, v̄, f̄, invariants, parameters, periodicity)
        new{typeof(v), typeof(f), typeof(v̄), typeof(f̄), typeof(invariants), typeof(parameters), typeof(periodicity)}(
                v, f, v̄, f̄, invariants, parameters, periodicity)
    end
end

PODE(v, f; v̄ = v, f̄ = f, invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity()) = PODE(v, f, v̄, f̄, invariants, parameters, periodicity)

GeometricBase.invariants(equation::PODE) = equation.invariants
GeometricBase.parameters(equation::PODE) = equation.parameters
GeometricBase.periodicity(equation::PODE) = equation.periodicity

hasvectorfield(::PODE) = true
hasinitialguess(::PODE{vType, ftype, <:Callable, <:Callable}) where {vType, ftype} = true
hasinitialguess(::PODE{vType, ftype, <:Nothing, <:Nothing}) where {vType, ftype} = false

function Base.show(io::IO, equation::PODE)
    print(io, "Partitioned Ordinary Differential Equation (PODE)", "\n")
    print(io, "\n")
    print(io, " with vector fields", "\n")
    print(io, "   v = ", equation.v, "\n")
    print(io, "   f = ", equation.f, "\n")
    print(io, "\n")
    print(io, " Invariants: \n")
    print(io, "   ", invariants(equation))
end

function initialstate(::PODE, t::InitialTime, ics::NamedTuple, params::OptionalParameters)
    (
        q = _statevariable(ics.q),
        p = _statevariable(ics.p),
    )
end

function initialstate(equ::PODE, q₀::InitialState, p₀::InitialState)
    initialstate(equ, (q = q₀, p = p₀))
end

function initialstate(equ::PODE, q₀::InitialStateVector, p₀::InitialStateVector)
    [initialstate(equ, q, p) for (q,p) in zip(q₀,p₀)]
end

function check_initial_conditions(::PODE, ics::NamedTuple)
    haskey(ics, :q) || return false
    haskey(ics, :p) || return false
    eltype(ics.q) == eltype(ics.p) || return false
    typeof(ics.q) == typeof(ics.p) || return false
    axes(ics.q) == axes(ics.p) || return false
    typeof(ics.q) <: StateVariable || return false
    typeof(ics.p) <: StateVariable || return false
    return true
end

function check_methods(equ::PODE, timespan, ics, params)
    applicable(equ.v, vectorfield(ics.q), timespan[begin], ics.q, ics.p, params) || return false
    applicable(equ.f, vectorfield(ics.p), timespan[begin], ics.q, ics.p, params) || return false
    equ.v̄ === nothing || applicable(equ.v̄, vectorfield(ics.q), timespan[begin], ics.q, ics.p, params) || return false
    equ.f̄ === nothing || applicable(equ.f̄, vectorfield(ics.p), timespan[begin], ics.q, ics.p, params) || return false
    return true
end

function GeometricBase.datatype(equ::PODE, ics::NamedTuple)
    @assert check_initial_conditions(equ, ics)
    return eltype(ics.q)
end

function GeometricBase.arrtype(equ::PODE, ics::NamedTuple)
    @assert check_initial_conditions(equ, ics)
    typeof(ics.q)
end

_get_v(equ::PODE, params) = (v, t, q, p) -> equ.v(v, t, q, p, params)
_get_f(equ::PODE, params) = (f, t, q, p) -> equ.f(f, t, q, p, params)
_get_v̄(equ::PODE, params) = (v, t, q, p) -> equ.v̄(v, t, q, p, params)
_get_f̄(equ::PODE, params) = (f, t, q, p) -> equ.f̄(f, t, q, p, params)
_get_invariant(::PODE, inv, params) = (t,q,p) -> inv(t, q, p, params)

_functions(equ::PODE) = (v = equ.v, f = equ.f)
_functions(equ::PODE, params::OptionalParameters) = (v = _get_v(equ, params), f = _get_f(equ, params))
_initialguess(equ::PODE) = (v = equ.v̄, f = equ.f̄)
_initialguess(equ::PODE, params::OptionalParameters) = (v = _get_v̄(equ, params), f = _get_f̄(equ, params))


@doc """
`PODEProblem`: Partitioned Ordinary Differential Equation Problem

$(pode_equations)

The dynamical variables ``(q,p)`` with initial conditions ``(q(t_{0}) = q_{0}, p(t_{0}) = p_{0})``
take values in ``\\mathbb{R}^{d} \\times \\mathbb{R}^{d}``.

### Constructors

```julia
PODEProblem(v, f, timespan, timestep, ics; kwargs...)
PODEProblem(v, f, timespan, timestep, q₀::StateVariable, p₀::StateVariable; kwargs...)
PODEProblem(v, f, timespan, timestep, q₀::AbstractArray, p₀::AbstractArray; kwargs...)
```
where `v` and `f` are the function computing the vector fields, 
`timespan` is the time interval `(t₀,t₁)` for the problem to be solved in,
`timestep` is the time step to be used in the simulation, and
`ics` is a `NamedTuple` with entries `q` and `p`.
The initial conditions `q₀` and `p₀` can also be prescribed
directly, with `StateVariable` an `AbstractArray{<:Number}`.
For the interfaces of the functions `v` and `f` see [`PODE`](@ref).

For possible keyword arguments see the documentation on [`EquationProblem`](@ref GeometricEquations.EquationProblem) subtypes.

### Function Definitions

$(pode_functions)

"""
const PODEProblem = EquationProblem{PODE}

function PODEProblem(v, f, timespan, timestep, ics...;
        invariants = NullInvariants(),
        parameters = NullParameters(),
        periodicity = NullPeriodicity(),
        v̄ = v, f̄ = f)
    equ = PODE(v, f, v̄, f̄, invariants, parameter_types(parameters), periodicity)
    EquationProblem(equ, timespan, timestep, initialstate(equ, ics...), parameters)
end

function GeometricBase.periodicity(prob::PODEProblem)
    (q = periodicity(equation(prob)), p = NullPeriodicity())
end


@doc """
`PODEEnsemble`: Partitioned Ordinary Differential Equation Ensemble

$(pode_equations)

The dynamical variables ``(q,p)`` take values in ``\\mathbb{R}^{d} \\times \\mathbb{R}^{d}``.

### Constructors

```julia
PODEEnsemble(v, f, timespan, timestep, ics::AbstractVector{<: NamedTuple}; kwargs...)
PODEEnsemble(v, f, timespan, timestep, q₀::AbstractVector{<: StateVariable}, p₀::AbstractVector{<: StateVariable}; kwargs...)
PODEEnsemble(v, f, timespan, timestep, q₀::AbstractVector{<: AbstractArray}, p₀::AbstractVector{<: AbstractArray}; kwargs...)
```
where `v` and `f` are the function computing the vector fields, 
`timespan` is the time interval `(t₀,t₁)` for the problem to be solved in,
`timestep` is the time step to be used in the simulation, and
`ics` is an `AbstractVector` of `NamedTuple`, each with entries `q` and `p`.
The initial conditions `q₀` and `p₀` can also be prescribed, each as an
`AbstractVector` of `StateVariable` or `AbstractArray{<:Number}`.
For the interfaces of the functions `v` and `f` see [`PODE`](@ref).

For possible keyword arguments see the documentation on [`EnsembleProblem`](@ref GeometricEquations.EnsembleProblem) subtypes.

"""
const PODEEnsemble  = EnsembleProblem{PODE}

function PODEEnsemble(v, f, timespan, timestep, ics...;
        invariants = NullInvariants(),
        parameters = NullParameters(),
        periodicity = NullPeriodicity(),
        v̄ = v, f̄ = f)
    equ = PODE(v, f, v̄, f̄, invariants, parameter_types(parameters), periodicity)
    EnsembleProblem(equ, timespan, timestep, initialstate(equ, ics...), parameters)
end
