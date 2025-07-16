
const hode_equations = raw"""
A canonical Hamiltonian system of equations is special case of a
partitioned ordinary differential equation,
```math
\begin{aligned}
\dot{q} (t) &= v(t, q(t), p(t)) , \\
\dot{p} (t) &= f(t, q(t), p(t)) ,
\end{aligned}
```
with vector fields ``v`` and ``f``, given by
```math
\begin{aligned}
v &=   \frac{\partial H}{\partial p} , &
f &= - \frac{\partial H}{\partial q} .
\end{aligned}
```
"""

const hode_functions = raw"""
The functions `v`, `f` and `hamiltonian` must have the interface
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

function hamiltonian(t, q, p, params)
    return ...
end
```
where `t` is the current time, `q` and `p` are the current solution vectors,
`v` and `f` are the vectors which hold the result of evaluating the vector
fields on `t`, `q` and `p`, and params is a `NamedTuple` of
additional parameters.
"""

@doc """
`HODE`: Hamiltonian Ordinary Differential Equation

$(hode_equations)

### Parameters

* `vType <: Callable`: type of `v`
* `fType <: Callable`: type of `f`
* `hamType <: Callable`: Hamiltonian type
* `invType <: OptionalInvariants`: invariants type
* `parType <: OptionalParameters`: parameters type
* `perType <: OptionalPeriodicity`: periodicity type

### Fields

* `v`: function computing the vector field ``v``
* `f`: function computing the vector field ``f``
* `hamiltonian`: function computing the Hamiltonian ``H`` (usually the total energy of the system)
* `invariants`: functions for the computation of invariants, either a `NamedTuple` containing the equation's invariants or `NullInvariants`
* `parameters`: type constraints for parameters, either a `NamedTuple` containing the equation's parameters or `NullParameters`
* `periodicity`: determines the periodicity of the state vector `q` for cutting periodic solutions, either a `AbstractArray` or `NullPeriodicity`

### Constructors

```julia
HODE(v, f, hamiltonian, invariants, parameters, periodicity)
HODE(v, f, hamiltonian; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
```

### Function Definitions

$(hode_functions)

"""
struct HODE{vType <: Callable, fType <: Callable,
    v̄Type <: OptionalCallable,
    f̄Type <: OptionalCallable,
    hamType <: Callable,
    invType <: OptionalInvariants,
    parType <: OptionalParameters,
    perType <: OptionalPeriodicity} <:
       AbstractEquationPODE{invType, parType, perType}
    v::vType
    f::fType
    v̄::v̄Type
    f̄::f̄Type

    hamiltonian::hamType
    invariants::invType
    parameters::parType
    periodicity::perType

    function HODE(v, f, v̄, f̄, hamiltonian, invariants, parameters, periodicity)
        _periodicity = promote_periodicity(periodicity)
        new{typeof(v), typeof(f), typeof(v̄), typeof(f̄), typeof(hamiltonian),
            typeof(invariants), typeof(parameters), typeof(_periodicity)}(
            v, f, v̄, f̄, hamiltonian, invariants, parameters, _periodicity)
    end
end

function HODE(v, f, hamiltonian; v̄ = v, f̄ = f, invariants = NullInvariants(),
        parameters = NullParameters(), periodicity = NullPeriodicity())
    HODE(v, f, v̄, f̄, hamiltonian, invariants, parameters, periodicity)
end

GeometricBase.invariants(equation::HODE) = equation.invariants
GeometricBase.parameters(equation::HODE) = equation.parameters
GeometricBase.periodicity(equation::HODE) = equation.periodicity

hasvectorfield(::HODE) = true
hashamiltonian(::HODE) = true
hasinitialguess(::HODE{vType, ftype, <:Callable, <:Callable}) where {vType, ftype} = true
hasinitialguess(::HODE{vType, ftype, <:Nothing, <:Nothing}) where {vType, ftype} = false

function Base.show(io::IO, equation::HODE)
    print(io, "Hamiltonian Ordinary Differential Equation (HODE)", "\n")
    print(io, "\n")
    print(io, " with vector fields", "\n")
    print(io, "   v = ", equation.v, "\n")
    print(io, "   f = ", equation.f, "\n")
    print(io, "\n")
    print(io, " Hamiltonian: H = ", equation.hamiltonian, "\n")
    print(io, "\n")
    print(io, " Invariants: \n")
    print(io, "   ", invariants(equation))
end

function initialstate(
        equ::HODE, t::InitialTime, ics::NamedTuple, params::OptionalParameters)
    (
        q = _statevariable(ics.q, periodicity(equ)),
        p = _statevariable(ics.p, NullPeriodicity())
    )
end

function initialstate(equ::HODE, q₀::InitialState, p₀::InitialState)
    initialstate(equ, (q = q₀, p = p₀))
end

function initialstate(equ::HODE, q₀::InitialStateVector, p₀::InitialStateVector)
    [initialstate(equ, q, p) for (q, p) in zip(q₀, p₀)]
end

function check_initial_conditions(::HODE, ics::NamedTuple)
    haskey(ics, :q) || return false
    haskey(ics, :p) || return false
    eltype(ics.q) == eltype(ics.p) || return false
    typeof(ics.q) == typeof(ics.p) || return false
    axes(ics.q) == axes(ics.p) || return false
    typeof(ics.q) <: StateVariable || return false
    typeof(ics.p) <: StateVariable || return false
    return true
end

function check_methods(equ::HODE, timespan, ics, params)
    applicable(equ.v, vectorfield(ics.q), timespan[begin], ics.q, ics.p, params) ||
        return false
    applicable(equ.f, vectorfield(ics.p), timespan[begin], ics.q, ics.p, params) ||
        return false
    applicable(equ.hamiltonian, timespan[begin], ics.q, ics.p, params) || return false
    equ.v̄ === nothing ||
        applicable(equ.v̄, vectorfield(ics.q), timespan[begin], ics.q, ics.p, params) ||
        return false
    equ.f̄ === nothing ||
        applicable(equ.f̄, vectorfield(ics.p), timespan[begin], ics.q, ics.p, params) ||
        return false
    return true
end

function GeometricBase.datatype(equ::HODE, ics::NamedTuple)
    @assert check_initial_conditions(equ, ics)
    return eltype(ics.q)
end

function GeometricBase.arrtype(equ::HODE, ics::NamedTuple)
    @assert check_initial_conditions(equ, ics)
    typeof(ics.q)
end

_get_v(equ::HODE, params) = (v, t, q, p) -> equ.v(v, t, q, p, params)
_get_f(equ::HODE, params) = (f, t, q, p) -> equ.f(f, t, q, p, params)
_get_v̄(equ::HODE, params) = (v, t, q, p) -> equ.v̄(v, t, q, p, params)
_get_f̄(equ::HODE, params) = (f, t, q, p) -> equ.f̄(f, t, q, p, params)
_get_h(equ::HODE, params) = (t, q, p) -> equ.hamiltonian(t, q, p, params)
_get_invariant(::HODE, inv, params) = (t, q, p) -> inv(t, q, p, params)

_functions(equ::HODE) = (v = equ.v, f = equ.f, h = equ.hamiltonian)
function _functions(equ::HODE, params::OptionalParameters)
    (v = _get_v(equ, params), f = _get_f(equ, params), h = _get_h(equ, params))
end
_initialguess(equ::HODE) = (v = equ.v̄, f = equ.f̄)
function _initialguess(equ::HODE, params::OptionalParameters)
    (v = _get_v̄(equ, params), f = _get_f̄(equ, params))
end

@doc """
`HODEProblem`: Hamiltonian Ordinary Differential Equation Problem

$(hode_equations)

The dynamical variables ``(q,p)`` with initial conditions ``(q(t_{0}) = q_{0}, p(t_{0}) = p_{0})``
take values in ``T^{*} Q \\simeq \\mathbb{R}^{d} \\times \\mathbb{R}^{d}``.

### Constructors

```julia
HODEProblem(v, f, hamiltonian, timespan, timestep, ics; kwargs...)
HODEProblem(v, f, hamiltonian, timespan, timestep, q₀::StateVariable, p₀::StateVariable; kwargs...)
HODEProblem(v, f, hamiltonian, timespan, timestep, q₀::AbstractArray, p₀::AbstractArray; kwargs...)
```
where `v` and `f` are the function computing the vector fields,
`hamiltonian` returns the value of the Hamiltonian (i.e. the total energy),
`timespan` is the time interval `(t₀,t₁)` for the problem to be solved in,
`timestep` is the time step to be used in the simulation, and
`ics` is a `NamedTuple` with entries `q` and `p`.
The initial conditions `q₀` and `p₀` can also be prescribed
directly, with `StateVariable` an `AbstractArray{<:Number}`.
For the interfaces of the functions `v`, `f` and `hamiltonian` see [`HODE`](@ref).

For possible keyword arguments see the documentation on [`EquationProblem`](@ref GeometricEquations.EquationProblem) subtypes.

### Function Definitions

$(hode_functions)

"""
const HODEProblem = EquationProblem{HODE}

function HODEProblem(v, f, hamiltonian, timespan, timestep, ics...;
        invariants = NullInvariants(),
        parameters = NullParameters(),
        periodicity = NullPeriodicity(),
        v̄ = v, f̄ = f)
    equ = HODE(
        v, f, v̄, f̄, hamiltonian, invariants, parameter_types(parameters), periodicity)
    EquationProblem(equ, timespan, timestep, initialstate(equ, ics...), parameters)
end

function GeometricBase.periodicity(prob::HODEProblem)
    (q = periodicity(equation(prob)), p = NullPeriodicity())
end

function compute_vectorfields!(vecfield, sol, prob::HODEProblem)
    initialguess(prob).v(vecfield.q, sol.t, sol.q, sol.p, parameters(prob))
    initialguess(prob).f(vecfield.p, sol.t, sol.q, sol.p, parameters(prob))
end

@doc """
`HODEEnsemble`: Hamiltonian Ordinary Differential Equation Ensemble

$(hode_equations)

The dynamical variables ``(q,p)`` take values in ``T^{*} Q \\simeq \\mathbb{R}^{d} \\times \\mathbb{R}^{d}``.

### Constructors

```julia
HODEEnsemble(v, f, hamiltonian, timespan, timestep, ics::AbstractVector{<: NamedTuple}; kwargs...)
HODEEnsemble(v, f, hamiltonian, timespan, timestep, q₀::AbstractVector{<: StateVariable}, p₀::AbstractVector{<: StateVariable}; kwargs...)
HODEEnsemble(v, f, hamiltonian, timespan, timestep, q₀::AbstractVector{<: AbstractArray}, p₀::AbstractVector{<: AbstractArray}; kwargs...)
```
where `v` and `f` are the function computing the vector fields,
`hamiltonian` returns the value of the Hamiltonian (i.e. the total energy),
`timespan` is the time interval `(t₀,t₁)` for the problem to be solved in,
`timestep` is the time step to be used in the simulation, and
`ics` is an `AbstractVector` of `NamedTuple`, each with entries `q` and `p`.
The initial conditions `q₀` and `p₀` can also be prescribed, each as an
`AbstractVector` of `StateVariable` or `AbstractArray{<:Number}`.
For the interfaces of the functions `v`, `f` and `hamiltonian` see [`HODE`](@ref).

For possible keyword arguments see the documentation on [`EnsembleProblem`](@ref GeometricEquations.EnsembleProblem) subtypes.

"""
const HODEEnsemble = EnsembleProblem{HODE}

function HODEEnsemble(v, f, hamiltonian, timespan, timestep, ics...;
        invariants = NullInvariants(),
        parameters = NullParameters(),
        periodicity = NullPeriodicity(),
        v̄ = v, f̄ = f)
    equ = HODE(
        v, f, v̄, f̄, hamiltonian, invariants, parameter_types(parameters), periodicity)
    EnsembleProblem(equ, timespan, timestep, initialstate(equ, ics...), parameters)
end
