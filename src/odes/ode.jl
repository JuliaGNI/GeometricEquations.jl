
const ode_equations = raw"""
Ordinary differential equations define an initial value problem of the form
```math
\dot{q} (t) = v(t, q(t)) ,
```
with vector field ``v``.
"""

const ode_functions = raw"""
The function `v` providing the vector field must have the interface
```julia
function v(v, t, q, params)
    v[1] = ...
    v[2] = ...
    ...
end
```
where `t` is the current time, `q` is the current solution vector, `v` is the
vector which holds the result of evaluating the vector field ``v`` on `t` and
`q`, and `params` is a `NamedTuple` of additional parameters on which the
vector field may depend.
"""


@doc """
`ODE`: Ordinary Differential Equation

$(ode_equations)

### Parameters

* `vType <: Callable`: type of `v`
* `invType <: OptionalInvariants`: invariants type
* `parType <: OptionalParameters`: parameters type
* `perType <: OptionalPeriodicity`: periodicity type

### Fields

* `v`: function computing the vector field
* `invariants`: functions for the computation of invariants, either a `NamedTuple` containing the equation's invariants or `NullInvariants`
* `parameters`: type constraints for parameters, either a `NamedTuple` containing the equation's parameters or `NullParameters`
* `periodicity`: determines the periodicity of the state vector `q` for cutting periodic solutions, either a `AbstractArray` or `NullPeriodicity`

### Constructors

```julia
ODE(v, invariants, parameters, periodicity)
ODE(v; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
```

### Function Definitions

$(ode_functions)

"""
struct ODE{vType <: Callable,
           invType <: OptionalInvariants,
           parType <: OptionalParameters,
           perType <: OptionalPeriodicity} <: AbstractEquationODE{invType, parType, perType}
    v::vType

    invariants::invType
    parameters::parType
    periodicity::perType

    function ODE(v, invariants, parameters, periodicity)
        @assert !isempty(methods(v))
        # @assert hasmethod(v, (Real, AbstractArray, AbstractArray, OptionalParameters))

        new{typeof(v), typeof(invariants), typeof(parameters), typeof(periodicity)}(v,
                                                                                    invariants,
                                                                                    parameters,
                                                                                    periodicity)
    end
end

function ODE(v; invariants = NullInvariants(), parameters = NullParameters(),
             periodicity = NullPeriodicity())
    ODE(v, invariants, parameters, periodicity)
end

GeometricBase.invariants(equation::ODE) = equation.invariants
GeometricBase.parameters(equation::ODE) = equation.parameters
GeometricBase.periodicity(equation::ODE) = equation.periodicity

hasvectorfield(::ODE) = true

function Base.show(io::IO, equation::ODE)
    print(io, "Ordinary Differential Equation (ODE)", "\n")
    print(io, "\n")
    print(io, " with vector field")
    print(io, "\n")
    print(io, "   v = ", equation.v, "\n")
    print(io, "\n")
    print(io, " Invariants: \n")
    print(io, "   ", invariants(equation))
end

function initialstate(::ODE, t::InitialTime, ics::NamedTuple, params::OptionalParameters)
    (
        q = _statevariable(ics.q),
    )
end

function initialstate(equ::ODE, q₀::InitialState)
    initialstate(equ, (q = q₀,))
end

function initialstate(equ::ODE, q₀::InitialStateVector)
    [initialstate(equ, q) for q in q₀]
end

function check_initial_conditions(::ODE, ics::NamedTuple)
    haskey(ics, :q) || return false
    typeof(ics.q) <: StateVariable || return false
    return true
end

function check_methods(equ::ODE, tspan, ics, params)
    applicable(equ.v, vectorfield(ics.q), tspan[begin], ics.q, params) || return false
    return true
end

function GeometricBase.datatype(equ::ODE, ics::NamedTuple)
    @assert check_initial_conditions(equ, ics)
    return eltype(ics.q)
end

function GeometricBase.arrtype(equ::ODE, ics::NamedTuple)
    @assert check_initial_conditions(equ, ics)
    typeof(ics.q)
end

_get_v(equ::ODE, params) = (v, t, q) -> equ.v(v, t, q, params)
_get_v̄(equ::ODE, params) = _get_v(equ, params)
_get_invariant(::ODE, inv, params) = (t, q) -> inv(t, q, params)

_functions(equ::ODE) = (v = equ.v,)
_functions(equ::ODE, params::OptionalParameters) = (v = _get_v(equ, params),)


@doc """
`ODEProblem`: Ordinary Differential Equation Problem

$(ode_equations)

The dynamical variables ``q`` with initial condition ``q_{0}`` take values in ``\\mathbb{R}^{d}``.

### Constructors

```julia
ODEProblem(v, tspan, tstep, ics::NamedTuple; kwargs...)
ODEProblem(v, tspan, tstep, q₀::StateVariable; kwargs...)
ODEProblem(v, tspan, tstep, q₀::AbstractArray; kwargs...)
```
where `v` is the function computing the vector field, 
`tspan` is the time interval `(t₀,t₁)` for the problem to be solved in,
`tstep` is the time step to be used in the simulation, and
`ics` is a `NamedTuple` with entry `q` of type `StateVariable`.
The initial condition `q₀` can also be prescribed directly, as a
`StateVariable` or an `AbstractArray{<:Number}`.

For possible keyword arguments see the documentation on [`EquationProblem`](@ref GeometricEquations.EquationProblem) subtypes.

### Function Definitions

$(ode_functions)

"""
const ODEProblem = EquationProblem{ODE}

function ODEProblem(v, tspan, tstep, ics...;
                    invariants = NullInvariants(),
                    parameters = NullParameters(),
                    periodicity = NullPeriodicity())
    equ = ODE(v, invariants, parameter_types(parameters), periodicity)
    EquationProblem(equ, tspan, tstep, initialstate(equ, ics...), parameters)
end

GeometricBase.periodicity(prob::ODEProblem) = (q = periodicity(equation(prob)),)


@doc """
`ODEEnsemble`: Ordinary Differential Equation Ensemble

$(ode_equations)

The dynamical variables ``q`` take values in ``\\mathbb{R}^{d}``.

### Constructors

```julia
ODEEnsemble(v, tspan, tstep, ics::AbstractVector{<: NamedTuple}; kwargs...)
ODEEnsemble(v, tspan, tstep, q₀::AbstractVector{<: StateVariable}; kwargs...)
ODEEnsemble(v, tspan, tstep, q₀::AbstractVector{<: AbstractArray}; kwargs...)
```
where `v` is the function computing the vector field, 
`tspan` is the time interval `(t₀,t₁)` for the problem to be solved in,
`tstep` is the time step to be used in the simulation, and 
`ics` is an `AbstractVector` of `NamedTuple`, each with an entry `q`
of type `StateVariable`.
The initial condition `q₀` can also be prescribed as an `AbstractVector`
of `StateVariable` or `AbstractArray{<:Number}`.
For the interface of the function `v` see [`ODE`](@ref).

For possible keyword arguments see the documentation on [`EnsembleProblem`](@ref GeometricEquations.EnsembleProblem) subtypes.

"""
const ODEEnsemble = EnsembleProblem{ODE}

function ODEEnsemble(v, tspan, tstep, ics...;
                    invariants = NullInvariants(),
                    parameters = NullParameters(),
                    periodicity = NullPeriodicity())
    equ = ODE(v, invariants, parameter_types(parameters), periodicity)
    EnsembleProblem(equ, tspan, tstep, initialstate(equ, ics...), parameters)
end
