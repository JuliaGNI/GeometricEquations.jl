
const dele_equations = raw"""
Discrete Euler-Lagrange equations define an initial value problem of the form
```math
D_1 L_d (q_{n}, q_{n+1}) + D_2 L_d (q_{n-1}, q_{n}) = 0 ,
```
with discrete Lagrangian ``L_d``.
"""

const dele_functions = raw"""
The function `Ld` providing the discrete Lagrangian must have the interface
```julia
function Ld(t_{n}, t_{n+1}, q_{n}, q_{n+1}, params)
    return ...
end
```
where `t_{n}` and `t_{n+1}` are the time of the previous and next step, `q_{n}`
and `q_{n+1}` are the solution vectors of the previous and next step, and
`params` is a `NamedTuple` of additional parameters on which the Lagrangian may
depend.
The derivatives of the discrete Lagrangian with respect to its first and second
argument, `D_1 L_d` and `D_2 L_d`, respectively, must have the interfaces
```julia
function D1Ld(D, t_{n}, t_{n+1}, q_{n}, q_{n+1}, params)
    D1[1] = ...
    D1[2] = ...
    ...
end
```
and
```julia
function D2Ld(D2, t_{n}, t_{n+1}, q_{n}, q_{n+1}, params)
    D2[1] = ...
    D2[2] = ...
    ...
end
```
"""


@doc """
`DELE`: Discrete Euler-Lagrange Equation

$(dele_equations)

### Parameters

* `LdType <: Callable`: type of `Ld`
* `D1LdType <: Callable`: type of `D1Ld`
* `D2LdType <: Callable`: type of `D2Ld`
* `invType <: OptionalInvariants`: invariants type
* `parType <: OptionalParameters`: parameters type
* `perType <: OptionalPeriodicity`: periodicity type

### Fields

* `Ld`: function computing the discrete Lagrangian
* `D1Ld`: function computing the derivative of the discrete Lagrangian w.r.t. the first argument
* `D2Ld`: function computing the derivative of the discrete Lagrangian w.r.t. the second argument
* `invariants`: functions for the computation of invariants, either a `NamedTuple` containing the equation's invariants or `NullInvariants`
* `parameters`: type constraints for parameters, either a `NamedTuple` containing the equation's parameters or `NullParameters`
* `periodicity`: determines the periodicity of the state vector `q` for cutting periodic solutions, either a `AbstractArray` or `NullPeriodicity`

### Constructors

```julia
DELE(Ld, D1Ld, D2Ld, invariants, parameters, periodicity)
DELE(Ld, D1Ld, D2Ld; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
```

### Function Definitions

$(dele_functions)

"""
struct DELE{LdType <: Callable, D1LdType <: Callable, D2LdType <: Callable,
           invType <: OptionalInvariants,
           parType <: OptionalParameters,
           perType <: OptionalPeriodicity} <: DiscreteEquation{invType, parType, perType}

    Ld::LdType
    D1Ld::D1LdType
    D2Ld::D2LdType

    invariants::invType
    parameters::parType
    periodicity::perType

    function DELE(Ld, D1Ld, D2Ld, invariants, parameters, periodicity)
        @assert !isempty(methods(Ld))
        @assert !isempty(methods(D1Ld))
        @assert !isempty(methods(D2Ld))
        # @assert hasmethod(v, (Real, AbstractArray, AbstractArray, OptionalParameters))

        _periodicity = promote_periodicity(periodicity)

        new{typeof(Ld), typeof(D1Ld), typeof(D2Ld),
            typeof(invariants),
            typeof(parameters),
            typeof(_periodicity)
            }(Ld, D1Ld, D2Ld, invariants, parameters, _periodicity)
    end
end

function DELE(Ld, D1Ld, D2Ld;
             invariants = NullInvariants(),
             parameters = NullParameters(),
             periodicity = NullPeriodicity())
    DELE(Ld, D1Ld, D2Ld, invariants, parameters, periodicity)
end

GeometricBase.invariants(equation::DELE) = equation.invariants
GeometricBase.parameters(equation::DELE) = equation.parameters
GeometricBase.periodicity(equation::DELE) = equation.periodicity

hasvectorfield(::DELE) = true
hasinitialguess(::DELE) = false
# hasinitialguess(::DELE{LdType, D1LdType, D2LdType, <:Callable}) where {LdType, D1LdType, D2LdType} = true
# hasinitialguess(::DELE{LdType, D1LdType, D2LdType, <:Nothing}) where {LdType, D1LdType, D2LdType} = false

function Base.show(io::IO, equation::DELE)
    print(io, "Discrete Euler-Lagrange Equation (DELE)", "\n")
    print(io, "\n")
    print(io, " with discrete Lagrangian")
    print(io, "\n")
    print(io, "   Ld = ", equation.Ld, "\n")
    print(io, " and its derivatives")
    print(io, "\n")
    print(io, "   D1Ld = ", equation.D1Ld, "\n")
    print(io, "   D2Ld = ", equation.D2Ld, "\n")
    print(io, "\n")
    print(io, " Invariants: \n")
    print(io, "   ", invariants(equation))
end

function initialstate(equ::DELE, t::InitialTime, ics::NamedTuple, params::OptionalParameters)
    (
        q̄ = _statevariable(ics.q̄, periodicity(equ)),
        q = _statevariable(ics.q, periodicity(equ)),
    )
end

function initialstate(equ::DELE, q₀::InitialState, q₁::InitialState)
    initialstate(equ, (q̄ = q₀, q = q₁))
end

function initialstate(equ::DELE, q₀::InitialStateVector, q₁::InitialStateVector)
    [initialstate(equ, q̄, q) for (q̄,q) in zip(q₀,q₁)]
end

function check_initial_conditions(::DELE, ics::NamedTuple)
    haskey(ics, :q̄) || return false
    haskey(ics, :q) || return false
    typeof(ics.q̄) <: StateVariable || return false
    typeof(ics.q) <: StateVariable || return false
    return true
end

function check_methods(equ::DELE, tspan, ics, params)
    applicable(equ.Ld, tspan[begin], tspan[end], ics.q̄, ics.q, params) || return false
    applicable(equ.D1Ld, vectorfield(ics.q), tspan[begin], tspan[end], ics.q̄, ics.q, params) || return false
    applicable(equ.D2Ld, vectorfield(ics.q), tspan[begin], tspan[end], ics.q̄, ics.q, params) || return false
    # equ.v̄ === nothing || applicable(equ.v̄, vectorfield(ics.q), tspan[begin], ics.q, params) || return false
    return true
end

function GeometricBase.datatype(equ::DELE, ics::NamedTuple)
    @assert check_initial_conditions(equ, ics)
    return eltype(ics.q)
end

function GeometricBase.arrtype(equ::DELE, ics::NamedTuple)
    @assert check_initial_conditions(equ, ics)
    typeof(ics.q)
end

_functions(equ::DELE) = (Ld = equ.Ld, D1Ld = equ.D1Ld, D2Ld = equ.D2Ld)
# _initialguess(equ::DELE) = (v = equ.v̄,)


@doc """
`DELEProblem`: Discrete Euler-Lagrange Equation Problem

$(dele_equations)

The dynamical variables ``q`` with initial conditions ``q_{0}`` and ``q_{1}`` take values in ``\\mathbb{R}^{d}``.

### Constructors

```julia
DELEProblem(Ld, D1Ld, D2Ld, tspan, tstep, ics::NamedTuple; kwargs...)
DELEProblem(Ld, D1Ld, D2Ld, tspan, tstep, q₀::StateVariable, q₁::StateVariable; kwargs...)
DELEProblem(Ld, D1Ld, D2Ld, tspan, tstep, q₀::AbstractArray, q₁::AbstractArray; kwargs...)
```
where `Ld` is the function computing the discrete Lagrangian, `D1Ld` and `D2Ld`
are functions computing the derivative of `Ld` with respect to the first and second argument,
`tspan` is the time interval `(t₀,tₙ)` for the problem to be solved in,
`tstep` is the time step to be used in the simulation, and
`ics` is a `NamedTuple` with entries `q₀` and  `q₁` of type `StateVariable`.
The initial conditions `q₀` and  `q₁` can also be prescribed directly, as a
`StateVariable` or an `AbstractArray{<:Number}`.

For the interface of the functions `Ld`, `D1Ld` and `D2Ld` see [`DELE`](@ref).

For possible keyword arguments see the documentation on [`EquationProblem`](@ref GeometricEquations.EquationProblem) subtypes.

### Function Definitions

$(dele_functions)

"""
const DELEProblem = EquationProblem{DELE}

function DELEProblem(Ld, D1Ld, D2Ld, tspan, tstep, ics...;
                    invariants = NullInvariants(),
                    parameters = NullParameters(),
                    periodicity = NullPeriodicity())
    equ = DELE(Ld, D1Ld, D2Ld, invariants, parameter_types(parameters), periodicity)
    EquationProblem(equ, tspan, tstep, initialstate(equ, ics...), parameters)
end

GeometricBase.periodicity(prob::DELEProblem) = (q = periodicity(equation(prob)),)


@doc """
`DELEEnsemble`: Discrete Euler-Lagrange Equation Ensemble

$(dele_equations)

The dynamical variables ``q`` take values in ``\\mathbb{R}^{d}``.

### Constructors

```julia
DELEEnsemble(Ld, D1Ld, D2Ld, tspan, tstep, ics::AbstractVector{<: NamedTuple}; kwargs...)
DELEEnsemble(Ld, D1Ld, D2Ld, tspan, tstep, q₀::AbstractVector{<: StateVariable}, q₁::AbstractVector{<: StateVariable}; kwargs...)
DELEEnsemble(Ld, D1Ld, D2Ld, tspan, tstep, q₀::AbstractVector{<: AbstractArray}, q₁::AbstractVector{<: AbstractArray}; kwargs...)
```
where `Ld` is the function computing the discrete Lagrangian, `D1Ld` and `D2Ld`
are functions computing the derivative of `Ld` with respect to the first and second argument,
`tspan` is the time interval `(t₀,t₁)` for the problem to be solved in,
`tstep` is the time step to be used in the simulation, and
`ics` is an `AbstractVector` of `NamedTuple`, each with entries `q₀` and  `q₁` of type `StateVariable`.
The initial conditions `q₀` and  `q₁` can also be prescribed as two `AbstractVector`s
of `StateVariable` or `AbstractArray{<:Number}`.
For the interface of the functions `Ld`, `D1Ld` and `D2Ld` see [`DELE`](@ref).

For possible keyword arguments see the documentation on [`EnsembleProblem`](@ref GeometricEquations.EnsembleProblem) subtypes.

"""
const DELEEnsemble = EnsembleProblem{DELE}

function DELEEnsemble(Ld, D1Ld, D2Ld, tspan, tstep, ics...;
                    invariants = NullInvariants(),
                    parameters = NullParameters(),
                    periodicity = NullPeriodicity())
    equ = DELE(Ld, D1Ld, D2Ld, invariants, parameter_types(parameters), periodicity)
    EnsembleProblem(equ, tspan, tstep, initialstate(equ, ics...), parameters)
end
