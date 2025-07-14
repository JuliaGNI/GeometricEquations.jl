
const sode_equations = raw"""
Defines an initial value problem
```math
\dot{q} (t) = v(t, q(t)) , \qquad q(t_{0}) = q_{0} ,
```
with vector field ``v``, initial condition ``q_{0}`` and the solution
``q`` taking values in ``\mathbb{R}^{d}``. Here, the vector field ``v``
is given as a sum of vector fields
```math
v (t) = v_1 (t) + ... + v_r (t) .
```
"""

const sode_functions = raw"""
The functions `v_i` providing the vector field must have the interface
```julia
function v_i(v, t, q, params)
    v[1] = ...
    v[2] = ...
    ...
end
```
and the functions `q_i` providing the solutions must have the interface
```julia
function q_i(q₁, t₁, q₀, t₀, params)
    q₁[1] = q₀[1] + ...
    q₁[2] = q₀[2] + ...
    ...
end
```
where `t₀` is the current time, `q₀` is the current solution vector, `q₁` is the
new solution vector at time `t₁`, holding the result of computing one substep

The fact that the function `v` returns the solution and not just the vector
field for each substep increases the flexibility for the use of splitting
methods, e.g., it allows to use another integrator for solving substeps.
with the vector field ``v_i``.
"""

sode_equations_compatibility(v::Nothing, q::Nothing) = false
sode_equations_compatibility(v::Tuple, q::Nothing) = true
sode_equations_compatibility(v::Nothing, q::Tuple) = true
sode_equations_compatibility(v::Tuple, q::Tuple) = length(q) == length(v)

@doc """
`SODE`: Split Ordinary Differential Equation

$(sode_equations)

The dynamical variables ``q`` with initial condition ``q_{0}`` take values in ``\\mathbb{R}^{d}``.

### Parameters

* `vType <: Union{Tuple,Nothing}`: type of `v`
* `qType <: Union{Tuple,Nothing}`: type of `q`
* `invType <: OptionalInvariants`: invariants type
* `parType <: OptionalParameters`: parameters type
* `perType <: OptionalPeriodicity`: periodicity type

### Fields

* `v`: tuple of functions computing the vector fields for each substep
* `q`: tuple of functions computing the solutions for each substep
* `invariants`: functions for the computation of invariants, either a `NamedTuple` containing the equation's invariants or `NullInvariants`
* `parameters`: type constraints for parameters, either a `NamedTuple` containing the equation's parameters or `NullParameters`
* `periodicity`: determines the periodicity of the state vector `q` for cutting periodic solutions, either a `AbstractArray` or `NullPeriodicity`

### Constructors

```julia
SODE(v, invariants, parameters, periodicity)
SODE(v; invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity())
SODE(v, q, invariants, parameters, periodicity)
SODE(v, q; invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity())
```

### Function Definitions

$(sode_functions)

"""
struct SODE{vType <: Union{Tuple,Nothing}, qType <: Union{Tuple,Nothing}, v̄Type <: OptionalCallable,
            invType <: OptionalInvariants,
            parType <: OptionalParameters,
            perType <: OptionalPeriodicity} <: AbstractEquationODE{invType,parType,perType}

    v::vType
    q::qType
    v̄::v̄Type

    invariants::invType
    parameters::parType
    periodicity::perType

    function SODE(v, q, v̄, invariants, parameters, periodicity)
        @assert sode_equations_compatibility(v, q)
        _periodicity = promote_periodicity(periodicity)
        new{typeof(v), typeof(q), typeof(v̄), typeof(invariants), typeof(parameters), typeof(_periodicity)}(
                v, q, v̄, invariants, parameters, _periodicity)
    end
end

_sode_default_v̄(v, t, q, params) = nothing

SODE(v, q=nothing; v̄ = _sode_default_v̄, invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity()) = SODE(v, q, v̄, invariants, parameters, periodicity)

GeometricBase.invariants(equation::SODE) = equation.invariants
GeometricBase.parameters(equation::SODE) = equation.parameters
GeometricBase.periodicity(equation::SODE) = equation.periodicity

GeometricBase.nsteps(equ::SODE{<:Nothing}) = length(equ.q)
GeometricBase.nsteps(equ::SODE) = length(equ.v)

const SODEQT{QT,VT,invT,parT,perT} = SODE{VT,QT,invT,parT,perT} # type alias for dispatch on solution type parameter
const SODEVT{VT,QT,invT,parT,perT} = SODE{VT,QT,invT,parT,perT} # type alias for dispatch on vector field type parameter

hassolution(::SODEQT{<:Nothing}) = false
hassolution(::SODEQT{<:Tuple}) = true # && all(typeof(Q) <: Functiong for Q in equ.q)

hassolution(::SODEQT{<:Nothing}, i::Int) = false
hassolution(equ::SODEQT{<:Tuple}, i::Int) = i ≤ length(equ.q)# && typeof(equ.q[i]) <: Function

hasvectorfield(::SODEVT{<:Nothing}) = false
hasvectorfield(::SODEVT{<:Tuple}) = true # && all(typeof(V) <: Function for V in equ.v)

hasvectorfield(::SODEVT{<:Nothing}, i::Int) = false
hasvectorfield(equ::SODEVT{<:Tuple}, i::Int) = i ≤ length(equ.v)# && typeof(equ.v[i]) <: Function

hasinitialguess(::SODE{vType, qType, <:Callable}) where {vType,qType} = true
hasinitialguess(::SODE{vType, qType, <:Nothing}) where {vType,qType} = false

function Base.show(io::IO, equation::SODE)
    print(io, "Split Ordinary Differential Equation (SODE)", "\n")
    print(io, "\n")
    print(io, " with vector fields")
    print(io, "\n")
    print(io, "   v = ", equation.v, "\n")
    print(io, "\n")
    print(io, " and solutions")
    print(io, "\n")
    print(io, "   q = ", equation.q, "\n")
    print(io, "\n")
    print(io, " Invariants: \n")
    print(io, "   ", invariants(equation))
end

function initialstate(equ::SODE, t::InitialTime, ics::NamedTuple, params::OptionalParameters)
    (
        q = _statevariable(ics.q, periodicity(equ)),
    )
end

function initialstate(equ::SODE, q₀::InitialState)
    initialstate(equ, (q = q₀,))
end

function initialstate(equ::SODE, q₀::InitialStateVector)
    [initialstate(equ, q) for q in q₀]
end

function check_initial_conditions(::SODE, ics::NamedTuple)
    haskey(ics, :q) || return false
    typeof(ics.q) <: StateVariable || return false
    return true
end

function check_methods(equ::SODE, timespan, ics, params)
    if hasvectorfield(equ)
        for v in equ.v
            applicable(v, vectorfield(ics.q), timespan[begin], ics.q, params) || return false
        end
    end
    if hassolution(equ)
        for q in equ.q
            applicable(q, vectorfield(ics.q), timespan[end], ics.q, timespan[begin], params) || return false
        end
    end
    applicable(equ.v̄, vectorfield(ics.q), timespan[begin], ics.q, params) || return false
    return true
end

function GeometricBase.datatype(equ::SODE, ics::NamedTuple)
    @assert check_initial_conditions(equ, ics)
    return eltype(ics.q)
end

function GeometricBase.arrtype(equ::SODE, ics::NamedTuple)
    @assert check_initial_conditions(equ, ics)
    typeof(ics.q)
end

_get_v(equ::SODE, params) = Tuple((v, t, q) -> V(v, t, q, params) for V in equ.v)
_get_q(equ::SODE, params) = Tuple((q̄, t̄, q, t) -> Q(q̄, t̄, q, t, params) for Q in equ.q)
_get_v̄(equ::SODE, params) = (v, t, q) -> equ.v̄(v, t, q, params)
_get_v(::SODEVT{<:Nothing}, params) = nothing
_get_q(::SODEQT{<:Nothing}, params) = nothing
_get_invariant(::SODE, inv, params) = (t, q) -> inv(t, q, params)

_functions(equ::SODE) = (v = equ.v,)
_solutions(equ::SODE) = (q = equ.q,)
_functions(equ::SODE, params::OptionalParameters) = (v = _get_v(equ, params),)
_solutions(equ::SODE, params::OptionalParameters) = (q = _get_q(equ, params),)

_initialguess(equ::SODE) = (v = equ.v̄,)
_initialguess(equ::SODE, params::OptionalParameters) = (v = _get_v̄(equ, params),)


@doc """
`SODEProblem`: Split Ordinary Differential Equation Problem

$(sode_equations)

The dynamical variables ``q`` with initial condition ``q_{0}`` take values in ``\\mathbb{R}^{d}``.

### Constructors

```julia
SODEProblem(v, q, timespan, timestep, ics::NamedTuple; kwargs...)
SODEProblem(v, q, timespan, timestep, q₀::StateVariable; kwargs...)
SODEProblem(v, q, timespan, timestep, q₀::AbstractArray; kwargs...)
SODEProblem(v, timespan, timestep, ics::NamedTuple; kwargs...)
SODEProblem(v, timespan, timestep, q₀::StateVariable; kwargs...)
SODEProblem(v, timespan, timestep, q₀::AbstractArray; kwargs...)
```

where `v` is a tuple of functions computing the vector fields for each substep,
`q` is an optional tuple of functions computing the solution for each substep,
`timespan` is the time interval `(t₀,t₁)` for the problem to be solved in,
`timestep` is the time step to be used in the simulation, and
`ics` is a `NamedTuple` with entry `q`.
The initial condition `q₀` can also be prescribed directly, with
`StateVariable` an `AbstractArray{<:Number}`.

For possible keyword arguments see the documentation on [`EquationProblem`](@ref GeometricEquations.EquationProblem) subtypes.

### Function Definitions

$(sode_functions)

"""
const SODEProblem = EquationProblem{SODE}

function SODEProblem(v::Tuple, q::Union{Tuple, Nothing}, timespan::Tuple, timestep::Real, ics...;
        invariants = NullInvariants(),
        parameters = NullParameters(),
        periodicity = NullPeriodicity(),
        v̄ = _sode_default_v̄)
    equ = SODE(v, q, v̄, invariants, parameter_types(parameters), periodicity)
    EquationProblem(equ, timespan, timestep, initialstate(equ, ics...), parameters)
end

function SODEProblem(v, args...; kwargs...)
    SODEProblem(v, nothing, args...; kwargs...)
end

GeometricBase.nsteps(prob::SODEProblem) = nsteps(equation(prob))
GeometricBase.periodicity(prob::SODEProblem) = (q = periodicity(equation(prob)),)


@doc """
`SODEEnsemble`: Split Ordinary Differential Equation Ensemble

$(sode_equations)

The dynamical variables ``q`` take values in ``\\mathbb{R}^{d}``.

### Constructors

```julia
SODEEnsemble(v, q, timespan, timestep, ics::AbstractVector{<: NamedTuple}; kwargs...)
SODEEnsemble(v, q, timespan, timestep, q₀::AbstractVector{<: StateVariable}; kwargs...)
SODEEnsemble(v, q, timespan, timestep, q₀::AbstractVector{<: AbstractArray}; kwargs...)
SODEEnsemble(v, timespan, timestep, ics::AbstractVector{<: NamedTuple}; kwargs...)
SODEEnsemble(v, timespan, timestep, q₀::AbstractVector{<: StateVariable}; kwargs...)
SODEEnsemble(v, timespan, timestep, q₀::AbstractVector{<: AbstractArray}; kwargs...)
```
where `v` is a tuple of functions computing the vector fields for each substep,
`q` is an optional tuple of functions computing the solution for each substep,
`timespan` is the time interval `(t₀,t₁)` for the problem to be solved in,
`timestep` is the time step to be used in the simulation, and
`ics` is an `AbstractVector` of `NamedTuple`, each with an entry `q`
of type `StateVariable`.
The initial condition `q₀` can also be prescribed as an `AbstractVector`
of `StateVariable` or `AbstractArray{<:Number}`.
For the interface of the functions `v` and `q` see [`SODE`](@ref).

For possible keyword arguments see the documentation on [`EnsembleProblem`](@ref GeometricEquations.EnsembleProblem) subtypes.

"""
const SODEEnsemble  = EnsembleProblem{SODE}

function SODEEnsemble(v, q, timespan::Tuple, timestep::Real, ics...;
        invariants = NullInvariants(),
        parameters = NullParameters(),
        periodicity = NullPeriodicity(),
        v̄ = _sode_default_v̄)
    equ = SODE(v, q, v̄, invariants, parameter_types(parameters), periodicity)
    EnsembleProblem(equ, timespan, timestep, initialstate(equ, ics...), parameters)
end

function SODEEnsemble(v, args...; kwargs...)
    SODEEnsemble(v, nothing, args...; kwargs...)
end
