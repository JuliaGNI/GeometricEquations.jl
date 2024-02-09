
const sde_equations = raw"""
Defines a stochastic differential initial value problem
```math
\begin{aligned}
dq (t) &= v(t, q(t)) \, dt + B(t, q(t)) \circ dW , & q(t_{0}) &= q_{0} ,
\end{aligned}
```
with drift vector field ``v``, diffusion matrix ``B``,
initial conditions ``q_{0}``, the dynamical variable ``q``
taking values in ``\mathbb{R}^{d}``, and the m-dimensional Wiener process W
"""

const sde_functions = raw"""
The functions `v` and `B`, providing the drift vector field and diffusion matrix.
The function `v` must have the interface
```julia
function v(v, t, q, params)
    v[1] = ...
    v[2] = ...
    ...
end
```
where `t` is the current time, `q` is the current solution vector, `v` is the
vector which holds the result of evaluating the vector field ``v`` on `t` and
`q`, and params are additional parameters.
The function `B` should have a method with interface
```julia
function B(B, t, q, params)
    B[1,1] = ...
    ...
end
```
"""


@doc """
`SDE`: Stratonovich Stochastic Differential Equation

$(sde_equations)

### Parameters

* `vType <: Callable`: type of `v`
* `BType <: Callable`: type of `B`
* `invType <: OptionalInvariants`: invariants type
* `parType <: OptionalParameters`: parameters type
* `perType <: OptionalPeriodicity`: periodicity type

### Fields

* `v`:  function computing the deterministic vector field
* `B`:  function computing the d x m diffusion matrix
* `invariants`: functions for the computation of invariants, either a `NamedTuple` containing the equation's invariants or `NullInvariants`
* `parameters`: type constraints for parameters, either a `NamedTuple` containing the equation's parameters or `NullParameters`
* `periodicity`: determines the periodicity of the state vector `q` for cutting periodic solutions, either a `AbstractArray` or `NullPeriodicity`

### Constructors

```julia
SDE(v, B, invariants, parameters, periodicity)
SDE(v, B; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
```

$(sde_functions)

"""
struct SDE{vType <: Callable,
           BType <: Callable,
           nType <: AbstractStochasticProcess,
           invType <: OptionalInvariants,
           parType <: OptionalParameters,
           perType <: OptionalPeriodicity} <: AbstractEquationSDE{invType,parType,perType}

    v::vType
    B::BType

    noise::nType

    invariants::invType
    parameters::parType
    periodicity::perType

    function SDE(v, B, noise, invariants, parameters, periodicity)
        new{typeof(v), typeof(B), typeof(noise), typeof(invariants), typeof(parameters), typeof(periodicity)}(
                v, B, noise, invariants, parameters, periodicity)
    end
end

SDE(v, B, noise; invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity()) = SDE(v, B, noise, invariants, parameters, periodicity)

GeometricBase.invariants(equation::SDE) = equation.invariants
GeometricBase.parameters(equation::SDE) = equation.parameters
GeometricBase.periodicity(equation::SDE) = equation.periodicity

hasvectorfield(::SDE) = true

function Base.show(io::IO, equation::SDE)
    print(io, "Stratonovich Stochastic Differential Equation (SDE)", "\n")
    print(io, "\n")
    print(io, " with vector field")
    print(io, "\n")
    print(io, "   v = ", equation.v, "\n")
    print(io, "\n")
    print(io, " and diffusion matrix")
    print(io, "\n")
    print(io, "   B = ", equation.B, "\n")
    print(io, "\n")
    print(io, " Invariants: \n")
    print(io, "   ", invariants(equation))
end

function check_initial_conditions(::SDE, ics::NamedTuple)
    haskey(ics, :q) || return false
    return true
end

function check_methods(equ::SDE, tspan, ics, params)
    applicable(equ.v, zero(ics.q), tspan[begin], ics.q, params) || return false
    # TODO add missing methods
    return true
end

function GeometricBase.datatype(equ::SDE, ics::NamedTuple)
    @assert check_initial_conditions(equ, ics)
    return eltype(ics.q)
end

function GeometricBase.arrtype(equ::SDE, ics::NamedTuple)
    @assert check_initial_conditions(equ, ics)
    typeof(ics.q)
end

_get_v(equ::SDE, params) = (v, t, q) -> equ.v(v, t, q, params)
_get_B(equ::SDE, params) = (B, t, q) -> equ.B(B, t, q, params)
_get_v̄(equ::SDE, params) = _get_v(equ, params)
_get_invariant(::SDE, inv, params) = (t, q) -> inv(t, q, params)

_functions(equ::SDE) = (v = equ.v, B = equ.B)
_functions(equ::SDE, params::OptionalParameters) = (v = _get_v(equ, params), B = _get_B(equ, params))


@doc """
`SDEProblem`: Stratonovich Stochastic Differential Equation Problem

$(sde_equations)

### Constructors

```julia
SDEProblem(v, B, tspan, tstep, ics::NamedTuple; kwargs...)
SDEProblem(v, B, tspan, tstep, q₀::StateVariable; kwargs...)
```
where `v` is the function computing the vector field and `B` computes the diffusion matrix
`tspan` is the time interval `(t₀,t₁)` for the problem to be solved in,
`tstep` is the time step to be used in the simulation, and
`ics` is a `NamedTuple` with entry `q`.
The initial condition `q₀` can also be prescribed directly, with
`StateVariable` an `AbstractArray{<:Number}`.

For possible keyword arguments see the documentation on [`EquationProblem`](@ref GeometricEquations.EquationProblem) subtypes.

### Function Definitions

$(sde_functions)

"""
const SDEProblem = EquationProblem{SDE}

function SDEProblem(v, B, noise, tspan, tstep, ics::NamedTuple; invariants = NullInvariants(),
                    parameters = NullParameters(), periodicity = NullPeriodicity())
    equ = SDE(v, B, noise, invariants, parameter_types(parameters), periodicity)
    EquationProblem(equ, tspan, tstep, ics, parameters)
end

function SDEProblem(v, B, noise, tspan, tstep, q₀::StateVariable; kwargs...)
    ics = (q = q₀,)
    SDEProblem(v, B, noise, tspan, tstep, ics; kwargs...)
end

GeometricBase.periodicity(prob::SDEProblem) = (q = periodicity(equation(prob)),)


const SDEEnsemble = EnsembleProblem{SDE}
