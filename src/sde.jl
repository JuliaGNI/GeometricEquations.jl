
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

const sde_examples = raw"""

#### Example: Kubo Oscillator

```julia
function v(t, q, v, params)
    v[1] = + params.λ * q[2]
    v[2] = - params.λ * q[1]
end

function B(t, q, B, params)
    for j in axes(B, 2)
        B[1,j] = + params.ν * q[2]
        B[2,j] = - params.ν * q[1]
    end
end

tspan = (0.0, 1.0); Δt = 0.01; q₀ = [0.5, 0.0];
prob = SDEProblem(v, B, tspan, Δt, q₀; parameters = (λ=2., μ=1.))
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
           invType <: OptionalInvariants,
           parType <: OptionalParameters,
           perType <: OptionalPeriodicity} <: AbstractEquationSDE{invType,parType,perType}

    v::vType
    B::BType

    invariants::invType
    parameters::parType
    periodicity::perType

    function SDE(v, B, invariants, parameters, periodicity)
        new{typeof(v), typeof(B), typeof(invariants), typeof(parameters), typeof(periodicity)}(
                v, B, invariants, parameters, periodicity)
    end
end

SDE(v, B; invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity()) = SDE(v, B, invariants, parameters, periodicity)

GeometricBase.invariants(equation::SDE) = equation.invariants
GeometricBase.parameters(equation::SDE) = equation.parameters
GeometricBase.periodicity(equation::SDE) = equation.periodicity

hasvectorfield(::SDE) = true

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
SDEProblem(v, B, tspan, tstep, q₀::State; kwargs...)
```
where `v` is the function computing the vector field and `B` computes the diffusion matrix
`tspan` is the time interval `(t₀,t₁)` for the problem to be solved in,
`tstep` is the time step to be used in the simulation, and
`ics` is a `NamedTuple` with entry `q`.
The initial condition `q₀` can also be prescribed directly, with
`State` an `AbstractArray{<:Number}`.

For possible keyword arguments see the documentation on [`GeometricProblem`](@ref GeometricEquations.GeometricProblem) subtypes.

### Function Definitions

$(sde_functions)

"""
const SDEProblem = GeometricProblem{SDE}

function SDEProblem(v, B, tspan, tstep, ics::NamedTuple; invariants = NullInvariants(),
                    parameters = NullParameters(), periodicity = NullPeriodicity())
    equ = SDE(v, B, invariants, parameter_types(parameters), periodicity)
    GeometricProblem(equ, tspan, tstep, ics, parameters)
end

function SDEProblem(v, B, tspan, tstep, q₀::State; kwargs...)
    ics = (q = q₀,)
    SDEProblem(v, B, tspan, tstep, ics; kwargs...)
end

GeometricBase.periodicity(prob::SDEProblem) = (q = periodicity(equation(prob)),)
