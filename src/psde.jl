
const psde_equations = raw"""
A partitioned stochastic differential equations is an initial value problem of the form
```math
\begin{aligned}
dq (t) &= v(t, q(t), p(t)) \, dt + B(t, q(t), p(t)) \circ dW , & q(t_{0}) &= q_{0} , \\
dp (t) &= f(t, q(t), p(t)) \, dt + G(t, q(t), p(t)) \circ dW , & p(t_{0}) &= p_{0}
\end{aligned}
```
with the drift vector fields ``v`` and ``f``, diffusion matrices ``B`` and ``G``,
initial conditions ``q_{0}`` and ``p_{0}``, the dynamical variables ``(q,p)`` taking
values in ``\mathbb{R}^{d} \times \mathbb{R}^{d}``, and the m-dimensional Wiener process W
"""

const psde_functions = raw"""
The functions `v`, `f`, `B` and `G`, providing the drift vector fields and diffusion matrices, each take five arguments,
`v(v, t, q, p, params)`, `f(f, t, q, p, params)`, `B(B, t, q, p, params)` and `G(G, t, q, p, params)`,
where `t` is the current time, `(q, p)` is the current solution, and `v`, `f`, `B` and `G` are the variables which hold the result
of evaluating the vector fields ``v``, ``f`` and the matrices ``B``, ``G`` on `t` and `(q,p)`, and `params` are optional parameters.

The corresponding methods should have the following signatures:

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

function B(B, t, q, p, params)
    B[1,1] = ...
    ...
end

function G(G, t, q, p, params)
    G[1,1] = ...
    ...
end
```
"""


@doc """
`PSDE`: Stratonovich Partitioned Stochastic Differential Equation

$(psde_equations)

### Parameters

* `vType <: Callable`: type of `v`
* `fType <: Callable`: type of `f`
* `BType <: Callable`: type of `B`
* `GType <: Callable`: type of `G`
* `invType <: OptionalInvariants`: invariants type
* `parType <: OptionalParameters`: parameters type
* `perType <: OptionalPeriodicity`: periodicity type

### Fields

* `v`:  function computing the drift vector field for the position variable ``q``
* `f`:  function computing the drift vector field for the momentum variable ``p``
* `B`:  function computing the d x m diffusion matrix for the position variable ``q``
* `G`:  function computing the d x m diffusion matrix for the momentum variable ``p``
* `invariants`: functions for the computation of invariants, either a `NamedTuple` containing the equation's invariants or `NullInvariants`
* `parameters`: type constraints for parameters, either a `NamedTuple` containing the equation's parameters or `NullParameters`
* `periodicity`: determines the periodicity of the state vector `q` for cutting periodic solutions, either a `AbstractArray` or `NullPeriodicity`

### Constructors

```julia
PSDE(v, f, B, G, invariants, parameters, periodicity)
PSDE(v, f, B, G; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
```

$(psde_functions)

"""
struct PSDE{vType <: Callable,
            fType <: Callable,
            BType <: Callable,
            GType <: Callable,
            invType <: OptionalInvariants,
            parType <: OptionalParameters,
            perType <: OptionalPeriodicity} <: AbstractEquationPSDE{invType,parType,perType}

    v::vType
    f::fType
    B::BType
    G::GType

    invariants::invType
    parameters::parType
    periodicity::perType

    function PSDE(v, f, B, G, invariants, parameters, periodicity)
        new{typeof(v), typeof(f), typeof(B), typeof(G), typeof(invariants), typeof(parameters), typeof(periodicity)}(
                v, f, B, G, invariants, parameters, periodicity)
    end
end

PSDE(v, f, B, G; invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity()) = PSDE(v, f, B, G, invariants, parameters, periodicity)

GeometricBase.invariants(equation::PSDE) = equation.invariants
GeometricBase.parameters(equation::PSDE) = equation.parameters
GeometricBase.periodicity(equation::PSDE) = equation.periodicity

hasvectorfield(::PSDE) = true

function check_initial_conditions(::PSDE, ics::NamedTuple)
    haskey(ics, :q) || return false
    haskey(ics, :p) || return false
    eltype(ics.q) == eltype(ics.p) || return false
    typeof(ics.q) == typeof(ics.p) || return false
    axes(ics.q) == axes(ics.p) || return false
    return true
end

function check_methods(equ::PSDE, tspan, ics, params)
    applicable(equ.v, zero(ics.q), tspan[begin], ics.q, ics.p, params) || return false
    applicable(equ.f, zero(ics.p), tspan[begin], ics.q, ics.p, params) || return false
    return true
end

function GeometricBase.datatype(equ::PSDE, ics::NamedTuple)
    @assert check_initial_conditions(equ, ics)
    return eltype(ics.q)
end

function GeometricBase.arrtype(equ::PSDE, ics::NamedTuple)
    @assert check_initial_conditions(equ, ics)
    typeof(ics.q)
end

_get_v(equ::PSDE, params) = hasparameters(equ) ? (v, t, q, p) -> equ.v(v, t, q, p, params) : equ.v
_get_f(equ::PSDE, params) = hasparameters(equ) ? (f, t, q, p) -> equ.f(f, t, q, p, params) : equ.f
_get_B(equ::PSDE, params) = hasparameters(equ) ? (B, t, q, p) -> equ.B(B, t, q, p, params) : equ.B
_get_G(equ::PSDE, params) = hasparameters(equ) ? (G, t, q, p) -> equ.G(G, t, q, p, params) : equ.G
_get_v̄(equ::PSDE, params) = _get_v(equ, params)
_get_f̄(equ::PSDE, params) = _get_f(equ, params)
_get_invariant(::PSDE, inv, params) = (t, q, p) -> inv(t, q, p, params)

_functions(equ::PSDE) = (v = equ.v, f = equ.f, B = equ.B, G = equ.G)
_functions(equ::PSDE, params::OptionalParameters) = (v = _get_v(equ, params), f = _get_f(equ, params), B = _get_B(equ, params), G = _get_G(equ, params))


@doc """
`PSDEProblem`: Stratonovich Partitioned Stochastic Differential Equation Problem

$(psde_equations)

### Constructors

```julia
PSDEProblem(v, f, B, G, tspan, tstep, ics::NamedTuple; kwargs...)
PSDEProblem(v, f, B, G, tspan, tstep, q₀::State; p₀::State; kwargs...)
```
where `v` and `f` are the functions computing the vector field and `B` and `G`
compute the diffusion matrices,
`tspan` is the time interval `(t₀,t₁)` for the problem to be solved in,
`tstep` is the time step to be used in the simulation, and
`ics` is a `NamedTuple` with entry `q`.
The initial condition `q₀` can also be prescribed directly, with
`State` an `AbstractArray{<:Number}`.

For possible keyword arguments see the documentation on [`GeometricProblem`](@ref GeometricEquations.GeometricProblem) subtypes.

### Function Definitions

$(psde_functions)

"""
const PSDEProblem = GeometricProblem{PSDE}

function PSDEProblem(v, f, B, G, tspan, tstep, ics::NamedTuple;
                     invariants = NullInvariants(), parameters = NullParameters(),
                     periodicity = NullPeriodicity())
    equ = PSDE(v, f, B, G, invariants, parameter_types(parameters), periodicity)
    GeometricProblem(equ, tspan, tstep, ics, parameters)
end

function PSDEProblem(v, f, B, G, tspan, tstep, q₀::State, p₀::State; kwargs...)
    ics = (q = q₀, p = p₀)
    PSDEProblem(v, f, B, G, tspan, tstep, ics; kwargs...)
end

function GeometricBase.periodicity(prob::PSDEProblem)
    (q = periodicity(equation(prob)), p = NullPeriodicity())
end
