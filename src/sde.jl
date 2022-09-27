@doc raw"""
`SDE`: Stratonovich Stochastic Differential Equation

Defines a stochastic differential initial value problem
```math
\begin{aligned}
dq (t) &= v(t, q(t)) \, dt + B(t, q(t)) \circ dW , & q(t_{0}) &= q_{0} ,
\end{aligned}
```
with drift vector field ``v``, diffusion matrix ``B``,
initial conditions ``q_{0}``, the dynamical variable ``q``
taking values in ``\mathbb{R}^{d}``, and the m-dimensional Wiener process W

### Parameters

* `vType <: Callable`: type of `v`
* `BType <: Callable`: type of `B`
* `invType <: OptionalNamedTuple`: invariants type
* `parType <: OptionalNamedTuple`: parameters type
* `perType <: OptionalArray{AT}`: periodicity type

### Fields

* `v`:  function computing the deterministic vector field
* `B`:  function computing the d x m diffusion matrix
* `invariants`: functions for the computation of invariants, either a `NamedTuple` containing the equation's invariants or `NullInvariants`
* `parameters`: type constraints for parameters, either a `NamedTuple` containing the equation's parameters or `NullParameters`
* `periodicity`: determines the periodicity of the state vector `q` for cutting periodic solutions, either a `AbstractArray` or `NullPeriodicity`

The functions `v` and `B`, providing the drift vector field and diffusion matrix.
The function `v` must have the interface
```julia
    function v(t, q, v, params)
        v[1] = ...
        v[2] = ...
        ...
    end
```
where `t` is the current time, `q` is the current solution vector, `v` is the
vector which holds the result of evaluating the vector field ``v`` on `t` and
`q`, and params are additional parameters.
The function `B` should have two methods with interfaces
    ```julia
    function B(t, q, B::AbstractVector, params, col)
        ...
    end

    function B(t, q, B::AbstractMatrix, params)
        ...
    end
```
where the first method returns a single column of B and the second method returns
the whole matrix.

### Constructors

```julia
SDE(v, B, invariants, parameters, periodicity)
SDE(v, B; invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity())
```
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
    applicable(equ.v, tspan[begin], ics.q, zero(ics.q), params) || return false
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

_get_v(equ::SDE, params) = (t,q,v) -> equ.v(t, q, v, params)
_get_B(equ::SDE, params) = (t,q,v) -> equ.B(t, q, B, params)
_get_v̄(equ::SDE, params) = _get_v(equ, params)
_get_invariant(::SDE, inv, params) = (t,q) -> inv(t, q, params)

_functions(equ::SDE) = (v = equ.v, B = equ.B)
_functions(equ::SDE, params::OptionalParameters) = (v = _get_v(equ, params), B = _get_B(equ, params))
