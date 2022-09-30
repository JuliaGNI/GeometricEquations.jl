@doc raw"""
`ODE`: Ordinary Differential Equation

Defines an initial value problem
```math
\dot{q} (t) = v(t, q(t)) , \qquad q(t_{0}) = q_{0} ,
```
with vector field ``v``, initial condition ``q_{0}`` and the solution
``q`` taking values in ``\mathbb{R}^{d}``.

### Parameters

* `vType <: Callable`: type of `v`
* `invType <: OptionalNamedTuple`: invariants type
* `parType <: OptionalNamedTuple`: parameters type
* `perType <: OptionalArray{AT}`: periodicity type

### Fields

* `v`: function computing the vector field
* `invariants`: functions for the computation of invariants, either a `NamedTuple` containing the equation's invariants or `NullInvariants`
* `parameters`: type constraints for parameters, either a `NamedTuple` containing the equation's parameters or `NullParameters`
* `periodicity`: determines the periodicity of the state vector `q` for cutting periodic solutions, either a `AbstractArray` or `NullPeriodicity`

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
`q`, and `params` are additional parameters on which the vector field may depend.

### Constructors

```julia
ODE(v, invariants, parameters, periodicity)
ODE(v; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
```

"""
struct ODE{vType <: Callable,
           invType <: OptionalInvariants,
           parType <: OptionalParameters,
           perType <: OptionalPeriodicity} <: AbstractEquationODE{invType,parType,perType}

    v::vType

    invariants::invType
    parameters::parType
    periodicity::perType

    function ODE(v, invariants, parameters, periodicity)
        @assert !isempty(methods(v))
        # @assert hasmethod(v, (Real, AbstractArray, AbstractArray, OptionalParameters))

        new{typeof(v), typeof(invariants), typeof(parameters), typeof(periodicity)}(
                v, invariants, parameters, periodicity)
    end
end

ODE(v; invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity()) = ODE(v, invariants, parameters, periodicity)

GeometricBase.invariants(equation::ODE) = equation.invariants
GeometricBase.parameters(equation::ODE) = equation.parameters
GeometricBase.periodicity(equation::ODE) = equation.periodicity

hasvectorfield(::ODE) = true

function check_initial_conditions(::ODE, ics::NamedTuple)
    haskey(ics, :q) || return false
    return true
end

function check_methods(equ::ODE, tspan, ics, params)
    applicable(equ.v, zero(ics.q), tspan[begin], ics.q, params) || return false
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
