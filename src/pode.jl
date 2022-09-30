@doc raw"""
`PODE`: Partitioned Ordinary Differential Equation

Defines a partitioned initial value problem
```math
\begin{aligned}
\dot{q} (t) &= v(t, q(t), p(t)) , &
q(t_{0}) &= q_{0} , \\
\dot{p} (t) &= f(t, q(t), p(t)) , &
p(t_{0}) &= p_{0} ,
\end{aligned}
```
with vector fields ``v`` and ``f``, initial conditions ``(q_{0}, p_{0})`` and the solution
``(q,p)`` taking values in ``\mathbb{R}^{d} \times \mathbb{R}^{d}``.

### Parameters

* `vType <: Function`: type of `v`
* `fType <: Function`: type of `f`
* `invType <: OptionalNamedTuple`: invariants type
* `parType <: OptionalNamedTuple`: parameters type
* `perType <: OptionalArray{AT}`: periodicity type

### Fields

* `v`: function computing the vector field ``v``
* `f`: function computing the vector field ``f``
* `invariants`: functions for the computation of invariants, either a `NamedTuple` containing the equation's invariants or `NullInvariants`
* `parameters`: type constraints for parameters, either a `NamedTuple` containing the equation's parameters or `NullParameters`
* `periodicity`: determines the periodicity of the state vector `q` for cutting periodic solutions, either a `AbstractArray` or `NullPeriodicity`

The functions `v` and `f` must have the interface
```julia
    function v(v, t, q, p, params)
        v[1] = ...
        v[2] = ...
        ...
    end
```
and
```julia
    function f(f, t, q, p, params)
        f[1] = ...
        f[2] = ...
        ...
    end
```
where `t` is the current time, `q` and `p` are the current solution vectors,
`v` and `f` are the vectors which hold the result of evaluating the vector
fields ``v`` and ``f`` on `t`, `q` and `p`, and params are additional parameters.

### Constructors

```julia
PODE(v, f, invariants, parameters, periodicity)
PODE(v, f; invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity())
```


"""
struct PODE{vType <: Callable, fType <: Callable,
            invType <: OptionalInvariants,
            parType <: OptionalParameters,
            perType <: OptionalPeriodicity} <: AbstractEquationPODE{invType,parType,perType}

    v::vType
    f::fType

    invariants::invType
    parameters::parType
    periodicity::perType

    function PODE(v, f, invariants, parameters, periodicity)
        new{typeof(v), typeof(f), typeof(invariants), typeof(parameters), typeof(periodicity)}(
                v, f, invariants, parameters, periodicity)
    end
end

PODE(v, f; invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity()) = PODE(v, f, invariants, parameters, periodicity)

GeometricBase.invariants(equation::PODE) = equation.invariants
GeometricBase.parameters(equation::PODE) = equation.parameters
GeometricBase.periodicity(equation::PODE) = equation.periodicity

hasvectorfield(::PODE) = true

function check_initial_conditions(::PODE, ics::NamedTuple)
    haskey(ics, :q) || return false
    haskey(ics, :p) || return false
    eltype(ics.q) == eltype(ics.p) || return false
    typeof(ics.q) == typeof(ics.p) || return false
    axes(ics.q) == axes(ics.p) || return false
    return true
end

function check_methods(equ::PODE, tspan, ics, params)
    applicable(equ.v, zero(ics.q), tspan[begin], ics.q, ics.p, params) || return false
    applicable(equ.f, zero(ics.p), tspan[begin], ics.q, ics.p, params) || return false
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
_get_v̄(equ::PODE, params) = _get_v(equ, params)
_get_f̄(equ::PODE, params) = _get_f(equ, params)
_get_invariant(::PODE, inv, params) = (t,q,p) -> inv(t, q, p, params)

_functions(equ::PODE) = (v = equ.v, f = equ.f)
_functions(equ::PODE, params::OptionalParameters) = (v = _get_v(equ, params), f = _get_f(equ, params))
