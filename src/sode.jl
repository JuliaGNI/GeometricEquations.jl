@doc raw"""
`SODE`: Split Ordinary Differential Equation

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

### Parameters

* `vType <: Union{Tuple,Nothing}`: type of `v`
* `qType <: Union{Tuple,Nothing}`: type of `q`
* `invType <: OptionalNamedTuple`: invariants type
* `parType <: OptionalNamedTuple`: parameters type
* `perType <: OptionalArray{AT}`: periodicity type

### Fields

* `v`: tuple of functions computing the vector field
* `q`: tuple of functions computing the solution
* `invariants`: either a `NamedTuple` containing the equation's invariants or `nothing`
* `parameters`: either a `NamedTuple` containing the equation's parameters or `nothing`
* `periodicity`: determines the periodicity of the state vector `q` for cutting periodic solutions

The functions `v_i` providing the vector field must have the interface
```julia
    function v_i(t, q, v, params)
        v[1] = ...
        v[2] = ...
        ...
    end
```
and the functions `q_i` providing the solutions must have the interface
```julia
    function q_i(t, q₀, q₁, h, params)
        q₁[1] = q₀[1] + ...
        q₁[2] = q₀[2] + ...
        ...
    end
```
where `t` is the current time, `q₀` is the current solution vector, `q₁` is the
new solution vector which holds the result of computing one substep with the
vector field ``v_i`` on `t` and `q₀`, and `h` is the (sub-)timestep to compute
the update for.

### Constructors

```julia
SODE(v, invariants, parameters, periodicity)
SODE(v; invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity())
SODE(v, q, invariants, parameters, periodicity)
SODE(v, q; invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity())
```

"""
struct SODE{vType <: Union{Tuple,Nothing}, qType <: Union{Tuple,Nothing},
            invType <: OptionalInvariants,
            parType <: OptionalParameters,
            perType <: OptionalPeriodicity} <: AbstractEquationODE{invType,parType,perType}

    v::vType
    q::qType

    invariants::invType
    parameters::parType
    periodicity::perType

    function SODE(v, q, invariants, parameters, periodicity)
        new{typeof(v), typeof(q), typeof(invariants), typeof(parameters), typeof(periodicity)}(
                v, q, invariants, parameters, periodicity)
    end
end

SODE(v; invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity()) = SODE(v, nothing, invariants, parameters, periodicity)
SODE(v, q; invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity()) = SODE(v, q, invariants, parameters, periodicity)

GeometricBase.invariants(equation::SODE) = equation.invariants
GeometricBase.parameters(equation::SODE) = equation.parameters
GeometricBase.periodicity(equation::SODE) = equation.periodicity

const SODEQT{QT,VT,invT,parT,perT} = SODE{VT,QT,invT,parT,perT} # type alias for dispatch on solution type parameter
const SODEVT{VT,QT,invT,parT,perT} = SODE{VT,QT,invT,parT,perT} # type alias for dispatch on vector field type parameter

hassolution(::SODEQT{<:Nothing}) = false
hassolution(::SODEQT{<:Tuple}) = true # && all(typeof(Q) <: Functiong for Q in equ.q)

hassolution(::SODEQT{<:Nothing}, i) = false
hassolution(equ::SODEQT{<:Tuple}, i) = i ≤ length(equ.q)# && typeof(equ.q[i]) <: Function

hasvectorfield(::SODEVT{<:Nothing}) = false
hasvectorfield(::SODEVT{<:Tuple}) = true # && all(typeof(V) <: Function for V in equ.v)

hasvectorfield(::SODEVT{<:Nothing}, i) = false
hasvectorfield(equ::SODEVT{<:Tuple}, i) = i ≤ length(equ.v)# && typeof(equ.v[i]) <: Function

function check_initial_conditions(::SODE, ics::NamedTuple)
    haskey(ics, :q) || return false
    return true
end

function check_methods(equ::SODE, tspan, ics, params)
    if hasvectorfield(equ)
        for v in equ.v
            applicable(v, tspan[begin], ics.q, zero(ics.q), params) || return false
        end
    end
    if hassolution(equ)
        for q in equ.q
            applicable(q, tspan[begin], ics.q, zero(ics.q), params) || return false
        end
    end
    return true
end

function datatype(equ::SODE, ics::NamedTuple)
    @assert check_initial_conditions(equ, ics)
    return eltype(ics.q)
end

function arrtype(equ::SODE, ics::NamedTuple)
    @assert check_initial_conditions(equ, ics)
    typeof(ics.q)
end

_get_v(equ::SODE, params) = ((t,q,v) -> V(t, q, v, params) for V in equ.v)
_get_q(equ::SODE, params) = ((t,q̄,q,h) -> Q(t, q̄, q, h, params) for Q in equ.q)

_functions(equ::SODE) = (v = equ.v,)
_solutions(equ::SODE) = (q = equ.q,)
_functions(equ::SODE, params::OptionalParameters) = (v = _get_v(equ, params),)
_solutions(equ::SODE, params::OptionalParameters) = (q = _get_q(equ, params),)
