@doc raw"""
`HODE`: Hamiltonian Ordinary Differential Equation

Defines a Hamiltonian ordinary differential initial value problem, that is
a canonical Hamiltonian system of equations,
```math
\begin{aligned}
\dot{q} (t) &= v(t, q(t), p(t)) , & q(t_{0}) &= q_{0} , \\
\dot{p} (t) &= f(t, q(t), p(t)) , & p(t_{0}) &= p_{0} ,
\end{aligned}
```
with vector fields ``v`` and ``f``, given by
```math
\begin{aligned}
v &=   \frac{\partial H}{\partial p} , &
f &= - \frac{\partial H}{\partial q} ,
\end{aligned}
```
initial conditions ``(q_{0}, p_{0})`` and the dynamical variables ``(q,p)``
taking values in ``\mathbb{R}^{d} \times \mathbb{R}^{d}``.

### Parameters

* `DT <: Number`: data type
* `TT <: Real`: time step type
* `AT <: AbstractArray{DT}`: array type
* `vType <: Function`: type of `v`
* `fType <: Function`: type of `f`
* `PType <: Function`: type of `P`
* `hamType <: Function`: Hamiltonian type
* `invType <: OptionalNamedTuple`: invariants type
* `parType <: OptionalNamedTuple`: parameters type
* `perType <: OptionalArray{AT}`: periodicity type

### Fields

* `d`: dimension of dynamical variables ``q`` and ``p`` as well as the vector fields ``v`` and ``f``
* `v`: function computing the vector field ``v``
* `f`: function computing the vector field ``f``
* `P`: function computing the Poisson matrix ``P``
* `t₀`: initial time
* `q₀`: initial condition for dynamical variable ``q``
* `p₀`: initial condition for dynamical variable ``p``
* `hamiltonian`: function computing the Hamiltonian ``H``
* `invariants`: either a `NamedTuple` containing the equation's invariants or `nothing`
* `parameters`: either a `NamedTuple` containing the equation's parameters or `nothing`
* `periodicity`: determines the periodicity of the state vector `q` for cutting periodic solutions

### Constructors

```julia
HODE(v, f, poisson, t₀, q₀, p₀, hamiltonian, invariants, parameters, periodicity)

HODE(v, f, h, t₀, q₀::StateVector, p₀::StateVector; kwargs...)
HODE(v, f, h, q₀::StateVector, p₀::StateVector; kwargs...)
HODE(v, f, h, t₀, q₀::State, p₀::State; kwargs...)
HODE(v, f, h, q₀::State, p₀::State; kwargs...)
```

### Keyword arguments:

* `poisson = symplectic_matrix`
* `invariants = nothing`
* `parameters = nothing`
* `periodicity = nothing`

"""
struct HODE{vType <: Callable, fType <: Callable, 
            poiType <: Callable,
            hamType <: Callable,
            invType <: OptionalInvariants,
            parType <: OptionalParameters,
            perType <: OptionalPeriodicity} <: AbstractEquationPODE{invType,parType,perType}

    v::vType
    f::fType

    poisson::poiType
    hamiltonian::hamType
    invariants::invType
    parameters::parType
    periodicity::perType

    function HODE(v, f, poisson, hamiltonian, invariants, parameters, periodicity)
        new{typeof(v), typeof(f), typeof(poisson), typeof(hamiltonian),
            typeof(invariants), typeof(parameters), typeof(periodicity)}(
                v, f, poisson, hamiltonian, invariants, parameters, periodicity)
    end
end

HODE(v, f, poisson, hamiltonian; invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity()) = HODE(v, f, poisson, hamiltonian, invariants, parameters, periodicity)

GeometricBase.invariants(equation::HODE) = equation.invariants
GeometricBase.parameters(equation::HODE) = equation.parameters
GeometricBase.periodicity(equation::HODE) = equation.periodicity

hasvectorfield(::HODE) = true
hashamiltonian(::HODE) = true

function check_initial_conditions(::HODE, ics::NamedTuple)
    haskey(ics, :q) || return false
    haskey(ics, :p) || return false
    eltype(ics.q) == eltype(ics.p) || return false
    typeof(ics.q) == typeof(ics.p) || return false
    axes(ics.q) == axes(ics.p) || return false
    return true
end

function check_methods(equ::HODE, tspan, ics, params)
    applicable(equ.v, tspan[begin], ics.q, ics.p, zero(ics.q), params) || return false
    applicable(equ.f, tspan[begin], ics.q, ics.p, zero(ics.p), params) || return false
    applicable(equ.hamiltonian, tspan[begin], ics.q, ics.p, params) || return false
    return true
end

function datatype(equ::HODE, ics::NamedTuple)
    @assert check_initial_conditions(equ, ics)
    return eltype(ics.q)
end

function arrtype(equ::HODE, ics::NamedTuple)
    @assert check_initial_conditions(equ, ics)
    typeof(ics.q)
end

_get_v(equ::HODE, params) = (t,q,p,v) -> equ.v(t, q, p, v, params)
_get_f(equ::HODE, params) = (t,q,p,f) -> equ.f(t, q, p, f, params)
_get_v̄(equ::HODE, params) = _get_v(equ, params)
_get_f̄(equ::HODE, params) = _get_f(equ, params)
_get_h(equ::HODE, params) = (t,q,p) -> equ.hamiltonian(t, q, p, params)
_get_poisson(equ::HODE, params) = (t,q,p,ω) -> equ.poisson(t, q, p, ω, params)

_functions(equ::HODE) = (v = equ.v, f = equ.f, poisson = equ.poisson, h = equ.hamiltonian)
_functions(equ::HODE, params::OptionalParameters) = (v = _get_v(equ, params), f = _get_f(equ, params), poisson = _get_poisson(equ, params), h = _get_h(equ, params))
