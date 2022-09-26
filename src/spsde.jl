@doc raw"""
`SPSDE`: Stratonovich Split Partitioned Stochastic Differential Equation

Defines a partitioned stochastic differential initial value problem
```math
\begin{aligned}
dq (t) &=   v(t, q(t)) \, dt + B(t, q(t)) \circ dW , & q(t_{0}) &= q_{0} , \\
dp (t) &= [ f_1(t, q(t)) + f_2(t, q(t)) ] \, dt + [ G_1(t, q(t)) + G_2(t, q(t)) ] \circ dW , & p(t_{0}) &= p_{0} ,
\end{aligned}
```
with the drift vector fields ``v`` and ``f_i``, diffusion matrices ``B`` and ``G_i``,
initial conditions ``q_{0}`` and ``p_{0}``, the dynamical variables ``(q,p)`` taking
values in ``\mathbb{R}^{d} \times \mathbb{R}^{d}``, and the m-dimensional Wiener process W

### Parameters

* `DT <: Number`: data type
* `TT <: Real`: time step type
* `AT <: AbstractArray{DT}`: array type
* `vType <: Function`: type of `v`
* `f1Type <: Function`: type of `f1`
* `f2Type <: Function`: type of `f2`
* `BType <: Function`: type of `B`
* `G1Type <: Function`: type of `G1`
* `G2Type <: Function`: type of `G2`
* `pType <: Union{NamedTuple,Nothing}`: parameters type

### Fields

* `d`:  dimension of dynamical variable ``q`` and the vector fields ``vi``
* `m`:  dimension of the Wiener process
* `ni`: number of initial conditions
* `ns`: number of sample paths
* `v` :  function computing the drift vector field for the position variable ``q``
* `f1`:  function computing the drift vector field for the momentum variable ``p``
* `f2`:  function computing the drift vector field for the momentum variable ``p``
* `B` :  function computing the d x m diffusion matrix for the position variable ``q``
* `G1`:  function computing the d x m diffusion matrix for the momentum variable ``p``
* `G2`:  function computing the d x m diffusion matrix for the momentum variable ``p``
* `t₀`: initial time
* `q₀`: initial condition for dynamical variable ``q`` (may be a random variable itself)
* `p₀`: initial condition for dynamical variable ``p`` (may be a random variable itself)
* `parameters`: either a `NamedTuple` containing the equations parameters or `nothing`
* `periodicity`: determines the periodicity of the state vector `q` for cutting periodic solutions

The functions `v`, `f`, `B` and `G`, providing the drift vector fields and diffusion matrices, take four arguments,
`v(t, q, p, v)`, `f(t, q, p, f)`, `B(t, q, p,  B)` and `G(t, q, p, G)`, where `t` is the current time, `(q, p)` is the
current solution vector, and `v`, `f`, `B` and `G` are the variables which hold the result
of evaluating the vector fields ``v``, ``f`` and the matrices ``B``, ``G`` on `t` and `(q,p)`.

### Constructors

```julia
SPSDE(m, ns, v, f1, f2, B, G1, G2, t₀, q₀, p₀; parameters=NullParameters(), periodicity=zero(q₀[begin]))
SPSDE(m, ns, v, f1, f2, B, G1, G2, q₀::StateVector, p₀::StateVector; kwargs...) = SPSDE(m, ns, v, f1, f2, B, G1, G2, 0.0, q₀, p₀; kwargs...)
SPSDE(m, ns, v, f1, f2, B, G1, G2, t₀, q₀::State, p₀::State; kwargs...) = SPSDE(m, ns, v, f1, f2, B, G1, G2, t₀, [q₀], [p₀]; kwargs...)
SPSDE(m, ns, v, f1, f2, B, G1, G2, q₀::State, p₀::State; kwargs...) = SPSDE(m, ns, v, f1, f2, B, G1, G2, 0.0, q₀, p₀; kwargs...)
```

### Example

```julia
    function v(λ, t, q, v)
        v[1] = λ*q[1]
        v[2] = λ*q[2]
    end

    function B(μ, t, q, B)
        B[1] = μ*q[1]
        B[2] = μ*q[2]
    end

    t₀ = 0.
    q₀ = [1., 1.]
    λ  = 2.
    μ  = 1.

    v_sde = (t, q, v) -> v(λ, t, q, v)
    B_sde = (t, q, B) -> B(μ, t, q, B)

    sde = SDE(v_sde, B_sde, t₀, q₀)
```
"""
struct SPSDE{vType <: Callable,
             f1Type <: Callable,
             f2Type <: Callable,
             BType <: Callable,
             G1Type <: Callable,
             G2Type <: Callable,
             invType <: OptionalInvariants,
             parType <: OptionalParameters,
             perType <: OptionalPeriodicity} <: AbstractEquationPSDE{invType,parType,perType}
            
    v::vType
    f1::f1Type
    f2::f2Type
    B::BType
    G1::G1Type
    G2::G2Type

    invariants::invType
    parameters::parType
    periodicity::perType

    function SPSDE(v, f1, f2, B, G1, G2, invariants, parameters, periodicity)
        new{typeof(v), typeof(f1), typeof(f2), typeof(B), typeof(G1), typeof(G2), typeof(invariants), typeof(parameters), typeof(periodicity)}(
                v, f1, f2, B, G1, G2, invariants, parameters, periodicity)
    end
end

SPSDE(v, f1, f2, B, G1, G2; invariants=NullInvariants(), parameters=NullParameters(), periodicity=NullPeriodicity()) = SPSDE(v, f1, f2, B, G1, G2, invariants, parameters, periodicity)

GeometricBase.invariants(equation::SPSDE) = equation.invariants
GeometricBase.parameters(equation::SPSDE) = equation.parameters
GeometricBase.periodicity(equation::SPSDE) = equation.periodicity

hasvectorfield(::SPSDE) = true

function check_initial_conditions(::SPSDE, ics::NamedTuple)
    haskey(ics, :q) || return false
    haskey(ics, :p) || return false
    eltype(ics.q) == eltype(ics.p) || return false
    typeof(ics.q) == typeof(ics.p) || return false
    axes(ics.q) == axes(ics.p) || return false
    return true
end

function check_methods(equ::SPSDE, tspan, ics, params)
    applicable(equ.v, tspan[begin], ics.q, ics.p, zero(ics.q), params) || return false
    applicable(equ.f1, tspan[begin], ics.q, ics.p, zero(ics.p), params) || return false
    applicable(equ.f2, tspan[begin], ics.q, ics.p, zero(ics.p), params) || return false
    return true
end

function GeometricBase.datatype(equ::SPSDE, ics::NamedTuple)
    @assert check_initial_conditions(equ, ics)
    return eltype(ics.q)
end

function GeometricBase.arrtype(equ::SPSDE, ics::NamedTuple)
    @assert check_initial_conditions(equ, ics)
    typeof(ics.q)
end

_get_v(equ::SPSDE, params) = hasparameters(equ) ? (t,q,p,v) -> equ.v(t, q, p, v, params) : equ.v
_get_f1(equ::SPSDE, params) = hasparameters(equ) ? (t,q,p,f) -> equ.f1(t, q, p, f, params) : equ.f1
_get_f2(equ::SPSDE, params) = hasparameters(equ) ? (t,q,p,f) -> equ.f2(t, q, p, f, params) : equ.f2
_get_B(equ::SPSDE, params) = hasparameters(equ) ? (t,q,p,B) -> equ.B(t, q, p, B, params) : equ.B
_get_G1(equ::SPSDE, params) = hasparameters(equ) ? (t,q,p,G) -> equ.G1(t, q, p, G, params) : equ.G1
_get_G2(equ::SPSDE, params) = hasparameters(equ) ? (t,q,p,G) -> equ.G2(t, q, p, G, params) : equ.G2
_get_v̄(equ::SPSDE, params) = _get_v(equ, params)
_get_f̄(equ::SPSDE, params) = _get_f(equ, params)
_get_invariant(::SPSDE, inv, params) = (t,q,p) -> inv(t, q, p, params)

_functions(equ::SPSDE) = (v = equ.v, f1 = equ.f1, f2 = equ.f2, B = equ.B, G1 = equ.G1, G2 = equ.G2)
_functions(equ::SPSDE, params::OptionalParameters) = (
        v = _get_v(equ, params),
        f1 = _get_f1(equ, params),
        f2 = _get_f2(equ, params),
        B = _get_B(equ, params),
        G1 = _get_G1(equ, params),
        G2 = _get_G2(equ, params)
    )
