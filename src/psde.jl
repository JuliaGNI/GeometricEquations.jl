@doc raw"""
`PSDE`: Stratonovich Partitioned Stochastic Differential Equation

Defines a partitioned stochastic differential initial value problem
```math
\begin{aligned}
dq (t) &= v(t, q(t)) \, dt + B(t, q(t)) \circ dW , & q(t_{0}) &= q_{0} , \\
dp (t) &= f(t, q(t)) \, dt + G(t, q(t)) \circ dW , & p(t_{0}) &= p_{0}
\end{aligned}
```
with the drift vector fields ``v`` and ``f``, diffusion matrices ``B`` and ``G``,
initial conditions ``q_{0}`` and ``p_{0}``, the dynamical variables ``(q,p)`` taking
values in ``\mathbb{R}^{d} \times \mathbb{R}^{d}``, and the m-dimensional Wiener process W

### Parameters

* `DT <: Number`: data type
* `TT <: Real`: time step type
* `AT <: AbstractArray{DT}`: array type
* `vType <: Function`: type of `v`
* `fType <: Function`: type of `f`
* `BType <: Function`: type of `B`
* `GType <: Function`: type of `G`
* `pType <: Union{NamedTuple,Nothing}`: parameters type

### Fields

* `d`:  dimension of dynamical variable ``q`` and the vector field ``v``
* `m`:  dimension of the Wiener process
* `ns`: number of sample paths
* `v`:  function computing the drift vector field for the position variable ``q``
* `f`:  function computing the drift vector field for the momentum variable ``p``
* `B`:  function computing the d x m diffusion matrix for the position variable ``q``
* `G`:  function computing the d x m diffusion matrix for the momentum variable ``p``
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
PSDE(m, ns, v, f, B, G, t₀, q₀, p₀; parameters=NullParameters(), periodicity=zero(q₀[begin]))
PSDE(m, ns, v, f, B, G, q₀::StateVector, p₀::StateVector; kwargs...) = PSDE(m, ns, v, f, B, G, 0.0, q₀, p₀; kwargs...)
PSDE(m, ns, v, f, B, G, t₀, q₀::State, p₀::State; kwargs...) = PSDE(m, ns, v, f, B, G, t₀, [q₀], [p₀]; kwargs...)
PSDE(m, ns, v, f, B, G, q₀::State, p₀::State; kwargs...) = PSDE(m, ns, v, f, B, G, 0.0, q₀, p₀; kwargs...)
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
    applicable(equ.v, tspan[begin], ics.q, ics.p, zero(ics.q), params) || return false
    applicable(equ.f, tspan[begin], ics.q, ics.p, zero(ics.p), params) || return false
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

_get_v(equ::PSDE, params) = hasparameters(equ) ? (t,q,p,v) -> equ.v(t, q, p, v, params) : equ.v
_get_f(equ::PSDE, params) = hasparameters(equ) ? (t,q,p,f) -> equ.f(t, q, p, f, params) : equ.f
_get_B(equ::PSDE, params) = hasparameters(equ) ? (t,q,p,B) -> equ.B(t, q, p, B, params) : equ.B
_get_G(equ::PSDE, params) = hasparameters(equ) ? (t,q,p,G) -> equ.G(t, q, p, G, params) : equ.G
_get_v̄(equ::PSDE, params) = _get_v(equ, params)
_get_f̄(equ::PSDE, params) = _get_f(equ, params)
_get_invariant(::PSDE, inv, params) = (t,q,p) -> inv(t, q, p, params)

_functions(equ::PSDE) = (v = equ.v, f = equ.f, B = equ.B, G = equ.G)
_functions(equ::PSDE, params::OptionalParameters) = (v = _get_v(equ, params), f = _get_f(equ, params), B = _get_B(equ, params), G = _get_G(equ, params))
