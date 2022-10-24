
const spsde_equations = raw"""
Defines a partitioned stochastic differential initial value problem
```math
\begin{aligned}
dq (t) &=   v(t, q(t), p(t)) \, dt + B(t, q(t), p(t)) \circ dW , & q(t_{0}) &= q_{0} , \\
dp (t) &= [ f_1(t, q(t), p(t)) + f_2(t, q(t), p(t)) ] \, dt + [ G_1(t, q(t), p(t)) + G_2(t, q(t), p(t)) ] \circ dW , & p(t_{0}) &= p_{0} ,
\end{aligned}
```
with the drift vector fields ``v`` and ``f_i``, diffusion matrices ``B`` and ``G_i``,
initial conditions ``q_{0}`` and ``p_{0}``, the dynamical variables ``(q,p)`` taking
values in ``\mathbb{R}^{d} \times \mathbb{R}^{d}``, and the m-dimensional Wiener process W
"""

const spsde_functions = raw"""
The functions `v`, `f1`, `f2, `B`, `G1` and `G2`, providing the drift vector fields and diffusion matrices,
all take five arguments, `(out, t, q, p, params)`.

```julia
function v(v, t, q, p, params)
    v[1] = ...
    v[2] = ...
    ...
end

function f1(f, t, q, p, params)
    f[1] = ...
    f[2] = ...
    ...
end

function f2(f, t, q, p, params)
    f[1] = ...
    f[2] = ...
    ...
end

function B(B, t, q, p, params)
    B[1,1] = ...
    ...
end

function G1(G, t, q, p, params)
    G[1,1] = ...
    ...
end

function G2(G, t, q, p, params)
    G[1,1] = ...
    ...
end
```

where `t` is the current time, `(q,p)` is the current solution vector, and `v`, `f`, `B` and `G`
are the variables which hold the result of evaluating the vector fields ``v``, ``f`` and the
matrices ``B_i``, ``G_i`` on `(t,q,p)`.
"""

const spsde_examples = raw"""

#### Example: Kubo Oscillator

```julia
function v(v, t, q, p, params)
    v[1] = + params.λ * p[1]
end

function f1(f, t, q, p, params)
    f[1] = - params.λ * q[1]
end

function f2(f, t, q, p, params)
    f[1] = 0
end

function B(B, t, q, p, params)
    B[1,1] = params.ν * p[1]
end

function G1(G, t, q, p, params)
    G[1,1] = - params.ν * q[1]
end

function G2(G, t, q, p, params)
    G[1,1] = 0
end

tspan = (0.0, 1.0); Δt = 0.01; q₀ = [0.5]; p₀ = [0.0];
prob = PSDEProblem(v, f1, f2, B, G1, G2, tspan, Δt, q₀, p₀; parameters = (λ=2., μ=1.))
```
"""

@doc """
`SPSDE`: Stratonovich Split Partitioned Stochastic Differential Equation

$(spsde_equations)

### Parameters

* `vType <: Function`: type of `v`
* `f1Type <: Function`: type of `f1`
* `f2Type <: Function`: type of `f2`
* `BType <: Function`: type of `B`
* `G1Type <: Function`: type of `G1`
* `G2Type <: Function`: type of `G2`
* `invType <: OptionalInvariants`: invariants type
* `parType <: OptionalParameters`: parameters type
* `perType <: OptionalPeriodicity`: periodicity type

### Fields

* `v` :  function computing the drift vector field for the position variable ``q``
* `f1`:  function computing the drift vector field for the momentum variable ``p``
* `f2`:  function computing the drift vector field for the momentum variable ``p``
* `B` :  function computing the d x m diffusion matrix for the position variable ``q``
* `G1`:  function computing the d x m diffusion matrix for the momentum variable ``p``
* `G2`:  function computing the d x m diffusion matrix for the momentum variable ``p``
* `invariants`: functions for the computation of invariants, either a `NamedTuple` containing the equation's invariants or `NullInvariants`
* `parameters`: type constraints for parameters, either a `NamedTuple` containing the equation's parameters or `NullParameters`
* `periodicity`: determines the periodicity of the state vector `q` for cutting periodic solutions, either a `AbstractArray` or `NullPeriodicity`

### Constructors

```julia
SPSDE(v, f1, f2, B, G1, G2, invariants, parameters, periodicity)
SPSDE(v, f1, f2, B, G1, G2; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
```

$(spsde_functions)

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
    applicable(equ.v,  zero(ics.q), tspan[begin], ics.q, ics.p, params) || return false
    applicable(equ.f1, zero(ics.p), tspan[begin], ics.q, ics.p, params) || return false
    applicable(equ.f2, zero(ics.p), tspan[begin], ics.q, ics.p, params) || return false
    # TODO add missing methods
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

_get_v(equ::SPSDE, params)  = hasparameters(equ) ? (v, t, q, p) -> equ.v(v, t, q, p, params) : equ.v
_get_f1(equ::SPSDE, params) = hasparameters(equ) ? (f, t, q, p) -> equ.f1(f, t, q, p, params) : equ.f1
_get_f2(equ::SPSDE, params) = hasparameters(equ) ? (f, t, q, p) -> equ.f2(f, t, q, p, params) : equ.f2
_get_B(equ::SPSDE, params)  = hasparameters(equ) ? (B, t, q, p) -> equ.B(B, t, q, p, params) : equ.B
_get_G1(equ::SPSDE, params) = hasparameters(equ) ? (G, t, q, p) -> equ.G1(G, t, q, p, params) : equ.G1
_get_G2(equ::SPSDE, params) = hasparameters(equ) ? (G, t, q, p) -> equ.G2(G, t, q, p, params) : equ.G2
_get_v̄(equ::SPSDE, params)  = _get_v(equ, params)
_get_f̄(equ::SPSDE, params)  = _get_f(equ, params)
_get_invariant(::SPSDE, inv, params) = (t, q, p) -> inv(t, q, p, params)

_functions(equ::SPSDE) = (v = equ.v, f1 = equ.f1, f2 = equ.f2, B = equ.B, G1 = equ.G1, G2 = equ.G2)
_functions(equ::SPSDE, params::OptionalParameters) = (
        v = _get_v(equ, params),
        f1 = _get_f1(equ, params),
        f2 = _get_f2(equ, params),
        B = _get_B(equ, params),
        G1 = _get_G1(equ, params),
        G2 = _get_G2(equ, params)
    )


@doc """
`SPSDEProblem`: Stratonovich Split Partitioned Stochastic Differential Equation Problem

$(spsde_equations)

### Constructors

```julia
SPSDEProblem(v, f1, f2, B, G1, G2, tspan, tstep, ics::NamedTuple; kwargs...)
SPSDEProblem(v, f1, f2, B, G1, G2, tspan, tstep, q₀::State; p₀::State; kwargs...)
```
where `v` and `f` are the functions computing the vector field and `Bᵢ` and `Gᵢ`
compute the diffusion matrices,
`tspan` is the time interval `(t₀,t₁)` for the problem to be solved in,
`tstep` is the time step to be used in the simulation, and
`ics` is a `NamedTuple` with entry `q`.
The initial condition `q₀` can also be prescribed directly, with
`State` an `AbstractArray{<:Number}`.

For possible keyword arguments see the documentation on [`GeometricProblem`](@ref GeometricEquations.GeometricProblem) subtypes.

### Function Definitions

$(spsde_functions)

"""
const SPSDEProblem = GeometricProblem{SPSDE}

function SPSDEProblem(v, f1, f2, B, G1, G2, tspan, tstep, ics::NamedTuple;
                      invariants = NullInvariants(), parameters = NullParameters(),
                      periodicity = NullPeriodicity())
    equ = SPSDE(v, f1, f2, B, G1, G2, invariants, parameter_types(parameters), periodicity)
    GeometricProblem(equ, tspan, tstep, ics, parameters)
end

function SPSDEProblem(v, f1, f2, B, G1, G2, tspan, tstep, q₀::State, p₀::State; kwargs...)
    ics = (q = q₀, p = p₀)
    SPSDEProblem(v, f1, f2, B, G1, G2, tspan, tstep, ics; kwargs...)
end

function GeometricBase.periodicity(prob::SPSDEProblem)
    (q = periodicity(equation(prob)), p = NullPeriodicity())
end
