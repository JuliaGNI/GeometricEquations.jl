
const iode_equations = raw"""
Defines an implicit initial value problem
```math
\begin{aligned}
\dot{q} (t) &= v(t) , &
q(t_{0}) &= q_{0} , \\
\dot{p} (t) &= f(t, q(t), v(t)) , &
p(t_{0}) &= p_{0} , \\
p(t) &= ϑ(t, q(t), v(t))
\end{aligned}
```
with force field ``f``, the momentum defined by ``p``, initial conditions ``(q_{0}, p_{0})``
and the solution ``(q,p)`` taking values in ``\mathbb{R}^{d} \times \mathbb{R}^{d}``.
This is a special case of a differential algebraic equation with dynamical
variables ``(q,p)`` and algebraic variable ``v``, that is determined such that the constraint
``p(t) = ϑ(t, q(t), v(t))`` is satisfied.

Most integrators perform a projection step in order to enforce this constraint. To this end,
the system is extended to
```math
\begin{aligned}
\dot{q} (t) &= v(t) + λ(t) , &
q(t_{0}) &= q_{0} , \\
\dot{p} (t) &= f(t, q(t), v(t)) + g(t, q(t), v(t), λ(t)) , &
p(t_{0}) &= p_{0} , \\
p(t) &= ϑ(t, q(t), v(t)) , &
λ(t_{0}) &= λ_{0}
\end{aligned}
```
where the vector field defining the projection step is usually given as
```math
\begin{aligned}
g(t, q(t), v(t), λ(t)) &= λ(t) \cdot \nabla ϑ(t, q(t), v(t)) .
\end{aligned}
```
"""

const iode_constructors = raw"""
The functions `ϑ`, `f` and `g` compute the momentum and the vector fields, respectively.
"""

const iode_functions = raw"""
The functions `ϑ` and `f` must have the interface
```julia
    function ϑ(p, t, q, v)
        p[1] = ...
        p[2] = ...
        ...
    end
```
and
```julia
    function f(f, t, q, v)
        f[1] = ...
        f[2] = ...
        ...
    end
```
where `t` is the current time, `q` is the current solution vector, `v` is the
current velocity and `f` and `p` are the vectors which hold the result of
evaluating the functions ``f`` and ``ϑ`` on `t`, `q` and `v`.
In addition, the functions `g`, `v̄` and `f̄` are specified by
```julia
    function g(g, t, q, v, λ)
        g[1] = ...
        g[2] = ...
        ...
    end

    function v̄(v, t, q)
        v[1] = ...
        v[2] = ...
        ...
    end

    function f̄(f, t, q, v)
        f[1] = ...
        f[2] = ...
        ...
    end
```
The function `g` is used in projection methods that enforce ``p = ϑ(q)``.
The functions `v̄` and `f̄` are used for initial guesses in nonlinear implicit solvers.
"""

@doc """
`IODE`: Implicit Ordinary Differential Equation

$(iode_equations)

### Parameters

* `ϑType <: Callable`: type of `ϑ`
* `fType <: Callable`: type of `f`
* `gType <: Callable`: type of `g`
* `v̄Type <: Callable`: type of `v̄`
* `f̄Type <: Callable`: type of `f̄`
* `invType <: OptionalInvariants`: invariants type
* `parType <: OptionalParameters`: parameters type
* `perType <: OptionalPeriodicity`: periodicity type

### Fields

* `ϑ`: function determining the momentum
* `f`: function computing the vector field
* `g`: function determining the projection, given by ``\\nabla \\vartheta (t,q,v) \\cdot \\lambda``
* `v̄`: function computing an initial guess for the velocity field ``v`` (optional)
* `f̄`: function computing an initial guess for the force field ``f`` (optional)
* `invariants`: functions for the computation of invariants, either a `NamedTuple` containing the equation's invariants or `NullInvariants`
* `parameters`: type constraints for parameters, either a `NamedTuple` containing the equation's parameters or `NullParameters`
* `periodicity`: determines the periodicity of the state vector `q` for cutting periodic solutions, either a `AbstractArray` or `NullPeriodicity`

### Constructors

```julia
IODE(ϑ, f, g, v̄, f̄, invariants, parameters, periodicity)
IODE(ϑ, f, g; v̄ = _iode_default_v̄, f̄ = f, invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
```

where 

```julia
_iode_default_v̄(v, t, q, params) = nothing
```

$(iode_constructors)

### Function Definitions

$(iode_functions)

"""
struct IODE{ϑType <: Callable, fType <: Callable, gType <: Callable,
            v̄Type <: Callable, f̄Type <: Callable,
            invType <: OptionalInvariants,
            parType <: OptionalParameters,
            perType <: OptionalPeriodicity} <:
       AbstractEquationPODE{invType, parType, perType}
    ϑ::ϑType
    f::fType
    g::gType
    v̄::v̄Type
    f̄::f̄Type

    invariants::invType
    parameters::parType
    periodicity::perType

    function IODE(ϑ, f, g, v̄, f̄, invariants, parameters, periodicity)
        new{typeof(ϑ), typeof(f), typeof(g), typeof(v̄), typeof(f̄),
            typeof(invariants), typeof(parameters), typeof(periodicity)}(ϑ, f, g, v̄, f̄,
                                                                         invariants,
                                                                         parameters,
                                                                         periodicity)
    end
end

_iode_default_v̄(v, t, q, params) = nothing

function IODE(ϑ, f, g; invariants = NullInvariants(), parameters = NullParameters(),
              periodicity = NullPeriodicity(), v̄ = _iode_default_v̄, f̄ = f)
    IODE(ϑ, f, g, v̄, f̄, invariants, parameters, periodicity)
end

GeometricBase.invariants(equation::IODE) = equation.invariants
GeometricBase.parameters(equation::IODE) = equation.parameters
GeometricBase.periodicity(equation::IODE) = equation.periodicity

hasvectorfield(::IODE) = true

function check_initial_conditions(::IODE, ics::NamedTuple)
    haskey(ics, :q) || return false
    haskey(ics, :p) || return false
    haskey(ics, :λ) || return false
    eltype(ics.q) == eltype(ics.p) == eltype(ics.λ) || return false
    typeof(ics.q) == typeof(ics.p) == typeof(ics.λ) || return false
    axes(ics.q) == axes(ics.p) == axes(ics.λ) || return false
    return true
end

function check_methods(equ::IODE, tspan, ics::NamedTuple, params)
    applicable(equ.ϑ, zero(ics.p), tspan[begin], ics.q, zero(ics.q), params) || return false
    applicable(equ.f, zero(ics.p), tspan[begin], ics.q, zero(ics.q), params) || return false
    applicable(equ.g, zero(ics.p), tspan[begin], ics.q, zero(ics.q), ics.λ, params) || return false
    applicable(equ.v̄, zero(ics.q), tspan[begin], ics.q, params) || return false
    applicable(equ.f̄, zero(ics.p), tspan[begin], ics.q, vectorfield(ics.q), params) || return false
    return true
end

function GeometricBase.datatype(equ::IODE, ics::NamedTuple)
    @assert check_initial_conditions(equ, ics)
    return eltype(ics.q)
end

function GeometricBase.arrtype(equ::IODE, ics::NamedTuple)
    @assert check_initial_conditions(equ, ics)
    typeof(ics.q)
end

_get_ϑ(equ::IODE, params) = (ϑ, t, q, v) -> equ.ϑ(ϑ, t, q, v, params)
_get_f(equ::IODE, params) = (f, t, q, v) -> equ.f(f, t, q, v, params)
_get_g(equ::IODE, params) = (g, t, q, v, λ) -> equ.g(g, t, q, v, λ, params)
_get_v̄(equ::IODE, params) = (v, t, q) -> equ.v̄(v, t, q, params)
_get_f̄(equ::IODE, params) = (f, t, q, v) -> equ.f̄(f, t, q, v, params)
_get_invariant(::IODE, inv, params) = (t, q, v) -> inv(t, q, v, params)

_functions(equ::IODE) = (ϑ = equ.ϑ, f = equ.f, g = equ.g, v̄ = equ.v̄, f̄ = equ.f̄)
function _functions(equ::IODE, params::OptionalParameters)
    (ϑ = _get_ϑ(equ, params),
     f = _get_f(equ, params),
     g = _get_g(equ, params),
     v̄ = _get_v̄(equ, params),
     f̄ = _get_f̄(equ, params))
end

@doc """
`IODEProblem`: Implicit Ordinary Differential Equation Problem
 
$(iode_equations)

### Constructors

```julia
IODEProblem(ϑ, f, g, tspan, tstep, ics; kwargs...)
IODEProblem(ϑ, f, g, tspan, tstep, q₀::State, p₀::State, λ₀::State = zero(q₀); kwargs...)
```

$(iode_constructors)

`tspan` is the time interval `(t₀,t₁)` for the problem to be solved in,
`tstep` is the time step to be used in the simulation, and
`ics` is a `NamedTuple` with entries `q` and `p`.
The initial conditions `q₀` and `p₀` can also be prescribed
directly, with `State` an `AbstractArray{<:Number}`.

In addition to the standard keyword arguments for [`GeometricProblem`](@ref GeometricEquations.GeometricProblem) subtypes,
an `IODEProblem` accepts functions `v̄` and `f̄` for the computation of initial guesses for the vector fields with default
values `v̄ = _iode_default_v̄` and `f̄ = f`.

The function `g` should really be optional as it is not required for all but only for most
integrators, but for the time being it is required.

### Function Definitions

$(iode_functions)

"""
const IODEProblem = GeometricProblem{IODE}

function IODEProblem(ϑ, f, g, tspan, tstep, ics::NamedTuple; invariants = NullInvariants(),
                     parameters = NullParameters(), periodicity = NullPeriodicity(),
                     v̄ = _iode_default_v̄, f̄ = f)
    equ = IODE(ϑ, f, g, v̄, f̄, invariants, parameter_types(parameters), periodicity)
    GeometricProblem(equ, tspan, tstep, ics, parameters)
end

function IODEProblem(ϑ, f, g, tspan, tstep, q₀::State, p₀::State, λ₀::State = zero(q₀);
                     kwargs...)
    ics = (q = q₀, p = p₀, λ = λ₀)
    IODEProblem(ϑ, f, g, tspan, tstep, ics; kwargs...)
end

function GeometricBase.periodicity(prob::IODEProblem)
    (q = periodicity(equation(prob)), p = NullPeriodicity(), λ = NullPeriodicity())
end
