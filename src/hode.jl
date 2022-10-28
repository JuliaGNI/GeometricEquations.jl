
const hode_equations = raw"""
A canonical Hamiltonian system of equations is special case of a
partitioned ordinary differential equation,
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
taking values in ``T^{*} Q \simeq \mathbb{R}^{d} \times \mathbb{R}^{d}``.
"""

const hode_functions = raw"""
The functions `v`, `f` and `hamiltonian` must have the interface
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

    function hamiltonian(t, q, p, params)
        return ...
    end
```
where `t` is the current time, `q` and `p` are the current solution vectors,
`v` and `f` are the vectors which hold the result of evaluating the vector
fields on `t`, `q` and `p`, and params is a `NamedTuple` of
additional parameters.
"""

const hode_example = raw"""

#### Example: Harmonic Oscillator

As an example, let us consider the harmonic oscillator.
The dynamical equations are given by
```math
\begin{aligned}
\dot{q} (t) &= p (t) \\
\dot{p} (t) &= - k \, q(t) ,
\end{aligned}
```
which can also be written as 
```math
\begin{pmatrix}
\dot{q} (t) \\
\dot{p} (t) \\
\end{pmatrix} = \begin{pmatrix}
0 & 1 \\
-1 & 0 \\
\end{pmatrix}
\nabla H( q(t) , p(t) ) ,
\qquad
H(q,p) = \frac{p^2}{2} + k \, \frac{q^2}{2} ,
```
where $H$ is the Hamiltonian, i.e., the total energy of the system.

In order to create a `HODEProblem` for the harmonic oscillator, we need to write the following code:
```julia
function v(v, t, q, p, params)
    v[1] = p[1]
end

function f(f, t, q, p, params)
    f[1] = - params.k * q[1]
end

h(t, q, p, params) = p[1]^2 / 2 + params.k * q[1]^2 / 2

tspan = (0.0, 1.0)
tstep = 0.1
q₀ = [0.5]
p₀ = [0.0]

prob = HODEProblem(v, f, h, tspan, tstep, q₀, p₀; parameters = (k = 0.5,))
```
"""


@doc """
`HODE`: Hamiltonian Ordinary Differential Equation

$(hode_equations)

### Parameters

* `vType <: Callable`: type of `v`
* `fType <: Callable`: type of `f`
* `hamType <: Callable`: Hamiltonian type
* `invType <: OptionalInvariants`: invariants type
* `parType <: OptionalParameters`: parameters type
* `perType <: OptionalPeriodicity`: periodicity type

### Fields

* `v`: function computing the vector field ``v``
* `f`: function computing the vector field ``f``
* `hamiltonian`: function computing the Hamiltonian ``H`` (usually the total energy of the system)
* `invariants`: functions for the computation of invariants, either a `NamedTuple` containing the equation's invariants or `NullInvariants`
* `parameters`: type constraints for parameters, either a `NamedTuple` containing the equation's parameters or `NullParameters`
* `periodicity`: determines the periodicity of the state vector `q` for cutting periodic solutions, either a `AbstractArray` or `NullPeriodicity`

### Constructors

```julia
HODE(v, f, hamiltonian, invariants, parameters, periodicity)
HODE(v, f, hamiltonian; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
```

### Function Definitions

$(hode_functions)

"""
struct HODE{vType <: Callable, fType <: Callable,
            hamType <: Callable,
            invType <: OptionalInvariants,
            parType <: OptionalParameters,
            perType <: OptionalPeriodicity} <:
       AbstractEquationPODE{invType, parType, perType}
    v::vType
    f::fType

    hamiltonian::hamType
    invariants::invType
    parameters::parType
    periodicity::perType

    function HODE(v, f, hamiltonian, invariants, parameters, periodicity)
        new{typeof(v), typeof(f), typeof(hamiltonian),
            typeof(invariants), typeof(parameters), typeof(periodicity)}(v, f, hamiltonian,
                                                                         invariants,
                                                                         parameters,
                                                                         periodicity)
    end
end

function HODE(v, f, hamiltonian; invariants = NullInvariants(),
              parameters = NullParameters(), periodicity = NullPeriodicity())
    HODE(v, f, hamiltonian, invariants, parameters, periodicity)
end

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
    applicable(equ.v, zero(ics.q), tspan[begin], ics.q, ics.p, params) || return false
    applicable(equ.f, zero(ics.p), tspan[begin], ics.q, ics.p, params) || return false
    applicable(equ.hamiltonian, tspan[begin], ics.q, ics.p, params) || return false
    return true
end

function GeometricBase.datatype(equ::HODE, ics::NamedTuple)
    @assert check_initial_conditions(equ, ics)
    return eltype(ics.q)
end

function GeometricBase.arrtype(equ::HODE, ics::NamedTuple)
    @assert check_initial_conditions(equ, ics)
    typeof(ics.q)
end

_get_v(equ::HODE, params) = (v, t, q, p) -> equ.v(v, t, q, p, params)
_get_f(equ::HODE, params) = (f, t, q, p) -> equ.f(f, t, q, p, params)
_get_v̄(equ::HODE, params) = _get_v(equ, params)
_get_f̄(equ::HODE, params) = _get_f(equ, params)
_get_h(equ::HODE, params) = (t, q, p) -> equ.hamiltonian(t, q, p, params)
_get_invariant(::HODE, inv, params) = (t, q, p) -> inv(t, q, p, params)

_functions(equ::HODE) = (v = equ.v, f = equ.f, h = equ.hamiltonian)

function _functions(equ::HODE, params::OptionalParameters)
    (v = _get_v(equ, params), f = _get_f(equ, params), h = _get_h(equ, params))
end


@doc """
`HODEProblem`: Hamiltonian Ordinary Differential Equation Problem

$(hode_equations)

### Constructors

```julia
HODEProblem(v, f, hamiltonian, tspan, tstep, ics; kwargs...)
HODEProblem(v, f, hamiltonian, tspan, tstep, q₀::State, p₀::State; kwargs...)
```
where `v` and `f` are the function computing the vector fields, 
`hamiltonian` returns the value of the Hamiltonian (i.e. the total energy),
`tspan` is the time interval `(t₀,t₁)` for the problem to be solved in,
`tstep` is the time step to be used in the simulation, and
`ics` is a `NamedTuple` with entries `q` and `p`.
The initial conditions `q₀` and `p₀` can also be prescribed
directly, with `State` an `AbstractArray{<:Number}`.

### Keyword arguments:

* `invariants = NullInvariants()`
* `parameters = NullParameters()`
* `periodicity = NullPeriodicity()`

### Function Definitions

$(hode_functions)

"""
const HODEProblem = GeometricProblem{HODE}

function HODEProblem(v, f, hamiltonian, tspan, tstep, ics::NamedTuple;
                     invariants = NullInvariants(), parameters = NullParameters(),
                     periodicity = NullPeriodicity())
    equ = HODE(v, f, hamiltonian, invariants, parameter_types(parameters), periodicity)
    GeometricProblem(equ, tspan, tstep, ics, parameters)
end

function HODEProblem(v, f, hamiltonian, tspan, tstep, q₀::State, p₀::State; kwargs...)
    ics = (q = q₀, p = p₀)
    HODEProblem(v, f, hamiltonian, tspan, tstep, ics; kwargs...)
end

function GeometricBase.periodicity(prob::HODEProblem)
    (q = periodicity(equation(prob)), p = NullPeriodicity())
end
