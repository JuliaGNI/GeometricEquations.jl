
const hode_equations = raw"""
A canonical Hamiltonian system of equations is special case of a
partitioned ordinary differential equation,
```math
\begin{aligned}
\dot{q} (t) &= v(t, q(t), p(t)) , \\
\dot{p} (t) &= f(t, q(t), p(t)) ,
\end{aligned}
```
with vector fields ``v`` and ``f``, given by
```math
\begin{aligned}
v &=   \frac{\partial H}{\partial p} , &
f &= - \frac{\partial H}{\partial q} .
\end{aligned}
```
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

function Base.show(io::IO, equation::HODE)
    print(io, "Hamiltonian Ordinary Differential Equation (HODE)", "\n")
    print(io, "\n")
    print(io, " with vector fields", "\n")
    print(io, "   v = ", equation.v, "\n")
    print(io, "   f = ", equation.f, "\n")
    print(io, "\n")
    print(io, " Hamiltonian: H = ", equation.hamiltonian, "\n")
    print(io, "\n")
    print(io, " Invariants: \n")
    print(io, "   ", invariants(equation))
end

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

The dynamical variables ``(q,p)`` with initial conditions ``(q(t_{0}) = q_{0}, p(t_{0}) = p_{0})``
take values in ``T^{*} Q \\simeq \\mathbb{R}^{d} \\times \\mathbb{R}^{d}``.

### Constructors

```julia
HODEProblem(v, f, hamiltonian, tspan, tstep, ics; kwargs...)
HODEProblem(v, f, hamiltonian, tspan, tstep, q₀::StateVariable, p₀::StateVariable; kwargs...)
```
where `v` and `f` are the function computing the vector fields, 
`hamiltonian` returns the value of the Hamiltonian (i.e. the total energy),
`tspan` is the time interval `(t₀,t₁)` for the problem to be solved in,
`tstep` is the time step to be used in the simulation, and
`ics` is a `NamedTuple` with entries `q` and `p`.
The initial conditions `q₀` and `p₀` can also be prescribed
directly, with `StateVariable` an `AbstractArray{<:Number}`.
For the interfaces of the functions `v`, `f`, `poisson` and `hamiltonian` see [`HODE`](@ref).

For possible keyword arguments see the documentation on [`EquationProblem`](@ref GeometricEquations.EquationProblem) subtypes.

### Function Definitions

$(hode_functions)

"""
const HODEProblem = EquationProblem{HODE}

function HODEProblem(v, f, hamiltonian, tspan, tstep, ics::NamedTuple;
        invariants = NullInvariants(),
        parameters = NullParameters(),
        periodicity = NullPeriodicity())
    equ = HODE(v, f, hamiltonian, invariants, parameter_types(parameters), periodicity)
    EquationProblem(equ, tspan, tstep, ics, parameters)
end

function HODEProblem(v, f, hamiltonian, tspan, tstep, q₀::StateVariable, p₀::StateVariable; kwargs...)
    ics = (q = q₀, p = p₀)
    HODEProblem(v, f, hamiltonian, tspan, tstep, ics; kwargs...)
end

function GeometricBase.periodicity(prob::HODEProblem)
    (q = periodicity(equation(prob)), p = NullPeriodicity())
end


const HODEEnsemble  = EnsembleProblem{HODE}

function HODEEnsemble(v, f, hamiltonian, tspan, tstep, ics::AbstractVector{<:NamedTuple};
        invariants = NullInvariants(),
        parameters = NullParameters(),
        periodicity = NullPeriodicity())
    equ = HODE(v, f, hamiltonian, invariants, parameter_types(parameters), periodicity)
    EnsembleProblem(equ, tspan, tstep, ics, parameters)
end
