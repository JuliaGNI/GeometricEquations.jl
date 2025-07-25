
const hdae_equations = raw"""
A Hamiltonian differential algebraic is an initial value problem, that is
a canonical Hamiltonian system of equations subject to Dirac constraints,
```math
\begin{aligned}
\dot{q} (t) &= v(t, q(t), p(t)) + u(t, q(t), p(t), \lambda(t)) + \bar{u} (t, q(t), p(t), \mu(t)) , \\
\dot{p} (t) &= f(t, q(t), p(t)) + g(t, q(t), p(t), \lambda(t)) + \bar{g} (t, q(t), p(t), \mu(t)) , \\
0 &= \phi (t, q(t), p(t)) , \\
0 &= \psi (t, q(t), p(t), \dot{q}(t), \dot{p}(t)) ,
\end{aligned}
```
with vector fields ``v``, ``u``, ``\bar{u}`` and ``f``, ``g``, ``\bar{g}``,
primary constraint ``\phi(q,p)=0`` and secondary constraint ``\psi(q,p,\dot{q},\dot{p})=0``.
"""

const hdae_constructors = raw"""
The functions `v` and `f` compute the vector field, `u` and `g` compute the projections,
`ϕ` provides the algebraic constraint and `h` the Hamiltonian.
The functions `ψ`, `ū` and `ḡ` are optional and provide the secondary constraint, that
is the time derivative of the algebraic constraint, and the corresponding projection.
"""

const hdae_functions = raw"""
The functions `v`, `f`, `u`, `g`, `ϕ` and `h` must have the interface

```julia
function v(v, t, q, p, params)
    v[1] = ...
    v[2] = ...
    ...
end

function f(g, t, q, p, params)
    f[1] = ...
    f[2] = ...
    ...
end

function u(u, t, q, p, λ, params)
    u[1] = ...
    u[2] = ...
    ...
end

function g(g, t, q, p, λ, params)
    g[1] = ...
    g[2] = ...
    ...
end

function ϕ(ϕ, t, q, p, params)
    ϕ[1] = ...
end

function h(t, q, p, params)
    ...
end
```
where `t` is the current time, `q`, `p`, `λ` and `μ` are the current solution vectors,
`v`, `f`, `u` and `g` are the vectors which hold the result of evaluating the
vector fields ``v`` and ``f``, the projections on the primary constraint ``u`` and ``g``,
`ϕ` holds the algebraic constraint ``\phi``, and `h` returns the Hamiltonian of the system,
all evaluated on `t`, `q`, `p` and `λ`.

Some integrators also enforce the secondary constraint ``\psi`` and require
the following additional functions
```
function ū(u, t, q, p, μ, params)
    u[1] = ...
    u[2] = ...
    ...
end

function ḡ(g, t, q, p, μ, params)
    g[1] = ...
    g[2] = ...
    ...
end

function ψ(ψ, t, q, p, v, f, params)
    ψ[1] = ...
end
```
"""

@doc """
`HDAE`: Hamiltonian Differential Algebraic Equation

$(hdae_equations)

### Parameters

* `vType <: Callable`: type of `v`
* `fType <: Callable`: type of `f`
* `uType <: Callable`: type of `u`
* `gType <: Callable`: type of `g`
* `ϕType <: Callable`: type of `ϕ`
* `ūType <: Callable`: type of `ū`
* `ḡType <: Callable`: type of `ḡ`
* `ψType <: Callable`: type of `ψ`
* `v̄Type <: Callable`: type of `v̄`
* `f̄Type <: Callable`: type of `f̄`
* `hamType <: Callable`: Hamiltonian type
* `invType <: OptionalInvariants`: invariants type
* `parType <: OptionalParameters`: parameters type
* `perType <: OptionalPeriodicity`: periodicity type

### Fields

* `v`: function computing the Hamiltonian vector field ``v``
* `f`: function computing the Hamiltonian vector field ``f``
* `u`: function computing the projection for ``q``
* `g`: function computing the primary projection field ``g``
* `ϕ`: primary constraints
* `ū`: function computing the secondary projection field ``\\bar{u}`` (*optional*)
* `ḡ`: function computing the secondary projection field ``\\bar{g}`` (*optional*)
* `ψ`: secondary constraints (*optional*)
* `v̄`: function computing an initial guess for the velocity field ``v`` (*optional*, defaults to `v`)
* `f̄`: function computing an initial guess for the force field ``f`` (*optional*, defaults to `f`)
* `hamiltonian`: function computing the Hamiltonian ``H`` (usually the total energy of the system)
* `invariants`: functions for the computation of invariants, either a `NamedTuple` containing the equation's invariants or `NullInvariants`
* `parameters`: type constraints for parameters, either a `NamedTuple` containing the equation's parameters or `NullParameters`
* `periodicity`: determines the periodicity of the state vector `q` for cutting periodic solutions, either a `AbstractArray` or `NullPeriodicity`

### Constructors

```julia
HDAE(v, f, u, g, ϕ, ū, ḡ, ψ, h, v̄, f̄, invariants, parameters, periodicity)
HDAE(v, f, u, g, ϕ, ū, ḡ, ψ, h; kwargs...)
HDAE(v, f, u, g, ϕ, h; kwargs...)
```

$(hdae_constructors)

### Function Definitions

The functions are defined by

$(hdae_functions)

The `HDAE` is created by

```julia
equ = HDAE(v, f, u, g, ϕ, h)
```
or
```julia
equ = HDAE(v, f, u, g, ϕ, ū, ḡ, ψ, h)
```
"""
struct HDAE{vType <: Callable, fType <: Callable,
    uType <: Callable, gType <: Callable, ϕType <: Callable,
    ūType <: OptionalCallable, ḡType <: OptionalCallable, ψType <: OptionalCallable,
    v̄Type <: Callable, f̄Type <: Callable,
    hamType <: Callable,
    invType <: OptionalInvariants,
    parType <: OptionalParameters,
    perType <: OptionalPeriodicity} <:
       AbstractEquationPDAE{invType, parType, perType, ψType}
    v::vType
    f::fType
    u::uType
    g::gType
    ϕ::ϕType
    ū::ūType
    ḡ::ḡType
    ψ::ψType
    v̄::v̄Type
    f̄::f̄Type

    hamiltonian::hamType
    invariants::invType
    parameters::parType
    periodicity::perType

    function HDAE(v, f, u, g, ϕ, ū, ḡ, ψ, v̄, f̄, hamiltonian,
            invariants, parameters, periodicity)
        @assert !isempty(methods(v))
        @assert !isempty(methods(f))
        @assert !isempty(methods(u))
        @assert !isempty(methods(g))
        @assert !isempty(methods(ϕ))
        @assert !isempty(methods(ū)) || ū === nothing
        @assert !isempty(methods(ḡ)) || ḡ === nothing
        @assert !isempty(methods(ψ)) || ψ === nothing
        @assert !isempty(methods(v̄))
        @assert !isempty(methods(f̄))
        @assert !isempty(methods(hamiltonian))

        _periodicity = promote_periodicity(periodicity)

        new{typeof(v), typeof(f),
            typeof(u), typeof(g), typeof(ϕ),
            typeof(ū), typeof(ḡ), typeof(ψ),
            typeof(v̄), typeof(f̄),
            typeof(hamiltonian), typeof(invariants), typeof(parameters), typeof(_periodicity)}(
            v, f, u, g, ϕ, ū, ḡ, ψ, v̄, f̄, hamiltonian,
            invariants, parameters, _periodicity)
    end
end

function HDAE(v, f, u, g, ϕ, ū, ḡ, ψ, v̄, f̄, hamiltonian; invariants = NullInvariants(),
        parameters = NullParameters(), periodicity = NullPeriodicity())
    HDAE(v, f, u, g, ϕ, ū, ḡ, ψ, v̄, f̄, hamiltonian, invariants, parameters, periodicity)
end
function HDAE(v, f, u, g, ϕ, ū, ḡ, ψ, hamiltonian; v̄ = v, f̄ = f, kwargs...)
    HDAE(v, f, u, g, ϕ, ū, ḡ, ψ, v̄, f̄, hamiltonian; kwargs...)
end
function HDAE(v, f, u, g, ϕ, hamiltonian; kwargs...)
    HDAE(v, f, u, g, ϕ, nothing, nothing, nothing, hamiltonian; kwargs...)
end

GeometricBase.invariants(equation::HDAE) = equation.invariants
GeometricBase.parameters(equation::HDAE) = equation.parameters
GeometricBase.periodicity(equation::HDAE) = equation.periodicity

hasvectorfield(::HDAE) = true
hashamiltonian(::HDAE) = true
function hasinitialguess(::HDAE{
        vType, fType, uType, gType, ϕType, ūType, ḡType, ψType, <:Callable,
        <:Callable}) where {vType, fType, uType, gType, ϕType, ūType, ḡType, ψType}
    true
end

function Base.show(io::IO, equation::HDAE)
    print(io, "Hamiltonian Differential Algebraic Equation (HDAE)", "\n")
    print(io, "\n")
    print(io, " with vector fields")
    print(io, "\n")
    print(io, "   v = ", equation.v, "\n")
    print(io, "   f = ", equation.f, "\n")
    print(io, "   u = ", equation.u, "\n")
    print(io, "   g = ", equation.g, "\n")
    print(io, "   ū = ", equation.ū, "\n")
    print(io, "   ḡ = ", equation.ḡ, "\n")
    print(io, "\n")
    print(io, " and constraints")
    print(io, "\n")
    print(io, "   ϕ = ", equation.ϕ, "\n")
    print(io, "   ψ = ", equation.ψ, "\n")
    print(io, "\n")
    print(io, " Hamiltonian: H = ", equation.hamiltonian, "\n")
    print(io, "\n")
    print(io, " Invariants: \n")
    print(io, "   ", invariants(equation))
end

function initialstate(
        equ::HDAE, t::InitialTime, ics::NamedTuple, params::OptionalParameters)
    ics = haskey(ics, :μ) ? ics : merge(ics, (μ = zeroalgebraic(ics.λ),))

    (
        q = _statevariable(ics.q, periodicity(equ)),
        p = _statevariable(ics.p),
        λ = _algebraicvariable(ics.λ),
        μ = _algebraicvariable(ics.μ)
    )
end

function initialstate(equ::HDAE, q₀::InitialState, p₀::InitialState,
        λ₀::InitialAlgebraic, μ₀::InitialAlgebraic = zeroalgebraic(λ₀))
    initialstate(equ, (q = q₀, p = p₀, λ = λ₀, μ = μ₀))
end

function initialstate(equ::HDAE, q₀::InitialStateVector, p₀::InitialStateVector,
        λ₀::InitialAlgebraicVector, μ₀::InitialAlgebraicVector = zeroalgebraic(λ₀))
    [initialstate(equ, q, p, λ, μ) for (q, p, λ, μ) in zip(q₀, p₀, λ₀, μ₀)]
end

function check_initial_conditions(equ::HDAE, ics::NamedTuple)
    haskey(ics, :q) || return false
    haskey(ics, :p) || return false
    haskey(ics, :λ) || return false
    if hassecondary(equ)
        haskey(ics, :μ) || return false
        eltype(ics.λ) == eltype(ics.μ) || return false
        typeof(ics.λ) == typeof(ics.μ) || return false
        axes(ics.λ) == axes(ics.μ) || return false
    end
    return true
end

function check_methods(equ::HDAE, timespan, ics::NamedTuple, params)
    applicable(equ.v, zero(ics.q), timespan[begin], ics.q, ics.p, params) || return false
    applicable(equ.f, zero(ics.p), timespan[begin], ics.q, ics.p, params) || return false
    applicable(equ.u, zero(ics.q), timespan[begin], ics.q, ics.p, ics.λ, params) ||
        return false
    applicable(equ.g, zero(ics.p), timespan[begin], ics.q, ics.p, ics.λ, params) ||
        return false
    applicable(equ.ϕ, zero(ics.λ), timespan[begin], ics.q, ics.p, params) || return false
    applicable(equ.hamiltonian, timespan[begin], ics.q, ics.p, params) || return false
    equ.ū === nothing ||
        applicable(equ.ū, zero(ics.q), timespan[begin], ics.q, ics.p, ics.λ, params) ||
        return false
    equ.ḡ === nothing ||
        applicable(equ.ḡ, zero(ics.p), timespan[begin], ics.q, ics.p, ics.λ, params) ||
        return false
    equ.ψ === nothing ||
        applicable(equ.ψ, zero(ics.λ), timespan[begin], ics.q, ics.p,
            vectorfield(ics.q), vectorfield(ics.p), params) || return false
    equ.v̄ === nothing ||
        applicable(equ.v̄, zero(ics.q), timespan[begin], ics.q, ics.p, params) ||
        return false
    equ.f̄ === nothing ||
        applicable(equ.f̄, zero(ics.p), timespan[begin], ics.q, ics.p, params) ||
        return false
    return true
end

function GeometricBase.datatype(equ::HDAE, ics::NamedTuple)
    @assert check_initial_conditions(equ, ics)
    return eltype(ics.q)
end

function GeometricBase.arrtype(equ::HDAE, ics::NamedTuple)
    @assert check_initial_conditions(equ, ics)
    return typeof(ics.q)
end

_get_v(equ::HDAE, params) = (v, t, q, p) -> equ.v(v, t, q, p, params)
_get_f(equ::HDAE, params) = (f, t, q, p) -> equ.f(f, t, q, p, params)
_get_u(equ::HDAE, params) = (u, t, q, p, λ) -> equ.u(u, t, q, p, λ, params)
_get_g(equ::HDAE, params) = (g, t, q, p, λ) -> equ.g(g, t, q, p, λ, params)
_get_ϕ(equ::HDAE, params) = (ϕ, t, q, p) -> equ.ϕ(ϕ, t, q, p, params)
_get_ū(equ::HDAE, params) = (u, t, q, p, λ) -> equ.ū(u, t, q, p, λ, params)
_get_ḡ(equ::HDAE, params) = (g, t, q, p, λ) -> equ.ḡ(g, t, q, p, λ, params)
_get_ψ(equ::HDAE, params) = (ψ, t, q, p, v, f) -> equ.ψ(ψ, t, q, p, v, f, params)
_get_v̄(equ::HDAE, params) = (v, t, q, p) -> equ.v̄(v, t, q, p, params)
_get_f̄(equ::HDAE, params) = (f, t, q, p) -> equ.f̄(f, t, q, p, params)
_get_h(equ::HDAE, params) = (p, t, q) -> equ.hamiltonian(t, q, p, params)
_get_invariant(::HDAE, inv, params) = (t, q, p) -> inv(t, q, p, params)

function _functions(equ::HDAE)
    if hassecondary(equ)
        (v = equ.v, f = equ.f, u = equ.u, g = equ.g, ϕ = equ.ϕ,
            ū = equ.ū, ḡ = equ.ḡ, ψ = equ.ψ, h = equ.hamiltonian)
    else
        (v = equ.v, f = equ.f, u = equ.u, g = equ.g, ϕ = equ.ϕ, h = equ.hamiltonian)
    end
end

function _functions(equ::HDAE, params::OptionalParameters)
    if hassecondary(equ)
        (
            v = _get_v(equ, params),
            f = _get_f(equ, params),
            u = _get_u(equ, params),
            g = _get_g(equ, params),
            ϕ = _get_ϕ(equ, params),
            ū = _get_ū(equ, params),
            ḡ = _get_ḡ(equ, params),
            ψ = _get_ψ(equ, params),
            h = _get_h(equ, params)
        )
    else
        (
            v = _get_v(equ, params),
            f = _get_f(equ, params),
            u = _get_u(equ, params),
            g = _get_g(equ, params),
            ϕ = _get_ϕ(equ, params),
            h = _get_h(equ, params)
        )
    end
end

_initialguess(equ::HDAE) = (v = equ.v̄, f = equ.f̄)
function _initialguess(equ::HDAE, params::OptionalParameters)
    (v = _get_v̄(equ, params), f = _get_f̄(equ, params))
end

@doc """
`HDAEProblem`: Hamiltonian Differential Algebraic Equation

$(hdae_equations)

The dynamical variables ``(q,p)`` with initial conditions ``(q(t_{0}) = q_{0}, p(t_{0}) = p_{0})``
take values in ``\\mathbb{R}^{d} \\times \\mathbb{R}^{d}``. The algebraic variables ``(λ,μ)``
with initial condition ``(λ(t_{0}) = λ_{0}, μ(t_{0}) = μ_{0})`` take values in ``\\mathbb{R}^{m} \\times \\mathbb{R}^{m}``.

### Constructors

```julia
HDAEProblem(v, f, u, g, ϕ, ū, ḡ, ψ, h, timespan, timestep, ics::NamedTuple; kwargs...)
HDAEProblem(v, f, u, g, ϕ, ū, ḡ, ψ, h, timespan, timestep, q₀::StateVariable, p₀::StateVariable, λ₀::StateVariable, μ₀::StateVariable = zero(λ₀); kwargs...)
HDAEProblem(v, f, u, g, ϕ, h, timespan, timestep, ics::NamedTuple; kwargs...)
HDAEProblem(v, f, u, g, ϕ, h, timespan, timestep, q₀::StateVariable, p₀::StateVariable, λ₀::StateVariable; kwargs...)
```

$(hdae_constructors)

`timespan` is the time interval `(t₀,t₁)` for the problem to be solved in,
`timestep` is the time step to be used in the simulation, and
`ics` is a `NamedTuple` with entries `q`, `p`, `λ` and `μ`.
The initial conditions `q₀`, `p₀`, `λ₀` and `μ₀` can also be prescribed directly,
with `StateVariable` an `AbstractArray{<:Number}`.
For the interfaces of the functions `v`, `f`, `u`, `g`, `ϕ`, `ū`, `ḡ`, `ψ`, and `h` see [`HDAE`](@ref).

In addition to the standard keyword arguments for [`EquationProblem`](@ref GeometricEquations.EquationProblem) subtypes,
a `HDAEProblem` accepts functions `v̄` and `f̄` for the computation of initial guesses for the vector fields with default
values `v̄ = v` and `f̄ = f`.

### Function Definitions

$(hdae_functions)

With the above function definitions the `HDAEProblem` can be created by

```julia
timespan = (0.0, 1.0)
timestep = 0.1
q₀ = [1., 1.]
p₀ = [1., 0.]
λ₀ = [0.]
μ₀ = [0.]

prob = HDAEProblem(v, f, u, g, ϕ, h, timespan, timestep, q₀, p₀, λ₀)
```
or
```julia
prob = HDAEProblem(v, f, u, g, ϕ, ū, ḡ, ψ, h, timespan, timestep, q₀, p₀, λ₀, μ₀)
```
"""
const HDAEProblem = EquationProblem{HDAE}

function HDAEProblem(
        v, f, u, g, ϕ, ū, ḡ, ψ, hamiltonian, timespan::Tuple, timestep::Real, ics...;
        v̄ = v, f̄ = f, invariants = NullInvariants(),
        parameters = NullParameters(), periodicity = NullPeriodicity())
    equ = HDAE(v, f, u, g, ϕ, ū, ḡ, ψ, v̄, f̄, hamiltonian, invariants,
        parameter_types(parameters), periodicity)
    EquationProblem(equ, timespan, timestep, initialstate(equ, ics...), parameters)
end

function HDAEProblem(v, f, u, g, ϕ, hamiltonian, args...; kwargs...)
    HDAEProblem(v, f, u, g, ϕ, nothing, nothing, nothing, hamiltonian, args...; kwargs...)
end

function GeometricBase.periodicity(prob::HDAEProblem)
    (q = periodicity(equation(prob)), p = NullPeriodicity(), λ = NullPeriodicity(),
        μ = NullPeriodicity())
end

function compute_vectorfields!(vecfield, sol, prob::HDAEProblem)
    initialguess(prob).v(vecfield.q, sol.t, sol.q, sol.p, parameters(prob))
    initialguess(prob).f(vecfield.p, sol.t, sol.q, sol.p, parameters(prob))
end

const HDAEEnsemble = EnsembleProblem{HDAE}

function HDAEEnsemble(v, f, u, g, ϕ, ū, ḡ, ψ, hamiltonian, timespan::Tuple,
        timestep::Real, ics...; v̄ = v, f̄ = f,
        invariants = NullInvariants(),
        parameters = NullParameters(),
        periodicity = NullPeriodicity())
    equ = HDAE(v, f, u, g, ϕ, ū, ḡ, ψ, v̄, f̄, hamiltonian,
        invariants, parameter_types(parameters), periodicity)
    EnsembleProblem(equ, timespan, timestep, initialstate(equ, ics...), parameters)
end

function HDAEEnsemble(v, f, u, g, ϕ, hamiltonian, args...; kwargs...)
    HDAEEnsemble(v, f, u, g, ϕ, nothing, nothing, nothing, hamiltonian, args...; kwargs...)
end
