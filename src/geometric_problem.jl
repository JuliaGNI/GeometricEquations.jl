@doc raw"""
GeometricProblem: stores a GeometricEquation togehter with initial conditions, parameters, time span and time step size.

### Parameters

* `ST <: GeometricEquation`: super type, used for dispatch
* `DT <: Number`: data type
* `TT <: Real`: time step type
* `AT <: AbstractArray{DT}`: array type
* `equType <: GeometricEquation`: equation type
* `icsType <: paramsType`: initial conditions type
* `parType <: OptionalParameters`: parameters type

### Fields

* `equation`: reference to the parent equation object
* `tspan`: time span for problem `(t₀,t₁)`
* `tstep`: time step to be used in simulation
* `ics`: `NamedTuple` containing the initial conditions
* `parameters`: either a `NamedTuple` containing the equation's parameters or `NullParameters`

### Subtypes

The `GeometricProblem` type has various subtypes for the different equations types, that are defined e.g. via
```
const ODEProblem = GeometricProblem{ODE}
```
and provide convenience constructors to construct and equation and the corresponding problem in one step.

"""
struct GeometricProblem{superType<:GeometricEquation,dType<:Number,tType<:Real,arrayType<:AbstractArray{dType},
    equType<:GeometricEquation,icsType<:NamedTuple,paramsType<:OptionalParameters} <: AbstractProblem{dType,tType,arrayType}
    equation::equType
    tspan::Tuple{tType,tType}
    tstep::tType
    ics::icsType
    parameters::paramsType

    function GeometricProblem(equ, tspan, tstep, ics, parameters)
        _tspan = promote_tspan(tspan)
        _tspan, _tstep = promote_tspan_and_tstep(_tspan, tstep)

        @assert check_initial_conditions(equ, ics)
        @assert check_methods(equ, _tspan, ics, parameters)

        superType = eval(typeof(equ).name.name)
        tType = typeof(_tstep)
        dType = datatype(equ, ics)
        arrayType = arrtype(equ, ics)

        new{superType,dType,tType,arrayType,typeof(equ),typeof(ics),typeof(parameters)}(equ, _tspan, _tstep, ics, parameters)
    end
end

GeometricProblem(equ, tspan, tstep, ics; parameters = NullParameters()) = GeometricProblem(equ, tspan, tstep, ics, parameters)
GeometricProblem(equ, tspan, tstep, ics, ::Nothing) = GeometricProblem(equ, tspan, tstep, ics, NullParameters())

# Base.hash(prob::GeometricProblem, h::UInt) = hash(hash(prob.equation,
#                                                   hash(prob.tspan, hash(prob.tstep,
#                                                   hash(prov.ics, hash(prov.parameters, h))))))

# Base.:(==)(prob1::GeometricProblem, prob2::GeometricProblem) = (
#                                 prob1.equation   == prob2.equation
#                              && prob1.tspan      == prob2.tspan
#                              && prob1.tstep      == prob2.tstep
#                              && prob1.ics        == prob2.ics
#                              && prob1.parameters == prob2.parameters)

@inline GeometricBase.datatype(::GeometricProblem{ST,DT,TT,AT}) where {ST,DT,TT,AT} = DT
@inline GeometricBase.timetype(::GeometricProblem{ST,DT,TT,AT}) where {ST,DT,TT,AT} = TT
@inline GeometricBase.arrtype(::GeometricProblem{ST,DT,TT,AT}) where {ST,DT,TT,AT} = AT
@inline GeometricBase.equtype(::GeometricProblem{ST,DT,TT,AT}) where {ST,DT,TT,AT} = ST

@inline GeometricBase.equation(prob::GeometricProblem) = prob.equation
@inline GeometricBase.tspan(prob::GeometricProblem) = prob.tspan
@inline GeometricBase.tstep(prob::GeometricProblem) = prob.tstep

@inline GeometricBase.timestep(prob::GeometricProblem) = tstep(prob)
@inline GeometricBase.parameters(prob::GeometricProblem) = prob.parameters
@inline GeometricBase.nsamples(::GeometricProblem) = 1
@inline GeometricBase.ntime(prob::GeometricProblem) = Int(abs(div(tend(prob) - tbegin(prob), tstep(prob), RoundUp)))

@inline GeometricBase.functions(prob::GeometricProblem) = functions(equation(prob), parameters(prob))
@inline GeometricBase.solutions(prob::GeometricProblem) = solutions(equation(prob), parameters(prob))
@inline GeometricBase.invariants(prob::GeometricProblem) = invariants(equation(prob), parameters(prob))

initial_conditions(prob::GeometricProblem) = merge( (t = tbegin(prob),), prob.ics )

function Base.similar(prob::GeometricProblem, tspan, tstep = tstep(prob), ics = prob.ics, parameters = parameters(prob))
    GeometricProblem(equation(prob), tspan, tstep, ics, parameters)
end

function Base.similar(prob::GeometricProblem; tspan = tspan(prob), tstep = tstep(prob), ics = prob.ics, parameters = parameters(prob))
    similar(prob, tspan, tstep, ics, parameters)
end


@doc raw"""
`ODEProblem`: Ordinary Differential Equation Problem

Defines an initial value problem
```math
\dot{q} (t) = v(t, q(t)) , \qquad q(t_{0}) = q_{0} ,
```
with vector field ``v`` and initial condition ``q_{0}``.

### Constructors

```julia
ODE(v, tspan, tstep, ics::NamedTuple; kwargs...)
ODE(v, tspan, tstep, q₀::State; kwargs...)
```
where `ics` is a `NamedTuple` with entry `q`.

### Keyword arguments:

* `invariants = NullInvariants()`
* `parameters = NullParameters()`
* `periodicity = NullPeriodicity()`

"""
const ODEProblem = GeometricProblem{ODE}

function ODEProblem(v, tspan, tstep, ics::NamedTuple; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    equ = ODE(v, invariants, parameter_types(parameters), periodicity)
    GeometricProblem(equ, tspan, tstep, ics, parameters)
end

function ODEProblem(v, tspan, tstep, q₀::State; kwargs...)
    ics = (q = q₀,)
    ODEProblem(v, tspan, tstep, ics; kwargs...)
end

GeometricBase.periodicity(prob::ODEProblem) = (q = periodicity(equation(prob)),)


@doc raw"""
`PODEProblem`: Partitioned Ordinary Differential Equation Problem

Defines a partitioned initial value problem
```math
\begin{aligned}
\dot{q} (t) &= v(t, q(t), p(t)) , &
q(t_{0}) &= q_{0} , \\
\dot{p} (t) &= f(t, q(t), p(t)) , &
p(t_{0}) &= p_{0} ,
\end{aligned}
```
with vector fields ``v`` and ``f`` and initial conditions ``(q_{0}, p_{0})``.


### Constructors

```julia
PODE(v, f, tspan, tstep, ics; kwargs...)
PODE(v, f, tspan, tstep, q₀::State, p₀::State; kwargs...)
```
where `ics` is a `NamedTuple` with entries `q` and `p`.

### Keyword arguments:

* `invariants = NullInvariants()`
* `parameters = NullParameters()`
* `periodicity = NullPeriodicity()`

"""
const PODEProblem = GeometricProblem{PODE}

function PODEProblem(v, f, tspan, tstep, ics::NamedTuple; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    equ = PODE(v, f, invariants, parameter_types(parameters), periodicity)
    GeometricProblem(equ, tspan, tstep, ics, parameters)
end

function PODEProblem(v, f, tspan, tstep, q₀::State, p₀::State; kwargs...)
    ics = (q = q₀, p = p₀)
    PODEProblem(v, f, tspan, tstep, ics; kwargs...)
end

GeometricBase.periodicity(prob::PODEProblem) = (q = periodicity(equation(prob)), p = NullPeriodicity())


@doc raw"""
`HODEProblem`: Hamiltonian Ordinary Differential Equation Problem

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

### Constructors

```julia
HODEProblem(v, f, poisson, hamiltonian, tspan, tstep, ics; kwargs...)
HODEProblem(v, f, poisson, hamiltonian, tspan, tstep, q₀::State, p₀::State; kwargs...)
```
where `ics` is a `NamedTuple` with entries `q` and `p`.

### Keyword arguments:

* `invariants = NullInvariants()`
* `parameters = NullParameters()`
* `periodicity = NullPeriodicity()`

"""
const HODEProblem = GeometricProblem{HODE}

function HODEProblem(v, f, hamiltonian, tspan, tstep, ics::NamedTuple; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    equ = HODE(v, f, hamiltonian, invariants, parameter_types(parameters), periodicity)
    GeometricProblem(equ, tspan, tstep, ics, parameters)
end

function HODEProblem(v, f, hamiltonian, tspan, tstep, q₀::State, p₀::State; kwargs...)
    ics = (q = q₀, p = p₀)
    HODEProblem(v, f, hamiltonian, tspan, tstep, ics; kwargs...)
end

GeometricBase.periodicity(prob::HODEProblem) = (q = periodicity(equation(prob)), p = NullPeriodicity())


@doc raw"""
`IODEProblem`: Implicit Ordinary Differential Equation Problem

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
with momentum defined by ``\vartheta``, force field ``f``, projection field ``g``
and initial conditions ``(q_{0}, p_{0}, \lambda_{0})``.
This is a special case of a differential algebraic equation with dynamical
variables ``(q,p)`` and algebraic variable ``v``
"""
const IODEProblem = GeometricProblem{IODE}

function IODEProblem(ϑ, f, g, tspan, tstep, ics::NamedTuple; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity(), v̄ = _iode_default_v̄, f̄ = f)
    equ = IODE(ϑ, f, g, v̄, f̄, invariants, parameter_types(parameters), periodicity)
    GeometricProblem(equ, tspan, tstep, ics, parameters)
end

function IODEProblem(ϑ, f, g, tspan, tstep, q₀::State, p₀::State, λ₀::State = zero(q₀); kwargs...)
    ics = (q = q₀, p = p₀, λ = λ₀)
    IODEProblem(ϑ, f, g, tspan, tstep, ics; kwargs...)
end

GeometricBase.periodicity(prob::IODEProblem) = (q = periodicity(equation(prob)), p = NullPeriodicity(), λ = NullPeriodicity())


@doc raw"""
`LODEProblem`: Lagrangian Ordinary Differential Equation Problem

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
with momentum ``p`` and force field ``f``, given by
```math
\begin{aligned}
p &= \frac{\partial L}{\partial v} , &
f &= \frac{\partial L}{\partial q} ,
\end{aligned}
```
and initial conditions ``(q_{0}, p_{0})``.
This is a special case of a differential algebraic equation with dynamical
variables ``(q,p)`` and algebraic variable ``v``.
"""
const LODEProblem = GeometricProblem{LODE}

function LODEProblem(ϑ, f, g, ω, l, tspan, tstep, ics::NamedTuple; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity(), v̄ = _lode_default_v̄, f̄ = f)
    equ = LODE(ϑ, f, g, ω, v̄, f̄, l, invariants, parameter_types(parameters), periodicity)
    GeometricProblem(equ, tspan, tstep, ics, parameters)
end

function LODEProblem(ϑ, f, g, ω, l, tspan, tstep, q₀::State, p₀::State, λ₀::State = zero(q₀); kwargs...)
    ics = (q = q₀, p = p₀, λ = λ₀)
    LODEProblem(ϑ, f, g, ω, l, tspan, tstep, ics; kwargs...)
end

GeometricBase.periodicity(prob::LODEProblem) = (q = periodicity(equation(prob)), p = NullPeriodicity(), λ = NullPeriodicity())


@doc raw"""
`SODEProblem`: Split Ordinary Differential Equation Problem

Defines an initial value problem
```math
\dot{q} (t) = v(t, q(t)) , \qquad q(t_{0}) = q_{0} ,
```
with vector field ``v`` and initial condition ``q_{0}``. Here, the vector field
``v`` is given as a sum of vector fields
```math
v (t) = v_1 (t) + ... + v_r (t) .
```

### Fields

* `equation`: reference to the parent equation object
* `functions`: all functions (e.g. vector fields) provided by the equation
* `tspan`: time span for problem `(t₀,t₁)`
* `tstep`: time step to be used in simulation
* `ics`: initial condition (NamedTuple containing the field `q`)
* `parameters`: either a `NamedTuple` containing the equation's parameters or `nothing`

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
    function q_i(t, q₀, q₁, h, params))
        q₁[1] = q₀[1] + ...
        q₁[2] = q₀[2] + ...
        ...
    end
```
where `t` is the current time, `q₀` is the current solution vector, `q₁` is the
new solution vector which holds the result of computing one substep with the
vector field ``v_i`` on `t` and `q₀`, and `h` is the (sub-)timestep to compute
the update for.

The fact that the function `v` returns the solution and not just the vector
field for each substep increases the flexibility for the use of splitting
methods, e.g., it allows to use another integrator for solving substeps.

### Constructors

```julia
SODE(v, q, t₀, q₀, invariants, parameters, periodicity)

SODE(v, q::Union{Tuple,Nothing}, t₀::Real, q₀::StateVector; kwargs...)
SODE(v, q::Union{Tuple,Nothing}, t₀::Real, q₀::State; kwargs...)
SODE(v, q::Union{Tuple,Nothing}, q₀::StateVector; kwargs...)
SODE(v, q::Union{Tuple,Nothing}, q₀::State; kwargs...)

SODE(v, t₀::Real, q₀::StateVector; kwargs...)
SODE(v, t₀::Real, q₀::State; kwargs...)
SODE(v, q₀::StateVector; kwargs...)
SODE(v, q₀::State; kwargs...)
```

"""
const SODEProblem = GeometricProblem{SODE}

function SODEProblem(v::Tuple, q::Union{Tuple,Nothing}, tspan, tstep, ics::NamedTuple; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    equ = SODE(v, q, invariants, parameter_types(parameters), periodicity)
    GeometricProblem(equ, tspan, tstep, ics, parameters)
end

function SODEProblem(v::Tuple, q::Union{Tuple,Nothing}, tspan, tstep, q₀::State; kwargs...)
    ics = (q = q₀,)
    SODEProblem(v, q, tspan, tstep, ics; kwargs...)
end

function SODEProblem(v, tspan, tstep, ics::NamedTuple; kwargs...)
    SODEProblem(v, nothing, tspan, tstep, ics; kwargs...)
end

function SODEProblem(v, tspan, tstep, q₀::State; kwargs...)
    ics = (q = q₀,)
    SODEProblem(v, tspan, tstep, ics; kwargs...)
end

GeometricBase.periodicity(prob::SODEProblem) = (q = periodicity(equation(prob)), )


@doc raw"""
`DAEProblem`: Differential Algebraic Equation Problem

Defines a differential algebraic initial value problem
```math
\begin{aligned}
\dot{q} (t) &= v(t, q(t)) + u(t, q(t), \lambda(t)) , & q(t_{0}) &= q_{0} , \\
0 &= \phi (t, q(t), \lambda(t)) , & \lambda(t_{0}) &= \lambda_{0} ,
\end{aligned}
```
with vector field ``v``, projection ``u``, algebraic constraint ``\phi=0``,
and initial conditions ``q_{0}`` and ``\lambda_{0}``.
"""
const DAEProblem = GeometricProblem{DAE}

function DAEProblem(v, u, ϕ, ū, ψ, tspan, tstep, ics::NamedTuple; v̄ = v, invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    equ = DAE(v, u, ϕ, ū, ψ, v̄, invariants, parameter_types(parameters), periodicity)
    GeometricProblem(equ, tspan, tstep, ics, parameters)
end

function DAEProblem(v, u, ϕ, ū, ψ, tspan, tstep, q₀::State, λ₀::State, μ₀::State = zero(λ₀); kwargs...)
    ics = (q = q₀, λ = λ₀, μ = μ₀)
    DAEProblem(v, u, ϕ, ū, ψ, tspan, tstep, ics; kwargs...)
end

function DAEProblem(v, u, ϕ, tspan, tstep, ics::NamedTuple; kwargs...)
    DAEProblem(v, u, ϕ, nothing, nothing, tspan, tstep, ics; kwargs...)
end

function DAEProblem(v, u, ϕ, tspan, tstep, q₀::State, λ₀::State; kwargs...)
    ics = (q = q₀, λ = λ₀)
    DAEProblem(v, u, ϕ, tspan, tstep, ics; kwargs...)
end

GeometricBase.periodicity(prob::DAEProblem) = (q = periodicity(equation(prob)), λ = NullPeriodicity(), μ = NullPeriodicity())


@doc raw"""
`PDAEProblem`: Partitioned Differential Algebraic Equation Problem

Defines a partitioned differential algebraic initial value problem
```math
\begin{aligned}
\dot{q} (t) &= v(t, q(t), p(t)) + u(t, q(t), p(t), \lambda(t)) , & q(t_{0}) &= q_{0} , \\
\dot{p} (t) &= f(t, q(t), p(t)) + g(t, q(t), p(t), \lambda(t)) , & p(t_{0}) &= p_{0} , \\
0 &= \phi (t, q(t), p(t), \lambda(t)) , & \lambda(t_{0}) &= \lambda_{0} ,
\end{aligned}
```
with vector fields ``v`` and ``f``, projection ``u`` and ``g``,
algebraic constraint ``\phi=0``, and initial conditions
``(q_{0}, p_{0})`` and ``\lambda_{0}``.
"""
const PDAEProblem = GeometricProblem{PDAE}

function PDAEProblem(v, f, u, g, ϕ, ū, ḡ, ψ, tspan, tstep, ics::NamedTuple; v̄ = v, f̄ = f, invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    equ = PDAE(v, f, u, g, ϕ, ū, ḡ, ψ, v̄, f̄, invariants, parameter_types(parameters), periodicity)
    GeometricProblem(equ, tspan, tstep, ics, parameters)
end

function PDAEProblem(v, f, u, g, ϕ, ū, ḡ, ψ, tspan, tstep, q₀::State, p₀::State, λ₀::State, μ₀::State = zero(λ₀); kwargs...)
    ics = (q = q₀, p = p₀, λ = λ₀, μ = μ₀)
    PDAEProblem(v, f, u, g, ϕ, ū, ḡ, ψ, tspan, tstep, ics; kwargs...)
end

function PDAEProblem(v, f, u, g, ϕ, tspan, tstep, ics::NamedTuple; kwargs...)
    PDAEProblem(v, f, u, g, ϕ, nothing, nothing, nothing, tspan, tstep, ics; kwargs...)
end

function PDAEProblem(v, f, u, g, ϕ, tspan, tstep, q₀::State, p₀::State, λ₀::State; kwargs...)
    ics = (q = q₀, p = p₀, λ = λ₀)
    PDAEProblem(v, f, u, g, ϕ, tspan, tstep, ics; kwargs...)
end

GeometricBase.periodicity(prob::PDAEProblem) = (q = periodicity(equation(prob)), p = NullPeriodicity(), λ = NullPeriodicity())


@doc raw"""
`HDAE`: Hamiltonian Differential Algebraic Equation

Defines a Hamiltonian differential algebraic initial value problem, that is
a canonical Hamiltonian system of equations subject to Dirac constraints,
```math
\begin{aligned}
\dot{q} (t) &= v(t, q(t), p(t)) + u(t, q(t), p(t), \lambda(t)) + \bar{g}(t, q(t), p(t), \lambda(t), \gamma(t)) , & q(t_{0}) &= q_{0} , \\
\dot{p} (t) &= f(t, q(t), p(t)) + g(t, q(t), p(t), \lambda(t)) + \bar{f}(t, q(t), p(t), \lambda(t), \gamma(t)) , & p(t_{0}) &= p_{0} , \\
0 &= \phi (t, q(t), p(t)) , \\
0 &= \psi (t, q(t), p(t), \dot{q}(t), \dot{p}(t)) ,
\end{aligned}
```
with vector fields ``v``, ``u``, ``\bar{u}`` and ``f``, ``g``, ``\bar{g}``,
primary constraint ``\phi(q,p)=0`` and secondary constraint ``\psi(q,p,\lambda)=0``,
and initial conditions ``(q_{0}, p_{0})``.
"""
const HDAEProblem = GeometricProblem{HDAE}

function HDAEProblem(v, f, u, g, ϕ, ū, ḡ, ψ, hamiltonian, tspan, tstep, ics::NamedTuple; v̄ = v, f̄ = f, invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    equ = HDAE(v, f, u, g, ϕ, ū, ḡ, ψ, v̄, f̄, hamiltonian, invariants, parameter_types(parameters), periodicity)
    GeometricProblem(equ, tspan, tstep, ics, parameters)
end

function HDAEProblem(v, f, u, g, ϕ, ū, ḡ, ψ, hamiltonian, tspan, tstep, q₀::State, p₀::State, λ₀::State, μ₀::State = zero(λ₀); kwargs...)
    ics = (q = q₀, p = p₀, λ = λ₀, μ = μ₀)
    HDAEProblem(v, f, u, g, ϕ, ū, ḡ, ψ, hamiltonian, tspan, tstep, ics; kwargs...)
end

function HDAEProblem(v, f, u, g, ϕ, hamiltonian, tspan, tstep, ics::NamedTuple; kwargs...)
    HDAEProblem(v, f, u, g, ϕ, nothing, nothing, nothing, hamiltonian, tspan, tstep, ics; kwargs...)
end

function HDAEProblem(v, f, u, g, ϕ, hamiltonian, tspan, tstep, q₀::State, p₀::State, λ₀::State; kwargs...)
    ics = (q = q₀, p = p₀, λ = λ₀)
    HDAEProblem(v, f, u, g, ϕ, hamiltonian, tspan, tstep, ics; kwargs...)
end

GeometricBase.periodicity(prob::HDAEProblem) = (q = periodicity(equation(prob)), p = NullPeriodicity(), λ = NullPeriodicity(), μ = NullPeriodicity())


@doc raw"""
`IDAEProblem`: Implicit Differential Algebraic Equation Problem

Defines an implicit differential algebraic initial value problem
```math
\begin{aligned}
\dot{q} (t) &= v(t) + u(t, q(t), p(t), \lambda(t)) , & q(t_{0}) &= q_{0} , \\
\dot{p} (t) &= f(t, q(t), v(t)) + g(t, q(t), p(t), \lambda(t)) , & p(t_{0}) &= p_{0} , \\
p(t) &= \vartheta(t, q(t), v(t)) , && \\
0 &= \phi (t, q(t), p(t)) , & \lambda(t_{0}) &= \lambda_{0} ,
\end{aligned}
```
with force field ``f``, the momentum defined by ``p``, projection ``u`` and ``g``,
algebraic constraint ``\phi=0``, and initial conditions
``(q_{0}, p_{0})`` and ``\lambda_{0}``.
"""
const IDAEProblem = GeometricProblem{IDAE}

function IDAEProblem(ϑ, f, u, g, ϕ, ū, ḡ, ψ, tspan, tstep, ics::NamedTuple; v̄ = _idae_default_v̄, f̄ = f, invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    equ = IDAE(ϑ, f, u, g, ϕ, ū, ḡ, ψ, v̄, f̄, invariants, parameter_types(parameters), periodicity)
    GeometricProblem(equ, tspan, tstep, ics, parameters)
end

function IDAEProblem(ϑ, f, u, g, ϕ, ū, ḡ, ψ, tspan, tstep, q₀::State, p₀::State, λ₀::State, μ₀::State = zero(λ₀); kwargs...)
    ics = (q = q₀, p = p₀, λ = λ₀, μ = μ₀)
    IDAEProblem(ϑ, f, u, g, ϕ, ū, ḡ, ψ, tspan, tstep, ics; kwargs...)
end

function IDAEProblem(ϑ, f, u, g, ϕ, tspan, tstep, ics::NamedTuple; kwargs...)
    IDAEProblem(ϑ, f, u, g, ϕ, nothing, nothing, nothing, tspan, tstep, ics; kwargs...)
end

function IDAEProblem(ϑ, f, u, g, ϕ, tspan, tstep, q₀::State, p₀::State, λ₀::State; kwargs...)
    ics = (q = q₀, p = p₀, λ = λ₀)
    IDAEProblem(ϑ, f, u, g, ϕ, tspan, tstep, ics; kwargs...)
end

GeometricBase.periodicity(prob::IDAEProblem) = (q = periodicity(equation(prob)), p = NullPeriodicity(), λ = NullPeriodicity(), μ = NullPeriodicity())


@doc raw"""
`LDAEProblem`: Lagrangian Differential Algebraic Equation Problem

Defines a Lagrangian differential algebraic initial value problem, that is
a special implicit initial value problem
```math
\begin{aligned}
\dot{q} (t) &= v(t) + \lambda(t), &
q(t_{0}) &= q_{0} , \\
\dot{p} (t) &= f(t, q(t), v(t)) + g(t, q(t), \lambda(t)) + \bar{g} (t, q(t), \mu(t)) , &
p(t_{0}) &= p_{0} , \\
p(t) &= ϑ(t, q(t), v(t)) , \\
0 &= \phi (t, q(t), p(t)) , \\
0 &= \psi (t, q(t), p(t), \dot{q}(t), \dot{p}(t)) ,
\end{aligned}
```
with vector field ``f``, the momentum defined by ``p``, and initial conditions
``(q_{0}, p_{0}, \lambda_{0})``.
This is a special case of a differential algebraic equation with dynamical
variables ``(q,p)`` and algebraic variables ``v``, ``\lambda`` and ``\mu``.
"""
const LDAEProblem = GeometricProblem{LDAE}

function LDAEProblem(ϑ, f, u, g, ϕ, ū, ḡ, ψ, ω, lagrangian, tspan, tstep, ics::NamedTuple; v̄ = _ldae_default_v̄, f̄ = f, invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    equ = LDAE(ϑ, f, u, g, ϕ, ū, ḡ, ψ, ω, v̄, f̄, lagrangian, invariants, parameter_types(parameters), periodicity)
    GeometricProblem(equ, tspan, tstep, ics, parameters)
end

function LDAEProblem(ϑ, f, u, g, ϕ, ū, ḡ, ψ, ω, lagrangian, tspan, tstep, q₀::State, p₀::State, λ₀::State, μ₀::State = zero(λ₀); kwargs...)
    ics = (q = q₀, p = p₀, λ = λ₀, μ = μ₀)
    LDAEProblem(ϑ, f, u, g, ϕ, ū, ḡ, ψ, ω, lagrangian, tspan, tstep, ics; kwargs...)
end

function LDAEProblem(ϑ, f, u, g, ϕ, ω, lagrangian, tspan, tstep, ics::NamedTuple; kwargs...)
    LDAEProblem(ϑ, f, u, g, ϕ, nothing, nothing, nothing, ω, lagrangian, tspan, tstep, ics; kwargs...)
end

function LDAEProblem(ϑ, f, u, g, ϕ, ω, lagrangian, tspan, tstep, q₀::State, p₀::State, λ₀::State; kwargs...)
    ics = (q = q₀, p = p₀, λ = λ₀)
    LDAEProblem(ϑ, f, u, g, ϕ, ω, lagrangian, tspan, tstep, ics; kwargs...)
end

GeometricBase.periodicity(prob::LDAEProblem) = (q = periodicity(equation(prob)), p = NullPeriodicity(), λ = NullPeriodicity(), μ = NullPeriodicity())


"""

"""
const SPDAEProblem = GeometricProblem{SPDAE}

GeometricBase.periodicity(prob::SPDAEProblem) = (q = periodicity(equation(prob)), p = NullPeriodicity(), λ = NullPeriodicity())


"""

"""
const SDEProblem = GeometricProblem{SDE}

function SDEProblem(v, B, tspan, tstep, ics::NamedTuple; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    equ = SDE(v, B, invariants, parameter_types(parameters), periodicity)
    GeometricProblem(equ, tspan, tstep, ics, parameters)
end

function SDEProblem(v, B, tspan, tstep, q₀::State; kwargs...)
    ics = (q = q₀,)
    SDEProblem(v, B, tspan, tstep, ics; kwargs...)
end

GeometricBase.periodicity(prob::SDEProblem) = (q = periodicity(equation(prob)), )


"""

"""
const PSDEProblem = GeometricProblem{PSDE}

function PSDEProblem(v, f, B, G, tspan, tstep, ics::NamedTuple; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    equ = PSDE(v, f, B, G, invariants, parameter_types(parameters), periodicity)
    GeometricProblem(equ, tspan, tstep, ics, parameters)
end

function PSDEProblem(v, f, B, G, tspan, tstep, q₀::State, p₀::State; kwargs...)
    ics = (q = q₀, p = p₀)
    PSDEProblem(v, f, B, G, tspan, tstep, ics; kwargs...)
end

GeometricBase.periodicity(prob::PSDEProblem) = (q = periodicity(equation(prob)), p = NullPeriodicity())


"""

"""
const SPSDEProblem = GeometricProblem{SPSDE}

function SPSDEProblem(v, f1, f2, B, G1, G2, tspan, tstep, ics::NamedTuple; invariants = NullInvariants(), parameters = NullParameters(), periodicity = NullPeriodicity())
    equ = SPSDE(v, f1, f2, B, G1, G2, invariants, parameter_types(parameters), periodicity)
    GeometricProblem(equ, tspan, tstep, ics, parameters)
end

function SPSDEProblem(v, f1, f2, B, G1, G2, tspan, tstep, q₀::State, p₀::State; kwargs...)
    ics = (q = q₀, p = p₀)
    SPSDEProblem(v, f1, f2, B, G1, G2, tspan, tstep, ics; kwargs...)
end

GeometricBase.periodicity(prob::SPSDEProblem) = (q = periodicity(equation(prob)), p = NullPeriodicity())


# Union types for problems of similar kind

const AbstractProblemODE = Union{ODEProblem,SODEProblem}
const AbstractProblemDAE = Union{DAEProblem}
const AbstractProblemSDE = Union{SDEProblem}
const AbstractProblemPODE = Union{PODEProblem,HODEProblem,IODEProblem,LODEProblem}
const AbstractProblemPDAE = Union{PDAEProblem,HDAEProblem,IDAEProblem,LDAEProblem}
const AbstractProblemPSDE = Union{PSDEProblem,SPSDEProblem}
