@doc raw"""
EquationProblem: stores a GeometricEquation together with initial conditions, parameters, time span and time step size.

### Parameters

* `ST <: GeometricEquation`: super type, used for dispatch
* `DT <: Number`: data type
* `TT <: Real`: time step type
* `AT <: AbstractArray{DT}`: array type of state variable
* `equType <: GeometricEquation`: equation type
* `functionsType <: NamedTuple`: types of all function methods
* `solutionsType <: NamedTuple`: types of all solution methods
* `icsType <: NamedTuple`: types of all initial conditions
* `parType <: OptionalParameters`: parameters type

### Fields

* `equation`: reference to the parent equation object holding the vector fields, etc.
* `functions`: methods for all vector fields, etc., that define the problem
* `solutions`: methods for all solutions, etc., if defined
* `timespan`: time span for problem `(t₀,t₁)`
* `timestep`: time step to be used in simulation
* `ics`: `NamedTuple` containing the initial conditions, must contain one field for each state variable
* `parameters`: either a `NamedTuple` containing the equation's parameters or `NullParameters` indicating that the equation does not have any parameters

### Subtypes

The `EquationProblem` type has various subtypes for the different equations types, that are defined e.g. via
```julia
const ODEProblem = EquationProblem{ODE}
```
and provide convenience constructors to construct an equation and the corresponding problem in one step, e.g.,
```julia
ODEProblem(v, timespan, timestep, ics::NamedTuple; kwargs...)
ODEProblem(v, timespan, timestep, q₀::StateVariable; kwargs...)
```

All problem subtypes take the following keyword arguments:

* `invariants = NullInvariants()`
* `parameters = NullParameters()`
* `periodicity = NullPeriodicity()`

If not set to their corresponding Null types, the user needs to pass a `NamedTuple` whose values are

* functions for invariants,
* arbitrary data structures for parameters,
* the same data structure as the solution for periodicity.

The latter should be zero everywhere, except for those components, that are periodic, i.e.,
whose value are supposed to stay within a range `(0, max)`. Support for ranges starting
with other values than zero is currently missing but can be added if demand arises.

"""
struct EquationProblem{superType <: GeometricEquation, dType <: Number, tType <: Real,
    arrayType <: AbstractArray{dType},
    equType <: GeometricEquation,
    functionsType <: NamedTuple,
    solutionsType <: NamedTuple,
    iguessType <: NamedTuple,
    icsType <: NamedTuple,
    paramsType <: OptionalParameters} <:
       GeometricProblem{superType, dType, tType}
    equation::equType
    functions::functionsType
    solutions::solutionsType
    initialguess::iguessType
    timespan::Tuple{tType, tType}
    timestep::tType
    ics::icsType
    parameters::paramsType

    function EquationProblem(equ, timespan, timestep, ics, parameters)
        _timespan = promote_timespan(timespan)
        _timespan, _timestep = promote_timespan_and_timestep(_timespan, timestep)
        _ics = initialstate(equ, _timespan[begin], ics, parameters)

        @assert check_initial_conditions(equ, _ics)
        @assert check_methods(equ, _timespan, _ics, parameters)
        check_noise(equ, _timespan, _timestep)

        superType = eval(typeof(equ).name.name)
        tType = typeof(_timestep)
        dType = datatype(equ, _ics)
        arrayType = arrtype(equ, _ics)

        funcs = functions(equ)
        sols = solutions(equ)
        iguess = initialguess(equ)

        new{superType, dType, tType, arrayType, typeof(equ), typeof(funcs),
            typeof(sols), typeof(iguess), typeof(_ics), typeof(parameters)
        }(equ, funcs, sols, iguess, _timespan, _timestep, _ics, parameters)
    end
end

function EquationProblem(equ, timespan, timestep, ics, ::Nothing)
    EquationProblem(equ, timespan, timestep, ics, NullParameters())
end

function EquationProblem(equ, timespan, timestep, ics; parameters = NullParameters())
    EquationProblem(equ, timespan, timestep, ics, parameters)
end

# Base.hash(prob::EquationProblem, h::UInt) = hash(hash(prob.equation,
#                                                   hash(prob.timespan, hash(prob.timestep,
#                                                   hash(prov.ics, hash(prov.parameters, h))))))

function Base.:(==)(prob1::EquationProblem, prob2::EquationProblem)
    (
        prob1.equation == prob2.equation
        && prob1.functions == prob2.functions
        && prob1.solutions == prob2.solutions
        && prob1.initialguess == prob2.initialguess
        && prob1.timespan == prob2.timespan
        && prob1.timestep == prob2.timestep
        && prob1.ics == prob2.ics
        && prob1.parameters == prob2.parameters)
end

@inline GeometricBase.datatype(::EquationProblem{
    ST, DT, TT, AT}) where {ST, DT, TT, AT} = DT
@inline GeometricBase.timetype(::EquationProblem{
    ST, DT, TT, AT}) where {ST, DT, TT, AT} = TT
@inline GeometricBase.arrtype(::EquationProblem{ST, DT, TT, AT}) where {ST, DT, TT, AT} = AT
@inline GeometricBase.equtype(::EquationProblem{ST, DT, TT, AT}) where {ST, DT, TT, AT} = ST

@inline GeometricBase.equation(prob::EquationProblem) = prob.equation
@inline GeometricBase.timespan(prob::EquationProblem) = prob.timespan
@inline GeometricBase.timestep(prob::EquationProblem) = prob.timestep

@inline GeometricBase.functions(prob::EquationProblem) = prob.functions
@inline GeometricBase.solutions(prob::EquationProblem) = prob.solutions
@inline GeometricBase.initialguess(prob::EquationProblem) = prob.initialguess
@inline GeometricBase.parameters(prob::EquationProblem) = prob.parameters
@inline GeometricBase.nsamples(::EquationProblem) = 1

# A fixed-width integer, or a rational of them, in the next wider type; any other value as it is.
_widen(x::Union{Base.BitInteger, Rational{<:Base.BitInteger}}) = widen(x)
_widen(x) = x

"""
    ntimesteps(timespan, timestep)

The number of whole steps `timestep` in `timespan = (t₀, t₁)`, with `t₀`, `t₁` and `timestep`
promoted to one type `T`. Floating-point times are converted to `W = promote_type(T, Float64)`,
so that `t₁ - t₀` does not round in a narrower type, and every expression below is evaluated in
`W`: a count `n` whose end `t₀ + n * timestep` lies within `8 * eps(T) * max(abs(t₀), abs(t₁))` of
`t₁` fills the span, on either side of `t₁`, so a quotient `(t₁ - t₀) / timestep` that rounds just
below or above an integer neither loses nor adds a step. Any other count rounds toward zero, and
the last step ends before `t₁`. Exact time types, such as integers and rationals, count the whole
steps exactly; fixed-width integers and their rationals count in `widen(T)`, so that `t₁ - t₀`
does not overflow. The count is the magnitude, whatever the sign of `timestep`. A zero `timestep`
throws an `ArgumentError`.

For a time step that points from `t₀` to `t₁`, this is the rule of `TimeSeries(t₀, t₁, timestep)`
in GeometricSolutions, so the step count of a problem is the number of steps of its solution's
time grid.
"""
@inline function ntimesteps(timespan, timestep)
    iszero(timestep) && throw(ArgumentError("time step cannot be zero"))
    t₀, t₁, h = promote(timespan[begin], timespan[end], timestep)
    h isa AbstractFloat ||
        return abs(Int(div(_widen(t₁) - _widen(t₀), _widen(h), RoundToZero)))
    W = promote_type(typeof(h), Float64)
    x₀, x₁, y = W(t₀), W(t₁), W(h)
    q = (x₁ - x₀) / y
    n = round(q)
    fills = abs(x₀ + n * y - x₁) ≤ 8 * eps(typeof(h)) * max(abs(x₀), abs(x₁))
    return abs(Int(fills ? n : trunc(q)))
end

@inline function GeometricBase.ntime(prob::EquationProblem)
    ntimesteps(timespan(prob), timestep(prob))
end

@inline function GeometricBase.invariants(prob::EquationProblem)
    invariants(equation(prob))
end

function initial_conditions(prob::EquationProblem)
    merge((t = TimeVariable(initialtime(prob)),), prob.ics)
end

function initialstate(prob::EquationProblem)
    State(initialtime(prob), prob.ics)
end

function Base.show(io::IO, prob::EquationProblem)
    print(io, "Geometric Equation Problem for ", equation(prob))
    print(io, "\n\n")
    print(io, " Timespan: $(timespan(prob)) \n")
    print(io, " Timestep: $(timestep(prob)) \n")
    print(io, "\n")
    print(io, " Initial conditions: \n")
    print(io, "   ", initial_conditions(prob))
    print(io, "\n\n")
    print(io, " Parameters: \n")
    print(io, "   ", parameters(prob))
end

function Base.similar(
        prob::EquationProblem, timespan, timestep = timestep(prob), ics = prob.ics,
        parameters = parameters(prob))
    EquationProblem(equation(prob), timespan, timestep,
        initialstate(equation(prob), ics...), parameters)
end

function Base.similar(
        prob::EquationProblem; timespan = timespan(prob), timestep = timestep(prob),
        ics = prob.ics, parameters = parameters(prob))
    similar(prob, timespan, timestep, initialstate(equation(prob), ics...), parameters)
end
