@doc raw"""
GeometricProblem: stores a GeometricEquation together with initial conditions, parameters, time span and time step size.

### Parameters

* `ST <: GeometricEquation`: super type, used for dispatch
* `DT <: Number`: data type
* `TT <: Real`: time step type
* `AT <: AbstractArray{DT}`: array type of state variable
* `equType <: GeometricEquation`: equation type
* `icsType <: paramsType`: initial conditions type
* `parType <: OptionalParameters`: parameters type

### Fields

* `equation`: reference to the parent equation object holding the vector fields, etc.
* `tspan`: time span for problem `(t₀,t₁)`
* `tstep`: time step to be used in simulation
* `ics`: `NamedTuple` containing the initial conditions, must contain one field for each state variable
* `parameters`: either a `NamedTuple` containing the equation's parameters or `NullParameters` indicating that the equation does not have any parameters

### Subtypes

The `GeometricProblem` type has various subtypes for the different equations types, that are defined e.g. via
```julia
const ODEProblem = GeometricProblem{ODE}
```
and provide convenience constructors to construct an equation and the corresponding problem in one step, e.g.,
```julia
ODEProblem(v, tspan, tstep, ics::NamedTuple; kwargs...)
ODEProblem(v, tspan, tstep, q₀::State; kwargs...)
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
struct GeometricProblem{superType <: GeometricEquation, dType <: Number, tType <: Real,
                        arrayType <: AbstractArray{dType},
                        equType <: GeometricEquation,
                        functionsType <: NamedTuple,
                        solutionsType <: NamedTuple,
                        icsType <: NamedTuple,
                        paramsType <: OptionalParameters} <: AbstractProblem
    equation::equType
    functions::functionsType
    solutions::solutionsType
    tspan::Tuple{tType, tType}
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

        funcs = functions(equ, parameters)
        sols = solutions(equ, parameters)

        new{superType, dType, tType, arrayType, typeof(equ), typeof(funcs), typeof(sols), typeof(ics), typeof(parameters)
            }(equ, funcs, sols, _tspan, _tstep, ics, parameters)
    end
end

function GeometricProblem(equ, tspan, tstep, ics, ::Nothing)
    GeometricProblem(equ, tspan, tstep, ics, NullParameters())
end

function GeometricProblem(equ, tspan, tstep, ics; parameters = NullParameters())
    GeometricProblem(equ, tspan, tstep, ics, parameters)
end

# Base.hash(prob::GeometricProblem, h::UInt) = hash(hash(prob.equation,
#                                                   hash(prob.tspan, hash(prob.tstep,
#                                                   hash(prov.ics, hash(prov.parameters, h))))))

# Base.:(==)(prob1::GeometricProblem, prob2::GeometricProblem) = (
#                                 prob1.equation   == prob2.equation
#                              && prob1.tspan      == prob2.tspan
#                              && prob1.tstep      == prob2.tstep
#                              && prob1.ics        == prob2.ics
#                              && prob1.parameters == prob2.parameters)

@inline GeometricBase.datatype(::GeometricProblem{ST, DT, TT, AT}) where {ST, DT, TT, AT} = DT
@inline GeometricBase.timetype(::GeometricProblem{ST, DT, TT, AT}) where {ST, DT, TT, AT} = TT
@inline GeometricBase.arrtype(::GeometricProblem{ST, DT, TT, AT}) where {ST, DT, TT, AT} = AT
@inline GeometricBase.equtype(::GeometricProblem{ST, DT, TT, AT}) where {ST, DT, TT, AT} = ST

@inline GeometricBase.equation(prob::GeometricProblem) = prob.equation
@inline GeometricBase.tspan(prob::GeometricProblem) = prob.tspan
@inline GeometricBase.tstep(prob::GeometricProblem) = prob.tstep

@inline GeometricBase.timestep(prob::GeometricProblem) = tstep(prob)
@inline GeometricBase.functions(prob::GeometricProblem) = prob.functions
@inline GeometricBase.solutions(prob::GeometricProblem) = prob.solutions
@inline GeometricBase.parameters(prob::GeometricProblem) = prob.parameters
@inline GeometricBase.nsamples(::GeometricProblem) = 1

@inline function GeometricBase.ntime(prob::GeometricProblem)
    Int(abs(div(tend(prob) - tbegin(prob), tstep(prob), RoundUp)))
end

@inline function GeometricBase.invariants(prob::GeometricProblem)
    invariants(equation(prob), parameters(prob))
end

initial_conditions(prob::GeometricProblem) = merge((t = tbegin(prob),), prob.ics)

function Base.similar(prob::GeometricProblem, tspan, tstep = tstep(prob), ics = prob.ics,
                      parameters = parameters(prob))
    GeometricProblem(equation(prob), tspan, tstep, ics, parameters)
end

function Base.similar(prob::GeometricProblem; tspan = tspan(prob), tstep = tstep(prob),
                      ics = prob.ics, parameters = parameters(prob))
    similar(prob, tspan, tstep, ics, parameters)
end
