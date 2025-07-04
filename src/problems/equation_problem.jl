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
* `tspan`: time span for problem `(t₀,t₁)`
* `tstep`: time step to be used in simulation
* `ics`: `NamedTuple` containing the initial conditions, must contain one field for each state variable
* `parameters`: either a `NamedTuple` containing the equation's parameters or `NullParameters` indicating that the equation does not have any parameters

### Subtypes

The `EquationProblem` type has various subtypes for the different equations types, that are defined e.g. via
```julia
const ODEProblem = EquationProblem{ODE}
```
and provide convenience constructors to construct an equation and the corresponding problem in one step, e.g.,
```julia
ODEProblem(v, tspan, tstep, ics::NamedTuple; kwargs...)
ODEProblem(v, tspan, tstep, q₀::StateVariable; kwargs...)
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
                        paramsType <: OptionalParameters} <: GeometricProblem{superType}
    equation::equType
    functions::functionsType
    solutions::solutionsType
    initialguess::iguessType
    tspan::Tuple{tType, tType}
    tstep::tType
    ics::icsType
    parameters::paramsType

    function EquationProblem(equ, tspan, tstep, ics, parameters)
        _tspan = promote_tspan(tspan)
        _tspan, _tstep = promote_tspan_and_tstep(_tspan, tstep)
        _ics = initialstate(equ, _tspan[begin], ics, parameters)

        @assert check_initial_conditions(equ, _ics)
        @assert check_methods(equ, _tspan, _ics, parameters)

        superType = eval(typeof(equ).name.name)
        tType = typeof(_tstep)
        dType = datatype(equ, _ics)
        arrayType = arrtype(equ, _ics)

        funcs = functions(equ)
        sols = solutions(equ)
        iguess = initialguess(equ)

        new{superType, dType, tType, arrayType, typeof(equ), typeof(funcs), typeof(sols), typeof(iguess), typeof(_ics), typeof(parameters)
            }(equ, funcs, sols, iguess, _tspan, _tstep, _ics, parameters)
    end
end

function EquationProblem(equ, tspan, tstep, ics, ::Nothing)
    EquationProblem(equ, tspan, tstep, ics, NullParameters())
end

function EquationProblem(equ, tspan, tstep, ics; parameters = NullParameters())
    EquationProblem(equ, tspan, tstep, ics, parameters)
end

# Base.hash(prob::EquationProblem, h::UInt) = hash(hash(prob.equation,
#                                                   hash(prob.tspan, hash(prob.tstep,
#                                                   hash(prov.ics, hash(prov.parameters, h))))))

Base.:(==)(prob1::EquationProblem, prob2::EquationProblem) = (
                                prob1.equation   == prob2.equation
                             && prob1.functions  == prob2.functions
                             && prob1.solutions  == prob2.solutions
                             && prob1.initialguess == prob2.initialguess
                             && prob1.tspan      == prob2.tspan
                             && prob1.tstep      == prob2.tstep
                             && prob1.ics        == prob2.ics
                             && prob1.parameters == prob2.parameters)

@inline GeometricBase.datatype(::EquationProblem{ST, DT, TT, AT}) where {ST, DT, TT, AT} = DT
@inline GeometricBase.timetype(::EquationProblem{ST, DT, TT, AT}) where {ST, DT, TT, AT} = TT
@inline GeometricBase.arrtype(::EquationProblem{ST, DT, TT, AT}) where {ST, DT, TT, AT} = AT
@inline GeometricBase.equtype(::EquationProblem{ST, DT, TT, AT}) where {ST, DT, TT, AT} = ST

@inline GeometricBase.equation(prob::EquationProblem) = prob.equation
@inline GeometricBase.timespan(prob::EquationProblem) = prob.tspan
@inline GeometricBase.timestep(prob::EquationProblem) = prob.tstep

@inline GeometricBase.functions(prob::EquationProblem) = prob.functions
@inline GeometricBase.solutions(prob::EquationProblem) = prob.solutions
@inline GeometricBase.initialguess(prob::EquationProblem) = prob.initialguess
@inline GeometricBase.parameters(prob::EquationProblem) = prob.parameters
@inline GeometricBase.nsamples(::EquationProblem) = 1

@inline function GeometricBase.ntime(prob::EquationProblem)
    Int(abs(div(finaltime(prob) - initialtime(prob), timestep(prob), RoundUp)))
end

@inline function GeometricBase.invariants(prob::EquationProblem)
    invariants(equation(prob))
end

initial_conditions(prob::EquationProblem) = merge((t = initialtime(prob),), prob.ics)

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

function Base.similar(prob::EquationProblem, tspan, tstep = timestep(prob), ics = prob.ics,
                      parameters = parameters(prob))
    EquationProblem(equation(prob), tspan, tstep, initialstate(equation(prob), ics...), parameters)
end

function Base.similar(prob::EquationProblem; tspan = timespan(prob), tstep = timestep(prob),
                      ics = prob.ics, parameters = parameters(prob))
    similar(prob, tspan, tstep, initialstate(equation(prob), ics...), parameters)
end
