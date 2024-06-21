"""
EnsembleProblem: stores a GeometricEquation together with multiple sets of initial conditions and/or parameters, a time span for integration and a time step size.

An `EnsembleProblem` is initialized by providing a `GeometricEquation`, an integration time span (typically a tuple with two values, start and end time, respectively), a timestep, and one of the three following options:
* a vector of initial conditions and a vector of parameter sets, with both vectors having the same length,
* a vector of initial conditions and a single set of parameters,
* a single initial condition and a vector of parameter sets.

Each initial condition is a `NamedTuple` that contains one field for each state variable.
Each parameter set is either a `NamedTuple` containing the equation's parameters or `NullParameters` indicating that the equation does not have any parameters.

The different constructors then generate two vectors, one for the initial conditions and one for the parameters, where each pair of the corresponding entries defines one problem.
In the first case, the respective constructor just checks if both vectors are of the same size.
In the second case, the respective constructor creates a parameter vector, where each entry holds the same parameter set.
In the third case, the respective constructor creates an initial condition vector, where each entry holds the same initial conditions.
One may be inclined to think that the first and second constructor lead to a waste of memory, but in reality the respective vectors only hold references to the same initial conditons or parameters.
Thus the data is not actually duplicated.

Each pair of initial conditions and parameters is referred to as sample.
The methods
```
length(::EnsembleProblem)
nsamples(::EnsembleProblem)
```
return the number of samples in a `EnsembleProblem`.
A single initial condition or parameter set can be retrieved by the methods
```
initial_condition(::EnsembleProblem, i)
parameter(::EnsembleProblem, i)
```
where `i` is the index of the sample.
Typically, however, e.g. when integrating all samples of an `EnsembleProblem` with GeometricIntegrators, it is more convenient to retrieve the corresponding `GeometricProblem` via the method
```
problem(::EnsembleProblem, i)
```
The `EnsembleProblem` also allows to iterate over all samples, e.g.
```
for problem in ensemble
    # ...
    # integrate problem
    # ...
end
```
where `ensemble` is an `EnsembleProblem` and `problem` is the corresponding `GeometricProblem`.


### Parameters

* `ST <: GeometricEquation`: super type, used for dispatch
* `DT <: Number`: data type
* `TT <: Real`: time step type
* `AT <: AbstractArray{DT}`: array type of state variable
* `equType <: GeometricEquation`: equation type
* `functionsType <: NamedTuple`: types of all function methods
* `solutionsType <: NamedTuple`: types of all solution methods
* `icsType <: AbstractVector{<:NamedTuple}`: types of all initial conditions 
* `parType <: AbstractVector{<:OptionalParameters}`: parameters type

### Fields

* `equation`: reference to the parent equation object holding the vector fields, etc.
* `functions`: methods for all vector fields, etc., that define the problem
* `solutions`: methods for all solutions, etc., if defined
* `tspan`: time span for problem `(t₀,t₁)`
* `tstep`: time step to be used in simulation
* `ics`: vector of `NamedTuple` containing the initial conditions, each `NamedTuple` must contain one field for each state variable
* `parameters`: vector of either `NamedTuple` containing the equation's parameters or `NullParameters` indicating that the equation does not have any parameters

### Constructors

The `EnsembleProblem` provides the following constructors:
```
EnsembleProblem(equ, tspan, tstep, ics::AbstractVector{<:NamedTuple}, parameters::AbstractVector{<:OptionalParameters})
EnsembleProblem(equ, tspan, tstep, ics::AbstractVector{<:NamedTuple}, parameters::OptionalParameters=NullParameters())
EnsembleProblem(equ, tspan, tstep, ics::NamedTuple, parameters::AbstractVector{<:OptionalParameters})
EnsembleProblem(equ, tspan, tstep, ics, ::Nothing) = 
    EnsembleProblem(equ, tspan, tstep, ics, NullParameters())
EnsembleProblem(equ, tspan, tstep, ics; parameters = NullParameters()) = 
    EnsembleProblem(equ, tspan, tstep, ics, parameters)
```

* `equ` is a subtype of `GeometricEquation`
* `tspan` is a tuple `(t₀,t₁)` of the integration time span with `t₀` the start time and `t₁` the end time
* `tstep` is the time step, typically a value of some `AbstractFloat` subtype
* `ics` are the initial conditions, either a single set or a vector of multiple sets
* `parameters` are the static parameters of the problem, either a single set or a vector of multiple sets

"""
struct EnsembleProblem{superType <: GeometricEquation, dType <: Number, tType <: Real,
                 arrayType <: AbstractArray{dType}, 
                 equType <: GeometricEquation,
                 functionsType <: NamedTuple,
                 solutionsType <: NamedTuple,
                 icsType <: AbstractVector{<:NamedTuple},
                 paramsType <: AbstractVector{<:OptionalParameters}} <: GeometricProblem{superType}
    equation::equType
    functions::functionsType
    solutions::solutionsType
    tspan::Tuple{tType,tType}
    tstep::tType
    ics::icsType
    parameters::paramsType
end

function EnsembleProblem(equ::equType, tspan, tstep, ics::AbstractVector{<:NamedTuple}, parameters::AbstractVector{<:OptionalParameters}) where {equType}
    @assert axes(ics) == axes(parameters)

    _tspan = promote_tspan(tspan)
    _tspan, _tstep = promote_tspan_and_tstep(_tspan, tstep)
    _ics = [initialstate(equ, _tspan[begin], ic, param) for (ic,param) in zip(ics,parameters)]

    for i in eachindex(_ics)
        @assert check_initial_conditions(equ, _ics[i])
    end

    @assert check_methods(equ, _tspan, _ics[begin], parameters[begin])
    @assert axes(parameters) == axes(ics)

    length(_ics) == 1 && @warn("You created an EnsembleProblem with a single initial condition and a single set of parameters. You probably want to create a GeometricProblem instead.")

    superType = eval(typeof(equ).name.name)
    tType = typeof(_tstep)
    dType = datatype(equ, _ics[begin])
    arrayType = arrtype(equ, _ics[begin])

    funcs = functions(equ)
    sols = solutions(equ)

    EnsembleProblem{superType, dType, tType, arrayType, equType, typeof(funcs), typeof(sols), typeof(_ics), typeof(parameters)}(equ, funcs, sols, _tspan, _tstep, _ics, parameters)
end

function EnsembleProblem(equ, tspan, tstep, ics::AbstractVector{<:NamedTuple}, parameters::OptionalParameters = NullParameters())
    _params = similar(ics, typeof(parameters))

    for i in eachindex(_params)
        _params[i] = parameters
    end

    EnsembleProblem(equ, tspan, tstep, ics, _params)
end

function EnsembleProblem(equ, tspan, tstep, ics::NamedTuple, parameters::AbstractVector{<:OptionalParameters})
    _ics = similar(parameters, typeof(ics))

    for i in eachindex(_ics)
        _ics[i] = ics
    end

    EnsembleProblem(equ, tspan, tstep, _ics, parameters)
end

function EnsembleProblem(equ, tspan, tstep, ics::NamedTuple, parameters::OptionalParameters)
    EnsembleProblem(equ, tspan, tstep, ics, [parameters])
end

function EnsembleProblem(equ, tspan, tstep, ics, ::Nothing)
    EnsembleProblem(equ, tspan, tstep, ics, NullParameters())
end

function EnsembleProblem(equ, tspan, tstep, ics; parameters = NullParameters())
    EnsembleProblem(equ, tspan, tstep, ics, parameters)
end

Base.:(==)(ens1::EnsembleProblem, ens2::EnsembleProblem) = (
                                ens1.equation   == ens2.equation
                             && ens1.functions  == ens2.functions
                             && ens1.solutions  == ens2.solutions
                             && ens1.tspan      == ens2.tspan
                             && ens1.tstep      == ens2.tstep
                             && ens1.ics        == ens2.ics
                             && ens1.parameters == ens2.parameters)

@inline GeometricBase.datatype(::EnsembleProblem{ST, DT, TT, AT}) where {ST, DT, TT, AT} = DT
@inline GeometricBase.timetype(::EnsembleProblem{ST, DT, TT, AT}) where {ST, DT, TT, AT} = TT
@inline GeometricBase.arrtype(::EnsembleProblem{ST, DT, TT, AT}) where {ST, DT, TT, AT} = AT
@inline GeometricBase.equtype(::EnsembleProblem{ST, DT, TT, AT}) where {ST, DT, TT, AT} = ST

@inline GeometricBase.equation(ge::EnsembleProblem) = ge.equation
@inline GeometricBase.tspan(ge::EnsembleProblem) = ge.tspan
@inline GeometricBase.tstep(ge::EnsembleProblem) = ge.tstep

@inline GeometricBase.timestep(ge::EnsembleProblem) = tstep(ge)
@inline GeometricBase.functions(ge::EnsembleProblem) = ge.functions
@inline GeometricBase.solutions(ge::EnsembleProblem) = ge.solutions
@inline GeometricBase.parameters(ge::EnsembleProblem) = ge.parameters


@inline GeometricBase.nsamples(ge::EnsembleProblem) = length(initial_conditions(ge))

initial_conditions(ge::EnsembleProblem) = ge.ics
initial_condition(ge::EnsembleProblem, i) = initial_conditions(ge)[i]
parameter(ge::EnsembleProblem, i) = parameters(ge)[i]

function problem(ge::EnsembleProblem, i)
    EquationProblem(equation(ge), tspan(ge), tstep(ge), initial_condition(ge, i), parameter(ge, i))
end

Base.length(ge::EnsembleProblem) = nsamples(ge)
Base.iterate(ge::EnsembleProblem, i=1) = i > nsamples(ge) ? nothing : (problem(ge, i), i+1)
Base.getindex(ge::EnsembleProblem, i::Union{Integer,Base.AbstractCartesianIndex}) = problem(ge, i)

initialstate(::GeometricEquation, ics::AbstractVector{<:NamedTuple}) = ics
