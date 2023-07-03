"""
GeometricEnsemble: stores a GeometricEquation together with multiple sets of initial conditions, parameters, time span and time step size.

"""
struct GeometricEnsemble{superType <: GeometricEquation, dType <: Number, tType <: Real,
                 arrayType <: AbstractArray{dType}, 
                 equType <: GeometricEquation,
                 functionsType <: NamedTuple,
                 solutionsType <: NamedTuple,
                 icsType <: AbstractVector{<:NamedTuple},
                 paramsType <: AbstractVector{<:OptionalParameters}} <: AbstractProblem
    equation::equType
    functions::functionsType
    solutions::solutionsType
    tspan::Tuple{tType,tType}
    tstep::tType
    ics::icsType
    parameters::paramsType
end

function GeometricEnsemble(equ::equType, tspan, tstep, ics::AbstractVector{<:NamedTuple}, parameters::AbstractVector{<:OptionalParameters}) where {equType}
    _tspan = promote_tspan(tspan)
    _tspan, _tstep = promote_tspan_and_tstep(_tspan, tstep)

    for ic in ics
        @assert check_initial_conditions(equ, ic)
        @assert typeof(ic) == typeof(ics[begin])
    end

    @assert check_methods(equ, _tspan, ics[begin], parameters)
    @assert axes(parameters) == axes(ics)

    superType = eval(typeof(equ).name.name)
    tType = typeof(_tstep)
    dType = datatype(equ, ics[begin])
    arrayType = arrtype(equ, ics[begin])

    funcs = functions(equ)
    sols = solutions(equ)

    GeometricEnsemble{superType, dType, tType, arrayType, equType, typeof(funcs), typeof(sols), typeof(ics), typeof(parameters)}(equ, funcs, sols, _tspan, _tstep, ics, parameters)
end

function GeometricEnsemble(equ, tspan, tstep, ics::AbstractVector{<:NamedTuple}, params::OptionalParameters=NullParameters())
    _params = similar(ics, typeof(params))

    for i in eachindex(_params)
        _params[i] = params
    end

    GeometricEnsemble(equ, tspan, tstep, ics, _params)
end

function GeometricEnsemble(equ, tspan, tstep, ics::NamedTuple, params::AbstractVector{<:OptionalParameters})
    _ics = similar(params, typeof(ics))

    for i in eachindex(_ics)
        _ics[i] = ics
    end

    GeometricEnsemble(equ, tspan, tstep, _ics, params)
end

function GeometricEnsemble(equ, tspan, tstep, ics, ::Nothing)
    GeometricEnsemble(equ, tspan, tstep, ics, NullParameters())
end

function GeometricEnsemble(equ, tspan, tstep, ics; parameters = NullParameters())
    GeometricEnsemble(equ, tspan, tstep, ics, parameters)
end


@inline GeometricBase.datatype(::GeometricEnsemble{ST, DT, TT, AT}) where {ST, DT, TT, AT} = DT
@inline GeometricBase.timetype(::GeometricEnsemble{ST, DT, TT, AT}) where {ST, DT, TT, AT} = TT
@inline GeometricBase.arrtype(::GeometricEnsemble{ST, DT, TT, AT}) where {ST, DT, TT, AT} = AT
@inline GeometricBase.equtype(::GeometricEnsemble{ST, DT, TT, AT}) where {ST, DT, TT, AT} = ST

@inline GeometricBase.equation(ge::GeometricEnsemble) = ge.equation
@inline GeometricBase.tspan(ge::GeometricEnsemble) = ge.tspan
@inline GeometricBase.tstep(ge::GeometricEnsemble) = ge.tstep

@inline GeometricBase.timestep(ge::GeometricEnsemble) = tstep(ge)
@inline GeometricBase.functions(ge::GeometricEnsemble) = ge.functions
@inline GeometricBase.solutions(ge::GeometricEnsemble) = ge.solutions
@inline GeometricBase.parameters(ge::GeometricEnsemble) = ge.parameters

initial_conditions(ge::GeometricEnsemble) = ge.ics

@inline GeometricBase.nsamples(ge::GeometricEnsemble) = length(initial_conditions(ge))

initial_condition(ge::GeometricEnsemble, i) = initial_conditions(ge)[i]
parameter(ge::GeometricEnsemble, i) = parameters(ge)[i]

function problem(ge::GeometricEnsemble, i)
    GeometricProblem(equation(ge), tspan(ge), tstep(ge), initial_condition(ge, i), parameter(ge, i))
end

Base.iterate(ge::GeometricEnsemble, i=1) = i > nsamples(ge) ? nothing : (problem(ge, i), i+1)
