"""

"""
struct GeometricEnsemble{superType <: GeometricEquation, dType <: Number, tType <: Real, arrayType <: AbstractArray{dType}, 
                 equType <: GeometricEquation, fType <: NamedTuple, icsType <: AbstractVector{<:NamedTuple}, paramsType <: AbstractVector{<:OptionalParameters}} <: AbstractProblem{dType, tType, arrayType}
    equation::equType
    functions::fType
    tspan::Tuple{tType,tType}
    tstep::tType
    ics::icsType
    parameters::paramsType
end

function GeometricEnsemble(equ::equType, tspan, tstep, ics::AbstractVector{<:NamedTuple}, parameters::AbstractVector{<:OptionalParameters}) where {equType}
    for ic in ics
        @assert check_initial_conditions(equ, ic)
        @assert typeof(ic) == typeof(ics[begin])
        @assert axes(ic) == axes(ics[begin])
    end

    @assert axes(parameters) == axes(ics)

    superType = eval(typeof(equ).name.name)
    funcs = get_functions(equ)

    _tspan = promote_tspan(tspan)
    _tspan, _tstep = promote_tspan_and_tstep(_tspan, tstep)

    tType = typeof(_tstep)
    dType = eltype(ics[begin])
    arrayType = typeof(ics[begin])

    GeometricEnsemble{superType, dType, tType, arrayType, equType, typeof(funcs), eltype(ics), eltype(parameters)}(equ, funcs, _tspan, _tstep, ics, parameters)
end

function GeometricEnsemble(equ, tspan, tstep, ics::AbstractVector{<:NamedTuple}, params::OptionalParameters=NullParameters())
    _params = similar(ics, typeof(parameters))

    for i in eachindex(_params)
        _params[i] = parameters
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


Base.eltype(::GeometricEnsemble{ST,DT,TT}) where {ST,DT,TT} = DT
timetype(::GeometricEnsemble{ST,DT,TT}) where {ST,DT,TT} = TT

GeometricBase.nsamples(ge::GeometricEnsemble) = length(ge.ics)

equation(ge::GeometricEnsemble) = ge.equation
functions(ge::GeometricEnsemble) = ge.functions
parameters(ge::GeometricEnsemble) = ge.parameters
tspan(ge::GeometricEnsemble) = ge.tspan
tstep(ge::GeometricEnsemble) = ge.tstep


const ODEEnsemble   = GeometricEnsemble{ODE}
const IODEEnsemble  = GeometricEnsemble{IODE}
const PODEEnsemble  = GeometricEnsemble{PODE}
const HODEEnsemble  = GeometricEnsemble{HODE}
const LODEEnsemble  = GeometricEnsemble{LODE}
const SODEEnsemble  = GeometricEnsemble{SODE}
const DAEEnsemble   = GeometricEnsemble{DAE}
const IDAEEnsemble  = GeometricEnsemble{IDAE}
const PDAEEnsemble  = GeometricEnsemble{PDAE}
const HDAEEnsemble  = GeometricEnsemble{HDAE}
const LDAEEnsemble  = GeometricEnsemble{LDAE}
const SPDAEEnsemble = GeometricEnsemble{SPDAE}
const SDEEnsemble   = GeometricEnsemble{SDE}
const PSDEEnsemble  = GeometricEnsemble{PSDE}
const SPSDEEnsemble = GeometricEnsemble{SPSDE}
