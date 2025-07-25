"""
`GeometricEquation{invType,parType,perType}` is the abstract type all equation types are derived from.

All equations should have fields for defining invariants, parameters and periodicity of the main state variable.
The types of these fields are stored in the following type parameters:

* `invType <: OptionalInvariants`: invariants type
* `parType <: OptionalParameters`: parameters type
* `perType <: OptionalPeriodicity`: periodicity type

The `Optional*` types are all unions of the respective `Null*` types and `NamedTuple` or `AbstractArray`, i.e.,

```julia
const OptionalInvariants = Union{NamedTuple, NullInvariants}
const OptionalParameters = Union{NamedTuple, NullParameters}
const OptionalPeriodicity = Union{AbstractArray, NullPeriodicity}
```

The `Null*` types are empty structs, merely used for dispatch and the traits `hasinvariants`, `hasparameters` and `hasperiodicity`.

"""
abstract type GeometricEquation{invType,parType,perType} end

abstract type DifferentialEquation{invType,parType,perType} <: GeometricEquation{invType,parType,perType} end
abstract type DifferentialAlgebraicEquation{invType,parType,perType,secType} <: GeometricEquation{invType,parType,perType} end
abstract type StochasticDifferentialEquation{invType,parType,perType} <: GeometricEquation{invType,parType,perType} end
abstract type DiscreteEquation{invType,parType,perType} <: GeometricEquation{invType,parType,perType} end

abstract type AbstractEquationODE{invType,parType,perType} <: GeometricEquation{invType,parType,perType} end
abstract type AbstractEquationDAE{invType,parType,perType,secType} <: DifferentialAlgebraicEquation{invType,parType,perType,secType} end
abstract type AbstractEquationSDE{invType,parType,perType} <: StochasticDifferentialEquation{invType,parType,perType} end
abstract type AbstractEquationPODE{invType,parType,perType} <: GeometricEquation{invType,parType,perType} end
abstract type AbstractEquationPDAE{invType,parType,perType,secType} <: DifferentialAlgebraicEquation{invType,parType,perType,secType} end
abstract type AbstractEquationPSDE{invType,parType,perType} <: StochasticDifferentialEquation{invType,parType,perType} end

const GEinvType{invType,parType,perType} = GeometricEquation{invType,parType,perType} # type alias for dispatch on invariants type parameter
const GEparType{parType,invType,perType} = GeometricEquation{invType,parType,perType} # type alias for dispatch on parameters type parameter
const GEperType{perType,invType,parType} = GeometricEquation{invType,parType,perType} # type alias for dispatch on periodicity type parameter

hasinvariants(::GEinvType{<:NullInvariants}) = false
hasinvariants(::GEinvType{<:NamedTuple}) = true

hasparameters(::GEparType{<:NullParameters}) = false
hasparameters(::GEparType{<:NamedTuple}) = true

hasperiodicity(::GEperType{<:NullPeriodicity}) = false
hasperiodicity(::GEperType{<:Tuple{AT,AT}}) where {AT <: AbstractArray} = true

getperiodicity(equ::GEperType{<:NullPeriodicity}) = missing
getperiodicity(equ::GEperType{<:Tuple{AT,AT}}) where {DT, AT <: AbstractArray{DT}} = BitArray(periodicity(equ)[begin][i] ≠ -DT(Inf) && periodicity(equ)[end][i] ≠ +DT(Inf) for i in eachindex(periodicity(equ)[begin], periodicity(equ)[end]))

const DAEsecType{secType,invType,parType,perType} = DifferentialAlgebraicEquation{invType,parType,perType,secType} # type alias for dispatch on secondary constraint type parameter

hassecondary(::DAEsecType{<:Nothing}) = false
hassecondary(::DAEsecType{<:Callable}) = true

_functions(equ::GeometricEquation) = error("_functions(::GeometricEquation) not implemented for ", typeof(equ), ".")
_solutions(equ::GeometricEquation) = error("_solutions(::GeometricEquation) not implemented for ", typeof(equ), ".")
_initialguess(equ::GeometricEquation) = error("_initialguess(::GeometricEquation) not implemented for ", typeof(equ), ".")
_invariants(equ::GeometricEquation) = error("_invariants(::GeometricEquation) not implemented for ", typeof(equ), ".")

_functions(equ::GeometricEquation, ::OptionalParameters) = error("_functions(::GeometricEquation, ::OptionalParameters) not implemented for ", typeof(equ), ".")
_solutions(equ::GeometricEquation, ::OptionalParameters) = error("_solutions(::GeometricEquation, ::OptionalParameters) not implemented for ", typeof(equ), ".")
_initialguess(equ::GeometricEquation, ::OptionalParameters) = error("_initialguess(::GeometricEquation, ::OptionalParameters) not implemented for ", typeof(equ), ".")
_invariants(equ::GeometricEquation, ::OptionalParameters) = error("_invariants(::GeometricEquation, ::OptionalParameters) not implemented for ", typeof(equ), ".")

function GeometricBase.functions(equ::GeometricEquation)
    if hasvectorfield(equ)
        return _functions(equ)
    else
        return NamedTuple()
    end
end

function GeometricBase.functions(equ::GeometricEquation, params::OptionalParameters)
    @assert check_parameters(equ, params)
    if hasvectorfield(equ)
        return _functions(equ, params)
    else
        return NamedTuple()
    end
end

function GeometricBase.solutions(equ::GeometricEquation)
    if hassolution(equ)
        return _solutions(equ)
    else
        return NamedTuple()
    end
end

function GeometricBase.solutions(equ::GeometricEquation, params::OptionalParameters)
    @assert check_parameters(equ, params)
    if hassolution(equ)
        return _solutions(equ, params)
    else
        return NamedTuple()
    end
end

function GeometricBase.initialguess(equ::GeometricEquation)
    if hasinitialguess(equ)
        return _initialguess(equ)
    else
        return NamedTuple()
    end
end

function GeometricBase.initialguess(equ::GeometricEquation, params::OptionalParameters)
    @assert check_parameters(equ, params)
    if hasinitialguess(equ)
        return _initialguess(equ, params)
    else
        return NamedTuple()
    end
end

function GeometricBase.invariants(equ::GeometricEquation, params::OptionalParameters)
    @assert check_parameters(equ, params)
    if hasinvariants(equ)
        return NamedTuple{keys(invariants(equ))}( (_get_invariant(equ, inv, params) for inv in invariants(equ)) )
    else
        return NamedTuple()
    end
end

GeometricBase.parameters(equ::GeometricEquation, args...) = error("parameters(::GeometricEquation) not implemented for ", typeof(equ), ".")
GeometricBase.periodicity(equ::GeometricEquation, args...) = error("periodicity(::GeometricEquation) not implemented for ", typeof(equ), ".")

hassolution(::GeometricEquation) = false
hasvectorfield(::GeometricEquation) = false
hasinitialguess(::GeometricEquation) = false
hasprimary(::GeometricEquation) = false
hassecondary(::GeometricEquation) = false

hasinvariants(::GeometricEquation) = false
hasparameters(::GeometricEquation) = false
hasperiodicity(::GeometricEquation) = false

hashamiltonian(::GeometricEquation) = false
haslagrangian(::GeometricEquation) = false

GeometricBase.datatype(equ::GeometricEquation, ics::NamedTuple) = error("datatype(::GeometricEquation, ::NamedTuple) not implemented for ", typeof(equ), ".")
GeometricBase.arrtype(equ::GeometricEquation, ics::NamedTuple) = error("arrtype(::GeometricEquation, ::NamedTuple) not implemented for ", typeof(equ), ".")

initialstate(::GeometricEquation, ics::NamedTuple) = ics
initialstate(::GeometricEquation, ::InitialTime, ics::NamedTuple, ::OptionalParameters) = ics

check_initial_conditions(equ::GeometricEquation, ics::NamedTuple) = error("check_initial_conditions(::GeometricEquation, ::NamedTuple) not implemented for ", typeof(equ), ".")
check_methods(equ::GeometricEquation, timespan, ics, params) = error("check_methods(::GeometricEquation, ::Tuple, ::NamedTuple, ::OptionalParameters) not implemented for ", typeof(equ), ".")

function check_parameters(equ::GeometricEquation, params::NamedTuple)
    typeof(parameters(equ)) <: NamedTuple || return false
    keys(parameters(equ)) == keys(params) || return false
    for key in keys(params)
        haskey(parameters(equ), key) || return false
        typeof(params[key]) <: parameters(equ)[key] || return false
    end
    return true
end

function check_parameters(equ::GeometricEquation, ::NullParameters)
    if parameters(equ) == NullParameters()
        return true
    else
        return false
    end
end

function check_parameters(::GeometricEquation, ::Any)
    error("Parameters are neither a NamedTuple nor NullParameters.")
end
