
abstract type GeometricEquation end

abstract type AbstractEquationODE  <: GeometricEquation end
abstract type AbstractEquationDAE  <: GeometricEquation end
abstract type AbstractEquationSDE  <: GeometricEquation end
abstract type AbstractEquationPODE <: GeometricEquation end
abstract type AbstractEquationPDAE <: GeometricEquation end
abstract type AbstractEquationPSDE <: GeometricEquation end


_functions(equ::GeometricEquation) = error("_functions(::GeometricEquation) not implemented for ", typeof(equ), ".")
_solutions(equ::GeometricEquation) = error("_solutions(::GeometricEquation) not implemented for ", typeof(equ), ".")
_functions(equ::GeometricEquation, ::OptionalParameters) = error("_functions(::GeometricEquation, ::OptionalParameters) not implemented for ", typeof(equ), ".")
_solutions(equ::GeometricEquation, ::OptionalParameters) = error("_solutions(::GeometricEquation, ::OptionalParameters) not implemented for ", typeof(equ), ".")

function functions(equ::GeometricEquation)
    if hasvectorfield(equ)
        return _functions(equ)
    else
        return NamedTuple()
    end
end

function functions(equ::GeometricEquation, params::OptionalParameters)
    @assert check_parameters(equ, params)
    if hasvectorfield(equ)
        return _functions(equ, params)
    else
        return NamedTuple()
    end
end

function solutions(equ::GeometricEquation)
    if hassolution(equ)
        return _solutions(equ)
    else
        return NamedTuple()
    end
end

function solutions(equ::GeometricEquation, params::OptionalParameters)
    @assert check_parameters(equ, params)
    if hassolution(equ)
        return _solutions(equ, params)
    else
        return NamedTuple()
    end
end

GeometricBase.invariants(equ::GeometricEquation) = error("invariants(::GeometricEquation) not implemented for ", typeof(equ), ".")
GeometricBase.parameters(equ::GeometricEquation) = error("parameters(::GeometricEquation) not implemented for ", typeof(equ), ".")
GeometricBase.periodicity(equ::GeometricEquation) = error("periodicity(::GeometricEquation) not implemented for ", typeof(equ), ".")

hassolution(::GeometricEquation) = false
hasvectorfield(::GeometricEquation) = false
hasprimary(::GeometricEquation) = false
hassecondary(::GeometricEquation) = false

hasinvariants(::GeometricEquation) = false
hasparameters(::GeometricEquation) = false
hasperiodicity(::GeometricEquation) = false

hashamiltonian(::GeometricEquation) = false
haslagrangian(::GeometricEquation) = false

datatype(equ::GeometricEquation, ics::NamedTuple) = error("datatype(::GeometricEquation, ::NamedTuple) not implemented for ", typeof(equ), ".")
arrtype(equ::GeometricEquation, ics::NamedTuple) = error("arrtype(::GeometricEquation, ::NamedTuple) not implemented for ", typeof(equ), ".")

check_initial_conditions(equ::GeometricEquation, ics::NamedTuple) = error("check_initial_conditions(::GeometricEquation, ::NamedTuple) not implemented for ", typeof(equ), ".")
check_methods(equ::GeometricEquation, tspan, ics, params) = error("check_methods(::GeometricEquation, ::Tuple, ::NamedTuple, ::OptionalParameters) not implemented for ", typeof(equ), ".")

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
