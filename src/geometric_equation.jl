
abstract type GeometricEquation{dType <: Number, tType <: Real} end

abstract type AbstractEquationODE{dType, tType} <: GeometricEquation{dType, tType} end
abstract type AbstractEquationPODE{dType, tType} <: GeometricEquation{dType, tType} end
abstract type AbstractEquationDAE{dType, tType} <: GeometricEquation{dType, tType} end
abstract type AbstractEquationPDAE{dType, tType} <: GeometricEquation{dType, tType} end
abstract type AbstractEquationSDE{dType, tType} <: GeometricEquation{dType, tType} end
abstract type AbstractEquationPSDE{dType, tType} <: GeometricEquation{dType, tType} end


Base.ndims(equ::GeometricEquation) = error("ndims() not implemented for ", typeof(equ), ".")
Base.axes(equ::GeometricEquation) = error("axes() not implemented for ", typeof(equ), ".")

GeometricBase.periodicity(equ::GeometricEquation) = error("periodicity() not implemented for ", typeof(equ), ".")

hassolution(::GeometricEquation) = false
hasvectorfield(::GeometricEquation) = false
hasprimary(::GeometricEquation) = false
hassecondary(::GeometricEquation) = false

hasinvariants(::GeometricEquation) = false
hasparameters(::GeometricEquation) = false
hasperiodicity(::GeometricEquation) = false

hashamiltonian(::GeometricEquation) = false
haslagrangian(::GeometricEquation) = false

functions(equ::GeometricEquation) = error("functions() not implemented for ", typeof(equ), ".")
solutions(equ::GeometricEquation) = error("solutions() not implemented for ", typeof(equ), ".")
