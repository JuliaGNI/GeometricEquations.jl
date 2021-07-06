
abstract type Equation{dType <: Number, tType <: Real} end

abstract type AbstractEquationODE{dType, tType} <: Equation{dType, tType} end
abstract type AbstractEquationPODE{dType, tType} <: Equation{dType, tType} end
abstract type AbstractEquationDAE{dType, tType} <: Equation{dType, tType} end
abstract type AbstractEquationPDAE{dType, tType} <: Equation{dType, tType} end
abstract type AbstractEquationSDE{dType, tType} <: Equation{dType, tType} end
abstract type AbstractEquationPSDE{dType, tType} <: Equation{dType, tType} end


Base.eltype(::Equation{DT,TT}) where {DT,TT} = DT
timetype(::Equation{DT,TT}) where {DT,TT} = TT

Base.ndims(equ::Equation) = error("ndims() not implemented for ", typeof(equ), ".")

periodicity(equ::Equation) = error("periodicity() not implemented for ", typeof(equ), ".")

hasinvariants(::Equation) = false
hasparameters(::Equation) = false
hasperiodicity(::Equation) = false
hashamiltonian(::Equation) = false

get_functions(equ::Equation) = error("get_functions() not implemented for ", typeof(equ), ".")
get_solutions(equ::Equation) = error("get_solutions() not implemented for ", typeof(equ), ".")


const OptionalArray{arrayType} = Union{Nothing, arrayType} where {arrayType <: AbstractArray}

const OptionalAbstractArray = Union{Nothing, AbstractArray}
const OptionalFunction      = Union{Nothing, Function}
const OptionalNamedTuple    = Union{Nothing, NamedTuple}


const State{DT <: Number} = AbstractArray{DT}
const StateVector{DT,VT} = VT where {DT, VT <: AbstractVector{<:State{DT}}}

Base.zero(X::ST) where {DT, VT, ST <: StateVector{DT,VT}} = VT[zero(x) for x in X]
