
abstract type AbstractProblem{dType <: Number, tType <: Real, arrayType <: AbstractArray{dType}} end

"Returns the parent equation object of the problem."
equation(prob::AbstractProblem) = error("equation() not implemented for ", typeof(prob), ".")

"Returns a NamedTuple containing all functions (e.g. vector fields) provided by the equation."
functions(prob::AbstractProblem) = functions(equation(prob))

"Returns a NamedTuple containing all solutions provided by the equation."
solutions(prob::AbstractProblem) = solutions(equation(prob))

"Returns a NamedTuple containing all invariants provided by the equation."
invariants(prob::AbstractProblem) = invariants(equation(prob))

tspan(prob::AbstractProblem) = error("tspan() not implemented for ", typeof(prob), ".")
tstep(prob::AbstractProblem) = error("tstep() not implemented for ", typeof(prob), ".")
tbegin(prob::AbstractProblem) = tspan(prob)[begin]
tend(prob::AbstractProblem) = tspan(prob)[end]

GeometricBase.parameters(prob::AbstractProblem) = error("parameters() not implemented for ", typeof(prob), ".")
GeometricBase.periodicity(prob::AbstractProblem) = periodicity(equation(prob))

hassolution(prob::AbstractProblem) = hassolution(equation(prob))
hasvectorfield(prob::AbstractProblem) = hasvectorfield(equation(prob))
hasprimary(prob::AbstractProblem) = hasprimary(equation(prob))
hassecondary(prob::AbstractProblem) = hassecondary(equation(prob))

hasinvariants(prob::AbstractProblem) = hasinvariants(equation(prob))
hasparameters(prob::AbstractProblem) = hasparameters(equation(prob))
hasperiodicity(prob::AbstractProblem) = hasperiodicity(equation(prob))

hashamiltonian(prob::AbstractProblem) = hashamiltonian(equation(prob))
haslagrangian(prob::AbstractProblem) = haslagrangian(equation(prob))
