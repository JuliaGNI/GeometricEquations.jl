"""
Abstract type that describes a generic interface for different problem types.
"""
abstract type GeometricProblem{superType} <: AbstractProblem end

"Returns the parent equation object of the problem."
GeometricBase.equation(prob::GeometricProblem) = error("equation() not implemented for ", typeof(prob), ".")

"Returns a NamedTuple containing all functions (e.g. vector fields) provided by the equation."
GeometricBase.functions(prob::GeometricProblem) = functions(equation(prob))

"Returns a NamedTuple containing all solutions provided by the equation."
GeometricBase.solutions(prob::GeometricProblem) = solutions(equation(prob))

"Returns a NamedTuple containing all invariants provided by the equation."
GeometricBase.invariants(prob::GeometricProblem) = invariants(equation(prob))

GeometricBase.timespan(prob::GeometricProblem) = error("timespan() not implemented for ", typeof(prob), ".")
GeometricBase.timestep(prob::GeometricProblem) = error("timestep() not implemented for ", typeof(prob), ".")

GeometricBase.parameters(prob::GeometricProblem) = error("parameters() not implemented for ", typeof(prob), ".")
GeometricBase.periodicity(prob::GeometricProblem) = periodicity(equation(prob))

hassolution(prob::GeometricProblem) = hassolution(equation(prob))
hasvectorfield(prob::GeometricProblem) = hasvectorfield(equation(prob))
hasinitialguess(prob::GeometricProblem) = hasinitialguess(equation(prob))
hasprimary(prob::GeometricProblem) = hasprimary(equation(prob))
hassecondary(prob::GeometricProblem) = hassecondary(equation(prob))

hasinvariants(prob::GeometricProblem) = hasinvariants(equation(prob))
hasparameters(prob::GeometricProblem) = hasparameters(equation(prob))
hasperiodicity(prob::GeometricProblem) = hasperiodicity(equation(prob))

hashamiltonian(prob::GeometricProblem) = hashamiltonian(equation(prob))
haslagrangian(prob::GeometricProblem) = haslagrangian(equation(prob))
