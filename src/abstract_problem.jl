abstract type AbstractProblem end

"Returns the parent equation object of the problem."
GeometricBase.equation(prob::AbstractProblem) = error("equation() not implemented for ", typeof(prob), ".")

"Returns a NamedTuple containing all functions (e.g. vector fields) provided by the equation."
GeometricBase.functions(prob::AbstractProblem) = functions(equation(prob))

"Returns a NamedTuple containing all solutions provided by the equation."
GeometricBase.solutions(prob::AbstractProblem) = solutions(equation(prob))

"Returns a NamedTuple containing all invariants provided by the equation."
GeometricBase.invariants(prob::AbstractProblem) = invariants(equation(prob))

GeometricBase.tspan(prob::AbstractProblem) = error("tspan() not implemented for ", typeof(prob), ".")
GeometricBase.tstep(prob::AbstractProblem) = error("tstep() not implemented for ", typeof(prob), ".")

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
