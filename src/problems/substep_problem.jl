
struct SubstepProblem{superType <: GeometricEquation, dType <: Number, tType <: Real,
        problemType <: EquationProblem,
        functionsType <: NamedTuple,
        solutionsType <: NamedTuple,
        cType <: Real} <: GeometricProblem{superType}

    problem::problemType
    functions::functionsType
    solutions::solutionsType
    coefficient::cType
    index::Int

    function SubstepProblem(problem::SODEProblem, c::Real, i::Int)
        keys = ()
        fncs = ()
        for (key, fnc) in pairs(functions(problem))
            keys = (keys..., key)
            fncs = (fncs..., fnc[i])
        end
        fnc_tuple = NamedTuple{keys}(fncs)

        keys = ()
        sols = ()
        for (key, sol) in pairs(solutions(problem))
            keys = (keys..., key)
            sols = (sols..., sol[i])
        end
        sol_tuple = NamedTuple{keys}(sols)

        new{SODE, datatype(problem), timetype(problem), typeof(problem), typeof(fnc_tuple), typeof(sol_tuple), typeof(c)}(problem, fnc_tuple, sol_tuple, c, i)
    end
end

problem(ssp::SubstepProblem) = ssp.problem
coefficient(ssp::SubstepProblem) = ssp.coefficient

@inline Base.parent(ssp::SubstepProblem) = problem(ssp)

@inline initial_conditions(ssp::SubstepProblem) = initial_conditions(problem(ssp))

@inline GeometricBase.functions(ssp::SubstepProblem) = ssp.functions
@inline GeometricBase.solutions(ssp::SubstepProblem) = ssp.solutions

@inline GeometricBase.datatype(ssp::SubstepProblem) = datatype(problem(ssp))
@inline GeometricBase.timetype(ssp::SubstepProblem) = timetype(problem(ssp))
@inline GeometricBase.arrtype(ssp::SubstepProblem) = arrtype(problem(ssp))
@inline GeometricBase.equtype(ssp::SubstepProblem) = equtype(problem(ssp))

@inline GeometricBase.ndims(ssp::SubstepProblem) = ndims(problem(ssp))
@inline GeometricBase.ntime(ssp::SubstepProblem) = ntime(problem(ssp))
@inline GeometricBase.timespan(ssp::SubstepProblem) = timespan(problem(ssp))
@inline GeometricBase.timestep(ssp::SubstepProblem) = coefficient(ssp) * timestep(problem(ssp))

@inline GeometricBase.equation(ssp::SubstepProblem) = equation(problem(ssp))
@inline GeometricBase.invariants(ssp::SubstepProblem) = invariants(problem(ssp))
@inline GeometricBase.parameters(ssp::SubstepProblem) = parameters(problem(ssp))
@inline GeometricBase.nsamples(ssp::SubstepProblem) = nsamples(problem(ssp))
