@doc raw"""
# Exponential Growth

"""
module ExponentialGrowth

    using ...GeometricEquations

    using Parameters

    export odeproblem, odeensemble, exact_solution

    const x₀ = [1.0]
    const Δt = 0.1
    const tbeg = 0.0
    const tend = 10.0

    const ics = [(q = rand(1),), (q = rand(1),), (q = rand(1),)]

    const k = 1.0

    const default_parameters = (k=k,)
    
    function vectorfield(v, t, x, params)
        @unpack k = params
        v .= k .* x
        nothing
    end

    function solution(x₁, t₁, x₀, t₀, params)
        @unpack k = params
        x₁ .= x₀ .* exp(k * (t₁-t₀))
        return x₁
    end

    function odeproblem(x₀ = x₀; parameters = default_parameters, tbegin = tbeg, tend = tend, Δt = Δt)
        ODEProblem(vectorfield, (tbegin, tend), Δt, x₀; parameters = parameters)
    end

    function odeensemble(ics = ics; parameters = default_parameters, tbegin = tbeg, tend = tend, Δt = Δt)
        ODEEnsemble(vectorfield, (tbegin, tend), Δt, ics; parameters = parameters)
    end

    function exact_solution(prob::ODEProblem)
        sol = GeometricSolution(prob)
        for n in eachtimestep(sol)
            solution(sol.q[n], sol.t[n], sol.q[0], sol.t[0], parameters(prob))
        end
        return sol
    end

    function exact_solution(probs::ODEEnsemble)
        sols = EnsembleSolution(probs)
        for (sol,prob) in zip(sols.s, probs)
            for n in eachtimestep(sol)
                solution(sol.q[n], sol.t[n], sol.q[0], sol.t[0], parameters(prob))
            end
        end
        return sols
    end

end
