module GeometricEquations

    using GeometricBase

    import Base: Callable

    import GeometricBase: datatype, timetype, arrtype, equtype
    import GeometricBase: tspan, tstep, tbegin, tend, timestep
    import GeometricBase: equation, equations, functions, solutions, invariants, parameters, periodicity
    import GeometricBase: ntime, nsamples, nconstraints, nsteps

    export NullInvariants, NullParameters, NullPeriodicity

    export OptionalAbstractArray, OptionalArray,
           OptionalFunction, OptionalNamedTuple,
           OptionalInvariants, OptionalParameters,
           OptionalPeriodicity
    
    export State, StateVector
    
    export GeometricEquation
    export AbstractEquationODE, AbstractEquationPODE,
           AbstractEquationDAE, AbstractEquationPDAE,
           AbstractEquationSDE, AbstractEquationPSDE

    export GeometricProblem, EquationProblem, SubstepProblem
    export EnsembleProblem

    export AbstractProblemODE, AbstractProblemPODE, AbstractProblemIODE, AbstractProblemHODE, AbstractProblemLODE,
           AbstractProblemDAE, AbstractProblemPDAE, AbstractProblemIDAE,
           AbstractProblemSDE, AbstractProblemPSDE, AbstractProblemSPSDE

    export ODE, IODE, PODE, HODE, LODE, SODE
    export DAE, IDAE, PDAE, HDAE, LDAE, SPDAE
    export SDE, PSDE, SPSDE

    export ODEProblem,  IODEProblem, PODEProblem,
           HODEProblem, LODEProblem, SODEProblem
    export DAEProblem,  IDAEProblem, PDAEProblem,
           HDAEProblem, LDAEProblem, SPDAEProblem
    export SDEProblem,  PSDEProblem, SPSDEProblem

    export ODEEnsemble,  IODEEnsemble, PODEEnsemble,
           HODEEnsemble, LODEEnsemble, SODEEnsemble
    export DAEEnsemble,  IDAEEnsemble, PDAEEnsemble,
           HDAEEnsemble, LDAEEnsemble, SPDAEEnsemble
    export SDEEnsemble,  PSDEEnsemble, SPSDEEnsemble

    export datatype, timetype, arrtype, equtype
    export tspan, tstep, tbegin, tend, timestep
    export problem, equation, equations, functions, solutions
    export invariants, parameters, periodicity
    export initial_conditions
    export ntime, nsamples, nconstraints, nsteps

    export hassolution, hasvectorfield, hasprimary, hassecondary,
           hasinvariants, hasparameters, hasperiodicity,
           hashamiltonian, haslagrangian

    
    include("utils.jl")

    include("geometric_equation.jl")
    include("geometric_problem.jl")

    include("problems/ensemble_problem.jl")
    include("problems/equation_problem.jl")

    include("odes/ode.jl")
    include("odes/hode.jl")
    include("odes/iode.jl")
    include("odes/lode.jl")
    include("odes/pode.jl")

    include("daes/dae.jl")
    include("daes/hdae.jl")
    include("daes/idae.jl")
    include("daes/ldae.jl")
    include("daes/pdae.jl")

    include("odes/sode.jl")
    include("daes/spdae.jl")

    include("sdes/sde.jl")
    include("sdes/psde.jl")
    include("sdes/spsde.jl")

    include("problems/substep_problem.jl")
    
    include("conversion.jl")

    
    # Union types for problems that share a common interface

    const AbstractProblemODE = Union{ODEProblem, SubstepProblem, DAEProblem}
    const AbstractProblemDAE = Union{DAEProblem}

    const AbstractProblemPODE = Union{PODEProblem, HODEProblem, PDAEProblem, HDAEProblem}
    const AbstractProblemIODE = Union{IODEProblem, LODEProblem, IDAEProblem, LDAEProblem}

    const AbstractProblemHODE = Union{HODEProblem, HDAEProblem}
    const AbstractProblemLODE = Union{LODEProblem, LDAEProblem}

    const AbstractProblemPDAE = Union{PDAEProblem, HDAEProblem}
    const AbstractProblemIDAE = Union{IDAEProblem, LDAEProblem}

    const AbstractProblemSDE = Union{SDEProblem}
    const AbstractProblemPSDE = Union{PSDEProblem}
    const AbstractProblemSPSDE = Union{SPSDEProblem}


    include("tests/Tests.jl")
    
end
