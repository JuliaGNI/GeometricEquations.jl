module GeometricEquations

    using GeometricBase

    import Base: Callable
    import GeometricBase: equation, equations, functions, solutions, invariants, parameters, periodicity

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

    export AbstractProblem, GeometricProblem

    export AbstractProblemODE, AbstractProblemPODE,
           AbstractProblemDAE, AbstractProblemPDAE

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

    export nsamples, nconstraints
    export equation, equations, functions, solutions
    export invariants, parameters, periodicity
    export initial_conditions

    export hassolution, hasvectorfield, hasprimary, hassecondary,
           hasinvariants, hasparameters, hasperiodicity,
           hashamiltonian, haslagrangian

    export tspan, tstep, tbegin, tend

    export datatype, timetype, arrtype, equtype
           
    
    include("utils.jl")

    include("geometric_equation.jl")

    include("ode.jl")
    include("hode.jl")
    include("iode.jl")
    include("lode.jl")
    include("pode.jl")

    include("dae.jl")
    include("hdae.jl")
    include("idae.jl")
    include("ldae.jl")
    include("pdae.jl")

    include("sode.jl")
    include("spdae.jl")

    include("sde.jl")
    include("psde.jl")
    include("spsde.jl")

    include("abstract_problem.jl")
    include("geometric_problem.jl")

    include("conversion.jl")
    
end
