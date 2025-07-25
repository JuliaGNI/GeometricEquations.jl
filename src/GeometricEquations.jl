module GeometricEquations

using GeometricBase

import Base: Callable

import GeometricBase: datatype, timetype, arrtype, equtype
import GeometricBase: initialtime, finaltime, timespan, timestep
import GeometricBase: equation, equations, functions, solutions, initialguess
import GeometricBase: invariants, parameters, periodicity
import GeometricBase: ntime, nsamples, nconstraints, nsteps
import GeometricBase: AbstractStateVariable, AbstractStochasticProcess

export NullInvariants, NullParameters, NullPeriodicity

export OptionalAbstractArray, OptionalArray,
       OptionalFunction, OptionalNamedTuple,
       OptionalInvariants, OptionalParameters,
       OptionalPeriodicity

export AlgebraicVariable, StateVariable

export GeometricEquation
export AbstractEquationODE, AbstractEquationPODE,
       AbstractEquationDAE, AbstractEquationPDAE,
       AbstractEquationSDE, AbstractEquationPSDE,
       AbstractEquationDELE

export GeometricProblem, EquationProblem, SubstepProblem
export EnsembleProblem

export AbstractProblemODE, AbstractProblemPODE, AbstractProblemIODE, AbstractProblemHODE,
       AbstractProblemLODE,
       AbstractProblemDAE, AbstractProblemPDAE, AbstractProblemIDAE,
       AbstractProblemSDE, AbstractProblemPSDE, AbstractProblemSPSDE,
       AbstractProblemDELE

export ODE, IODE, PODE, HODE, LODE, SODE
export DAE, IDAE, PDAE, HDAE, LDAE#, SPDAE
export SDE, PSDE, SPSDE
export DELE

export ODEProblem, IODEProblem, PODEProblem,
       HODEProblem, LODEProblem, SODEProblem
export DAEProblem, IDAEProblem, PDAEProblem,
       HDAEProblem, LDAEProblem#, SPDAEProblem
export SDEProblem, PSDEProblem, SPSDEProblem
export DELEProblem

export ODEEnsemble, IODEEnsemble, PODEEnsemble,
       HODEEnsemble, LODEEnsemble, SODEEnsemble
export DAEEnsemble, IDAEEnsemble, PDAEEnsemble,
       HDAEEnsemble, LDAEEnsemble#, SPDAEEnsemble
export SDEEnsemble, PSDEEnsemble, SPSDEEnsemble
export DELEEnsemble

export datatype, timetype, arrtype, equtype
export initialtime, finaltime, timespan, timestep
export problem, equation, equations, functions, solutions, initialguess
export invariants, parameters, periodicity, getperiodicity
export initial_conditions, compute_vectorfields!
export ntime, nsamples, nconstraints, nsteps

export hassolution, hasvectorfield, hasinitialguess, hasprimary, hassecondary,
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
# include("daes/spdae.jl")

include("sdes/sde.jl")
include("sdes/psde.jl")
include("sdes/spsde.jl")

include("discrete/dele.jl")

include("problems/substep_problem.jl")

include("conversion.jl")

# Union types for problems that share a common interface

const AbstractProblemODE{DT, TT} = Union{
    ODEProblem{DT, TT}, SubstepProblem{DT, TT}, DAEProblem{DT, TT}}
const AbstractProblemDAE{DT, TT} = Union{DAEProblem{DT, TT}}

const AbstractProblemPODE{DT, TT} = Union{
    PODEProblem{DT, TT}, HODEProblem{DT, TT}, PDAEProblem{DT, TT}, HDAEProblem{DT, TT}}
const AbstractProblemIODE{DT, TT} = Union{
    IODEProblem{DT, TT}, LODEProblem{DT, TT}, IDAEProblem{DT, TT}, LDAEProblem{DT, TT}}

const AbstractProblemHODE{DT, TT} = Union{HODEProblem{DT, TT}, HDAEProblem{DT, TT}}
const AbstractProblemLODE{DT, TT} = Union{LODEProblem{DT, TT}, LDAEProblem{DT, TT}}

const AbstractProblemPDAE{DT, TT} = Union{PDAEProblem{DT, TT}, HDAEProblem{DT, TT}}
const AbstractProblemIDAE{DT, TT} = Union{IDAEProblem{DT, TT}, LDAEProblem{DT, TT}}

const AbstractProblemSDE{DT, TT} = Union{SDEProblem{DT, TT}}
const AbstractProblemPSDE{DT, TT} = Union{PSDEProblem{DT, TT}}
const AbstractProblemSPSDE{DT, TT} = Union{SPSDEProblem{DT, TT}}

const AbstractProblemDELE{DT, TT} = Union{DELEProblem{DT, TT}}

include("tests/Tests.jl")

end
