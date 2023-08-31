```@meta
CurrentModule = GeometricEquations
```

# GeometricEquations

Documentation for [GeometricEquations](https://github.com/JuliaGNI/GeometricEquations.jl).

In *GeometricEquations.jl* we define three basic types of equations:
* ordinary differential equations (ODEs),
* differential algebraic equations (DAEs),
* stochastic differential equations (SDEs).

For each type, there are several subtypes
* standard equations ([`ODE`](@ref), [`DAE`](@ref), [`SDE`](@ref)),
* implicit equations ([`IODE`](@ref), [`IDAE`](@ref)),
* partitioned equations ([`PODE`](@ref), [`PDAE`](@ref), [`PSDE`](@ref)),
* Hamiltonian equations ([`HODE`](@ref), [`HDAE`](@ref)),
* Lagrangian equations ([`LODE`](@ref), [`LDAE`](@ref)),
* split equations ([`SODE`](@ref), [`SPDAE`](@ref), [`SPSDE`](@ref)).

Each equation holds a number of functions determining the vector field, constraints, and possibly additional information like periodicity, invariants and the Hamiltonian or Lagrangian.
In addition to each equation type, *GeometricEquations.jl* implements a corresponding problem type.
Each problem holds an equation, a time span `(t₀,t₁)` to integrate over, a time step to be used in the simulation, initial conditions and optionally parameters.

* The [`GeometricProblem`](@ref) type holds an equation together with initial conditions and parameters, the time span and the time step of a simulation.
* The [`GeometricEnsemble`](@ref) type holds an equation together with several initial conditions and/or parameters, etc.

GeometricEquations used to be part of [GeometricIntegrators](https://github.com/JuliaGNI/GeometricIntegrators.jl) and is primarily used to define equations and problems for GeometricIntegrators.
[GeometricProblems](https://github.com/JuliaGNI/GeometricProblems.jl) contains various predefined problems and [EulerLagrange](https://github.com/JuliaGNI/EulerLagrange.jl) can be used to generate code for equations from action principles.
