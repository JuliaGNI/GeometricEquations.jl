# GeometricEquations

[![Stable Docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaGNI.github.io/GeometricEquations.jl/stable)
[![Latest Docs](https://img.shields.io/badge/docs-latest-blue.svg)](https://JuliaGNI.github.io/GeometricEquations.jl/latest)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Build Status](https://github.com/JuliaGNI/GeometricEquations.jl/workflows/CI/badge.svg)](https://github.com/JuliaGNI/GeometricEquations.jl/actions?query=workflow:CI)
[![Coverage](https://codecov.io/gh/JuliaGNI/GeometricEquations.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaGNI/GeometricEquations.jl)

*GeometricEquations.jl* define different types of ordinary differential equations, differential algebraic and stochastic equations, that hold a number of functions determining the vector field, constraints, initial conditions, and possibly additional information like parameters, periodicity, invariants and the Hamiltonian or Lagrangian.


## Development

We are using git hooks, e.g., to enforce that all tests pass before pushing.
In order to activate these hooks, the following command must be executed once:
```
git config core.hooksPath .githooks
```
