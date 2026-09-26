using Aqua
using GeometricEquations
using Test

# Package-level quality assurance: type piracy, method ambiguities, stale and duplicated
# dependencies, undefined exports, unbound type parameters, `Project.toml` validity.
Aqua.test_all(
    GeometricEquations;
    undefined_exports = (broken = true,),  # issue #39: AbstractEquationDELE is exported, not defined
    deps_compat = (broken = true,)        # issue #40: Random has no [compat] entry
)
