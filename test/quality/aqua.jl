using Aqua
using GeometricEquations
using Test

# Package-level quality assurance: method ambiguities, unbound type parameters, undefined
# exports, the agreement between `Project.toml` and `test/Project.toml`, stale dependencies,
# `[compat]` bounds, type piracy and persistent tasks.
Aqua.test_all(
    GeometricEquations;
    undefined_exports = (broken = true,),  # issue #39: AbstractEquationDELE is exported, not defined
    deps_compat = (broken = true,)        # issue #40: Random has no [compat] entry
)
