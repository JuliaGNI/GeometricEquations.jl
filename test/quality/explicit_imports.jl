using ExplicitImports
using GeometricEquations
using Test

# Stale explicit imports, names imported from or qualified through a module that does not own
# them, and self-qualified accesses.
test_explicit_imports(
    GeometricEquations;
    # the package relies on `using` of its dependencies throughout
    no_implicit_imports = false,
    # several imported names are not public in their owner (GeometricBase, and `Base.Callable`)
    all_explicit_imports_are_public = false,
    # several qualified names are not public in their owner (GeometricBase, and
    # `Base.AbstractCartesianIndex`)
    all_qualified_accesses_are_public = false
)
