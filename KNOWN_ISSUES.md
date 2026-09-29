# Known issues

### K1 · `initial_multiplier` in `src/utils.jl` now has no caller in the package.

- **location:** `src/utils.jl`
- **evidence:** `SPDAE` was its only
  one, and since that file stopped being compiled in February 2024 it has in practice had none for
  two years. It is unexported and its own testset in `test/utils_tests.jl` still passes, so it is
  left in place rather than removed alongside `SPDAE`; whether it is wanted for a future
  constrained equation type is a separate decision.
- **kind:** dead code
- **found:** 2026-09-08

### K2 · `ntime(problem)` rounds up, and for some floating-point combinations of time span and time step it therefore reports one step more than the run actually needs: with `Δt = 0.01` over `(0.0, 0.1)` it gives 11 rather than 10, because `10Δt` is below the end time in binary floating point.

- **location:** —
- **evidence:** `Δt = 0.01` over `(0.0, 1.0)`, and every other combination checked, is exact. This
  predates the noise processes and is not addressed here, but it is now easier to run into: a
  `GridProcess` sized from the same `nt` the time span was built from will be rejected as one
  increment short in exactly those cases.
- **kind:** defect
- **found:** 2026-09-02; one issue with StochasticIntegrators `K1`
