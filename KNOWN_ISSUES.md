# Known issues

### K1 · `initial_multiplier` in `src/utils.jl` has no caller in the package.

- location: `src/utils.jl:80`
- evidence: no file in `src/` calls it; `git grep -n initial_multiplier -- src test` gives only
  its four method definitions in `src/utils.jl` (lines 80, 85, 90 and 95) and its tests in
  `test/utils.jl` (lines 2 and 40–43). It is unexported, and its testset "Utility Functions" in
  `test/utils.jl` passes (30 of 30). It stays in place rather than going with `SPDAE`; whether it
  is wanted for a future constrained equation type is a separate decision.
- kind: dead code
- found: 2026-09-08, in GeometricEquations #36, which removes `SPDAE`, its one caller; the file
  `src/daes/spdae.jl` is not compiled from February 2024 until that removal

### K2 · `ntime(problem)` rounds up, so for some floating-point time spans and time steps it reports one step more than the run needs.

- location: `src/problems/equation_problem.jl:138`
- evidence: with `Δt = 0.01` over `(0.0, 0.1)` it gives 11 rather than 10, because `10Δt` is
  below the end time in binary floating point: `10 * big(0.01) < big(0.1)` is `true`, and
  `div(0.1, 0.01, RoundUp)` is `11.0`, although the rounded product `10 * 0.01 == 0.1` is also
  `true`. `Δt = 0.01` over `(0.0, 1.0)`, and every other combination checked, is exact. A
  `GridProcess` sized from the same `nt` that built the time span is rejected as one increment
  short in exactly those cases, and the noise processes make this case easier to reach.
- kind: defect
- found: 2026-09-02; one issue with StochasticIntegrators `K1`; the rounding is older than the
  noise processes
