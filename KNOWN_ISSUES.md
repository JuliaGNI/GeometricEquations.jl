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
