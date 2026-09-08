# Release Notes

All notable changes to GeometricEquations.jl.

This package is pre-1.0, so *every* minor release is potentially breaking in the sense of
[SemVer](https://semver.org) for `0.x` versions. The sections below name what actually
changed, so that a compat-only bump can be told apart from a rename or a change in results.

This file was started on 2026-08-31 and deliberately holds no entries for what preceded it. 68
versions were released before it, the most recent `v0.21.2`, and none of them are written up
here: the record of that history is `git log` and the tags. It is named as a gap rather than
reconstructed, because a changelog assembled after the fact loses exactly the reasoning that
makes it worth keeping.

## [Unreleased] — targeting 0.21.4

### Changed

- Every tracked source file is now Unicode NFC-normalised. The package was largely NFD — exactly
  four graphemes are rewritten, `ū` (215 times), `ḡ` (166), `ṗ` (32) and `ẋ` (10), each stored as a
  base letter followed by a combining mark rather than as one codepoint — which is an artefact of
  macOS rather than a decision.

  For ten of the fifteen files this is invisible: Julia's parser normalises identifiers to NFC, so
  the compiled symbols were already precomposed and nothing about dispatch, field names or method
  resolution changes. What changes is that the source now matches what a keyboard, an editor
  search, a `grep` pattern or an automated replacement produces. In an NFD file a pattern typed in
  NFC matches nothing at all, silently, and that is the failure this removes.

  The one behavioural difference is in `Base.show`, in the five remaining files. `HDAE`, `IDAE`,
  `LDAE` and `PDAE` print the literals `"   ū = "` and `"   ḡ = "`, and `DAE` prints the first of
  them; string literals are *not* parser-normalised, so those nine lines previously emitted
  decomposed bytes. Rendered output is unchanged to a reader, but code comparing it byte-for-byte
  against a string written in NFC would have failed before and succeeds now. No package in this
  tree that depends on GeometricEquations makes such a comparison, `src` and `test` contain no
  `Symbol("…")` literal at all, and the package has no doctests.

  The diff is mechanical and can be checked as such: every changed file is exactly the NFC
  normalisation of its predecessor. `scripts/verify_nfc.jl` asserts the standing invariant that
  every tracked `.jl`, `.md` and `.toml` file satisfies `s == Unicode.normalize(s, :NFC)`. Note
  that this is not the same as having no combining marks — `q̇`, `v̄`, `f̄`, `x̄` and `t̄` have no
  precomposed codepoint and remain two codepoints under NFC.

### Removed

- The `SPDAE` split partitioned DAE type, together with `SPDAEProblem` and `SPDAEEnsemble`.

  No code written against 0.21.3 is affected. `src/daes/spdae.jl` had not been part of the package
  since February 2024, when its `include` and all three exports were commented out rather than
  repaired, so `using GeometricEquations` never defined `SPDAE` and no caller could name it. No
  integrator ever used it either: in GeometricIntegrators the type appeared only in an export list,
  a `get_invariants` type union and the definition itself, all of which went when the equations
  moved to this package.

  What remained was a file stale against the refactor that moved initial conditions out of the
  equation types. `SPDAE` still carried `d`, `m`, `t₀`, `q₀`, `p₀`, `λ₀` and `μ₀` as fields and
  defined `Base.similar` on the equation rather than `initialstate` — an interface no other
  equation type has had for two years. A half-finished `StateVector` → `StateVariable` rename had
  additionally collapsed two `similar` methods onto one signature and left the surviving
  constructors forwarding to methods that no longer existed. Reviving the type would mean redoing
  that migration for something with no consumer, so it is removed instead;
  `git log --follow -- src/daes/spdae.jl` has it if it is ever wanted back.

## [0.21.3] — 2026-09-02

Nothing existing is renamed or removed, so code written against 0.21.2 keeps working — with two
caveats that stop this from being purely additive. `noise` and `noisedims` are newly *exported*,
so a package that both `using`s GeometricEquations and `using`s another package exporting a
generic of the same name will now see an ambiguity it did not see before, and has to qualify.
And building a stochastic problem whose driving `GridProcess` is too short for the run is now an
error where it was previously accepted (see below); code relying on that acceptance was already
heading for an out-of-bounds read partway through the integration.

### New Features

- Concrete noise processes `WienerProcess(m)` and `GridProcess(ΔW, ΔZ)`, together with the
  accessors `noise` and `noisedims` on the stochastic equations and their problems.

  `SDE`, `PSDE` and `SPSDE` have carried a `noise::AbstractStochasticProcess` field since 0.21,
  but there was nothing to put in it: `AbstractStochasticProcess` is an empty marker, so every
  caller had to declare its own type just to name the noise, and no consumer could ask that type
  anything. `GeometricProblems` had a `KuboNoise` for exactly this reason. A stochastic integrator
  in particular could not size its increment vectors, because nothing said how many Wiener
  processes an equation was driven by.

  For a caller this means the noise dimension now lives with the problem rather than having to be
  passed to the integrator alongside it, and `noisedims(problem)` answers it. `WienerProcess(m)`
  says only *which* noise drives the equation; the integrator draws the increments, since only it
  knows whether the scheme needs strong or weak ones. `GridProcess` prescribes the increments
  instead, which is what makes a run reproducible, lets two schemes be compared on one sample
  path, and — with zero increments — reduces a stochastic problem to its deterministic drift.

- A `GridProcess` that is too short for the run it is attached to is now rejected when the problem
  is built, rather than sending an integrator off the end of `ΔW` partway through. `ntime(process)`
  reports how many steps a process prescribes increments for.

## Open Issues

- `initial_multiplier` in `src/utils.jl` now has no caller in the package. `SPDAE` was its only
  one, and since that file stopped being compiled in February 2024 it has in practice had none for
  two years. It is unexported and its own testset in `test/utils_tests.jl` still passes, so it is
  left in place rather than removed alongside `SPDAE`; whether it is wanted for a future
  constrained equation type is a separate decision.

- `ntime(problem)` rounds up, and for some floating-point combinations of time span and time step
  it therefore reports one step more than the run actually needs: with `Δt = 0.01` over
  `(0.0, 0.1)` it gives 11 rather than 10, because `10Δt` is below the end time in binary floating
  point. `Δt = 0.01` over `(0.0, 1.0)`, and every other combination checked, is exact. This
  predates the noise processes and is not addressed here, but it is now easier to run into: a
  `GridProcess` sized from the same `nt` the time span was built from will be rejected as one
  increment short in exactly those cases.
