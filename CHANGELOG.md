# Release Notes

All notable changes to GeometricEquations.jl.

This package is pre-1.0, so *every* minor release is potentially breaking in the sense of
[SemVer](https://semver.org) for `0.x` versions. The sections below name what actually
changed, so that a compat-only bump can be told apart from a rename or a change in results.

This file was started on 2026-08-31 and deliberately holds no entries. 68 versions were
released before it, the most recent `v0.21.2`, and none of them are written up here: the
record of that history is `git log` and the tags. It is named as a gap rather than
reconstructed, because a changelog assembled after the fact loses exactly the reasoning that
makes it worth keeping. The `[Unreleased]` target below is provisional — confirm it when the
first entry is written.

## [Unreleased] — targeting 0.21.3

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

### Bug Fixes

### Breaking Changes

## Open Issues

- `ntime(problem)` rounds up, and for some floating-point combinations of time span and time step
  it therefore reports one step more than the run actually needs: with `Δt = 0.01` over
  `(0.0, 0.1)` it gives 11 rather than 10, because `10Δt` is below the end time in binary floating
  point. `Δt = 0.01` over `(0.0, 1.0)`, and every other combination checked, is exact. This
  predates the noise processes and is not addressed here, but it is now easier to run into: a
  `GridProcess` sized from the same `nt` the time span was built from will be rejected as one
  increment short in exactly those cases.
