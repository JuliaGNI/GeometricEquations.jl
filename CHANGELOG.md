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

The provisional `0.22.0` target was lowered to a patch: everything below is additive. No existing
type, constructor or exported name changes, and code written against 0.21.2 keeps working.

### New Features

- Concrete noise processes `WienerProcess(m)` and `GridProcess(ΔW, ΔZ)`, in
  `src/sdes/processes.jl`, together with the accessors `noise` and `noisedims` on the stochastic
  equations and their problems.

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

### Bug Fixes

### Breaking Changes

## Open Issues
