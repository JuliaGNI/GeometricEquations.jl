using GeometricEquations.Tests
using Test


@test_nowarn Tests.ExponentialGrowth.odeproblem()

@test_nowarn Tests.HarmonicOscillator.odeproblem()
@test_nowarn Tests.HarmonicOscillator.hodeproblem()
@test_nowarn Tests.HarmonicOscillator.iodeproblem()
@test_nowarn Tests.HarmonicOscillator.lodeproblem()
@test_nowarn Tests.HarmonicOscillator.podeproblem()
@test_nowarn Tests.HarmonicOscillator.sodeproblem()

@test_nowarn Tests.HarmonicOscillator.daeproblem()
@test_nowarn Tests.HarmonicOscillator.hdaeproblem()
@test_nowarn Tests.HarmonicOscillator.idaeproblem()
@test_nowarn Tests.HarmonicOscillator.ldaeproblem()
@test_nowarn Tests.HarmonicOscillator.pdaeproblem()
