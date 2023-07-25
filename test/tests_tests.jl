using GeometricEquations.Tests
using Test


@test_nowarn Tests.ExponentialGrowth.odeproblem()
@test_nowarn Tests.ExponentialGrowth.odeensemble()
