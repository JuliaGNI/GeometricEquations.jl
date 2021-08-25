
using GeometricEquations
using Test

struct TestEquation{DT,TT} <: Equation{DT,TT} end

@test_throws ErrorException ndims(TestEquation{Float64,Float64}())
