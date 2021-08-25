
using GeometricEquations: function_dummy_v, initial_multiplier, symplectic_matrix
using Test

include("initial_conditions.jl")

@test function_dummy_v(t₀, q₀, λ₀) === nothing
@test function_dummy_v(t₀, q₀, p₀, λ₀) === nothing
zero_vec = [zeros(3) for i in 1:3]

@test initial_multiplier(zeros(3), zeros(3)) == zeros(3)
@test initial_multiplier(zeros(3), zero_vec) == zeros(3)
@test initial_multiplier(zero_vec, zeros(3)) == zero_vec
@test initial_multiplier(zero_vec, zero_vec) == zero_vec


Ω = symplectic_matrix(Float64, 4)

@test Ω == -Ω'
@test Ω == -inv(Ω)

@test symplectic_matrix(eltype(q₀), length(q₀)) == symplectic_matrix(t₀, q₀, p₀)
