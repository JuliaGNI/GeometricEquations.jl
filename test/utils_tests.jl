
using GeometricEquations: function_v_dummy, get_λ₀, symplectic_matrix
using Test

include("initial_conditions.jl")

zero_vec = [zeros(3) for i in 1:3]

@test function_v_dummy(t₀, q₀, p₀, λ₀) === nothing

@test get_λ₀(zeros(3), zeros(3)) == zeros(3)
@test get_λ₀(zeros(3), zero_vec) == zeros(3)
@test get_λ₀(zero_vec, zeros(3)) == zero_vec
@test get_λ₀(zero_vec, zero_vec) == zero_vec

Ω = symplectic_matrix(Float64, 4)

@test Ω == -Ω'
@test Ω == -inv(Ω)

@test symplectic_matrix(eltype(q₀), length(q₀)) == symplectic_matrix(t₀, q₀, p₀)
