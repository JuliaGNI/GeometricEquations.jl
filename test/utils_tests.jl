
using GeometricEquations: function_dummy_v, initial_multiplier, symplectic_matrix,
                          promote_tspan, promote_tspan_and_tstep, parameter_types
using Test

include("initial_conditions.jl")

@test function_dummy_v(t₀, q₀, λ₀) === nothing
@test function_dummy_v(t₀, q₀, p₀, λ₀) === nothing


@test promote_tspan(0) == (0,0)
@test promote_tspan(1) == (0,1)
@test promote_tspan(0.0) == (0.0,0.0)
@test promote_tspan(1.0) == (0.0,1.0)
@test promote_tspan((0,   1.0)) == (0.0,1.0)
@test promote_tspan((0.0, 1  )) == (0.0,1.0)
@test promote_tspan([0,   1.0]) == (0.0,1.0)
@test promote_tspan([0.0, 1  ]) == (0.0,1.0)
@test promote_tspan([0,   1  ]) == (0,  1  )
@test promote_tspan([0.0, 1.0]) == (0.0,1.0)

@test_throws ErrorException promote_tspan([0.0])
@test_throws ErrorException promote_tspan([0.0, 0.5, 1.0])


@test promote_tspan_and_tstep((0,   1.0), 0  ) == ((0.0,1.0), 0.0)
@test promote_tspan_and_tstep((0,   1.0), 0.0) == ((0.0,1.0), 0.0)
@test promote_tspan_and_tstep((0.0, 1  ), 0  ) == ((0.0,1.0), 0.0)
@test promote_tspan_and_tstep((0.0, 1  ), 0.0) == ((0.0,1.0), 0.0)
@test promote_tspan_and_tstep((0,   1  ), 0  ) == ((0,  1  ), 0  )
@test promote_tspan_and_tstep((0,   1  ), 0.0) == ((0.0,1.0), 0.0)
@test promote_tspan_and_tstep((0.0, 1.0), 0  ) == ((0.0,1.0), 0.0)
@test promote_tspan_and_tstep((0.0, 1.0), 0.0) == ((0.0,1.0), 0.0)


params = (a = 1.0, b = 0, c = :c, d = [0,1])
types  = (a = Float64, b = Int, c = Symbol, d = Vector{Int})
@test parameter_types(params) == types


zero_vec = [zeros(3) for i in 1:3]

@test initial_multiplier(zeros(3), zeros(3)) == zeros(3)
@test initial_multiplier(zeros(3), zero_vec) == zeros(3)
@test initial_multiplier(zero_vec, zeros(3)) == zero_vec
@test initial_multiplier(zero_vec, zero_vec) == zero_vec


Ω = symplectic_matrix(Float64, 4)

@test Ω == -Ω'
@test Ω == -inv(Ω)

@test symplectic_matrix(eltype(q₀), length(q₀)) == symplectic_matrix(t₀, q₀, p₀)
