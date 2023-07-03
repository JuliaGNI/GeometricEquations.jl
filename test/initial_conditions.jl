
using GeometricEquations: parameter_types


Δt = .1
t₀ = 0.
t₁ = 1.
q₀ = [1.]
p₀ = [1.]
x₀ = [1., 1.]
λ₀ = [0.]

qₛ = rand(1,3,3)
pₛ = rand(1,3,3)
xₛ = rand(2,3,3)

ode_params = (α=1,)
ode_param_types = parameter_types(ode_params)

sde_params = (λ = 2., μ = 1., noise_intensity = 0.1)
sde_param_types = parameter_types(sde_params)

x₁ₛ = [[ 0.5, 0.0],
        [ 0.0, 0.5],
        [-0.5, 0.0]]

q₁ₛ = [[0.5], [0.0], [-0.5]]
p₁ₛ = [[0.0], [0.5], [ 0.0]]
