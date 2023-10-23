
using GeometricEquations: parameter_types
using Random


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

ode_ics   = (q=x₀,)
pode_ics  = (q=q₀, p=p₀)
iode_ics  = (q=q₀, p=p₀, λ=λ₀)
hode_ics  = (q=q₀, p=p₀)
lode_ics  = (q=q₀, p=p₀, λ=λ₀)

dae_ics   = (q=x₀, λ=λ₀)
pdae_ics  = (q=q₀, p=p₀, λ=λ₀)
hdae_ics  = (q=q₀, p=p₀, λ=λ₀)
idae_ics  = (q=q₀, p=p₀, λ=λ₀)
ldae_ics  = (q=q₀, p=p₀, λ=λ₀)

dae_ics_full   = (q=x₀, λ=λ₀, μ=λ₀)
pdae_ics_full  = (q=q₀, p=p₀, λ=λ₀, μ=λ₀)
hdae_ics_full  = (q=q₀, p=p₀, λ=λ₀, μ=λ₀)
idae_ics_full  = (q=q₀, p=p₀, λ=λ₀, μ=λ₀)
ldae_ics_full  = (q=q₀, p=p₀, λ=λ₀, μ=λ₀)

sde_ics   = (q=x₀,)
psde_ics  = (q=q₀, p=p₀)
spsde_ics = (q=q₀, p=p₀)

ode_params = (α=1,)
ode_param_types = parameter_types(ode_params)

sde_params = (λ = 2., μ = 1., noise_intensity = 0.1)
sde_param_types = parameter_types(sde_params)

x₁ₛ = [[ 0.5, 0.0],
        [ 0.0, 0.5],
        [-0.5, 0.0]]

q₁ₛ = [[0.5], [0.0], [-0.5]]
p₁ₛ = [[0.0], [0.5], [ 0.0]]
