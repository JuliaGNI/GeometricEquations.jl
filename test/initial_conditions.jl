using GeometricEquations: parameter_types
using GeometricEquations: AlgebraicVariable, StateVariable
using Random


Δt = .1
t₀ = 0.
t₁ = 1.
q₀ = [1.]
p₀ = [1.]
x₀ = [1., 1.]
λ₀ = [0.]
μ₀ = zero(λ₀)

qₛ = rand(1,3,3)
pₛ = rand(1,3,3)
xₛ = rand(2,3,3)

ode_ics   = (q = StateVariable(x₀),)
pode_ics  = (q = StateVariable(q₀), p = StateVariable(p₀))
iode_ics  = (q = StateVariable(q₀), p = StateVariable(p₀), λ = AlgebraicVariable(λ₀))
hode_ics  = (q = StateVariable(q₀), p = StateVariable(p₀))
lode_ics  = (q = StateVariable(q₀), p = StateVariable(p₀), λ = AlgebraicVariable(λ₀))

ode_ics_raw   = (q = x₀,)
pode_ics_raw  = (q = q₀, p = p₀)
iode_ics_raw  = (q = q₀, p = p₀, λ = λ₀)
hode_ics_raw  = (q = q₀, p = p₀)
lode_ics_raw  = (q = q₀, p = p₀, λ = λ₀)

dae_ics   = (q = StateVariable(x₀), λ = AlgebraicVariable(λ₀))
pdae_ics  = (q = StateVariable(q₀), p = StateVariable(p₀), λ = AlgebraicVariable(λ₀))
hdae_ics  = (q = StateVariable(q₀), p = StateVariable(p₀), λ = AlgebraicVariable(λ₀))
idae_ics  = (q = StateVariable(q₀), p = StateVariable(p₀), λ = AlgebraicVariable(λ₀))
ldae_ics  = (q = StateVariable(q₀), p = StateVariable(p₀), λ = AlgebraicVariable(λ₀))

dae_ics_raw   = (q = x₀, λ = λ₀)
pdae_ics_raw  = (q = q₀, p = p₀, λ = λ₀)
hdae_ics_raw  = (q = q₀, p = p₀, λ = λ₀)
idae_ics_raw  = (q = q₀, p = p₀, λ = λ₀)
ldae_ics_raw  = (q = q₀, p = p₀, λ = λ₀)

dae_ics_full   = (q = StateVariable(x₀), λ = AlgebraicVariable(λ₀), μ = AlgebraicVariable(μ₀))
pdae_ics_full  = (q = StateVariable(q₀), p = StateVariable(p₀), λ = AlgebraicVariable(λ₀), μ = AlgebraicVariable(μ₀))
hdae_ics_full  = (q = StateVariable(q₀), p = StateVariable(p₀), λ = AlgebraicVariable(λ₀), μ = AlgebraicVariable(μ₀))
idae_ics_full  = (q = StateVariable(q₀), p = StateVariable(p₀), λ = AlgebraicVariable(λ₀), μ = AlgebraicVariable(μ₀))
ldae_ics_full  = (q = StateVariable(q₀), p = StateVariable(p₀), λ = AlgebraicVariable(λ₀), μ = AlgebraicVariable(μ₀))

sde_ics   = (q = StateVariable(x₀),)
psde_ics  = (q = StateVariable(q₀), p = StateVariable(p₀))
spsde_ics = (q = StateVariable(q₀), p = StateVariable(p₀))

sde_ics_raw   = (q = x₀,)
psde_ics_raw  = (q = q₀, p = p₀)
spsde_ics_raw = (q = q₀, p = p₀)

ode_params = (α=1,)
ode_param_types = parameter_types(ode_params)

sde_params = (λ = 2., μ = 1., noise_intensity = 0.1)
sde_param_types = parameter_types(sde_params)

x₁ₛ = [[ 0.5, 0.0],
        [ 0.0, 0.5],
        [-0.5, 0.0]]

q₁ₛ = [[0.5], [0.0], [-0.5]]
p₁ₛ = [[0.0], [0.5], [ 0.0]]
