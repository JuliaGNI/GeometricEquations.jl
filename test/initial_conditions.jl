using GeometricEquations: parameter_types
using GeometricEquations: AlgebraicVariable, StateVariable, TimeVariable
using Random


Δt = .1
t₀ = 0.
t₁ = 1.
q₀ = [1.]
p₀ = [1.]
v₀ = [0.]
x₀ = [1., 1.]
x₁ = [2., 2.]
λ₀ = zero(v₀)
μ₀ = zero(v₀)

qₛ = rand(1,3,3)
pₛ = rand(1,3,3)
xₛ = rand(2,3,3)

x_arr = [x₀, rand(length(x₀))]
q_arr = [q₀, rand(length(q₀))]
p_arr = [p₀, rand(length(p₀))]
v_arr = [v₀, zero(v₀)]
λ_arr = [λ₀, zero(λ₀)]

x_sva = [StateVariable(x_arr[1]), StateVariable(x_arr[2])]
q_sva = [StateVariable(q_arr[1]), StateVariable(q_arr[2])]
p_sva = [StateVariable(p_arr[1]), StateVariable(p_arr[2])]
v_sva = [AlgebraicVariable(v_arr[1]), AlgebraicVariable(v_arr[2])]
λ_sva = [AlgebraicVariable(λ_arr[1]), AlgebraicVariable(λ_arr[2])]

ode_ics   = (q = StateVariable(x₀),)
pode_ics  = (q = StateVariable(q₀), p = StateVariable(p₀))
iode_ics  = (q = StateVariable(q₀), p = StateVariable(p₀), v = AlgebraicVariable(v₀))
hode_ics  = (q = StateVariable(q₀), p = StateVariable(p₀))
lode_ics  = (q = StateVariable(q₀), p = StateVariable(p₀), v = AlgebraicVariable(v₀))

ode_ics_raw   = (q = x₀,)
pode_ics_raw  = (q = q₀, p = p₀)
iode_ics_raw  = (q = q₀, p = p₀, v = v₀)
hode_ics_raw  = (q = q₀, p = p₀)
lode_ics_raw  = (q = q₀, p = p₀, v = v₀)

ode_ics_arr   = (q = x_arr,)
pode_ics_arr  = (q = q_arr, p = p_arr)
iode_ics_arr  = (q = q_arr, p = p_arr, v = v_arr)
hode_ics_arr  = (q = q_arr, p = p_arr)
lode_ics_arr  = (q = q_arr, p = p_arr, v = v_arr)

ode_ics_sva   = (q = x_sva,)
pode_ics_sva  = (q = q_sva, p = p_sva)
iode_ics_sva  = (q = q_sva, p = p_sva, v = v_sva)
hode_ics_sva  = (q = q_sva, p = p_sva)
lode_ics_sva  = (q = q_sva, p = p_sva, v = v_sva)

ode_ics_tpl = [
    (q = StateVariable(x_arr[1]),), (q = StateVariable(x_arr[2]),)
]

pode_ics_tpl = [
    (q = StateVariable(q_arr[1]), p = StateVariable(p_arr[1])),
    (q = StateVariable(q_arr[2]), p = StateVariable(p_arr[2])),
]

iode_ics_tpl = [
    (q = StateVariable(q_arr[1]), p = StateVariable(p_arr[1]), v = AlgebraicVariable(v_arr[1])),
    (q = StateVariable(q_arr[2]), p = StateVariable(p_arr[2]), v = AlgebraicVariable(v_arr[2])),
]

hode_ics_tpl = [
    (q = StateVariable(q_arr[1]), p = StateVariable(p_arr[1])),
    (q = StateVariable(q_arr[2]), p = StateVariable(p_arr[2])),
]

lode_ics_tpl = [
    (q = StateVariable(q_arr[1]), p = StateVariable(p_arr[1]), v = AlgebraicVariable(v_arr[1])),
    (q = StateVariable(q_arr[2]), p = StateVariable(p_arr[2]), v = AlgebraicVariable(v_arr[2])),
]

dae_ics   = (q = StateVariable(x₀), λ = AlgebraicVariable(λ₀))
pdae_ics  = (q = StateVariable(q₀), p = StateVariable(p₀), λ = AlgebraicVariable(λ₀))
hdae_ics  = (q = StateVariable(q₀), p = StateVariable(p₀), λ = AlgebraicVariable(λ₀))
idae_ics  = (q = StateVariable(q₀), p = StateVariable(p₀), v = AlgebraicVariable(v₀), λ = AlgebraicVariable(λ₀))
ldae_ics  = (q = StateVariable(q₀), p = StateVariable(p₀), v = AlgebraicVariable(v₀), λ = AlgebraicVariable(λ₀))

dae_ics_raw   = (q = x₀, λ = λ₀)
pdae_ics_raw  = (q = q₀, p = p₀, λ = λ₀)
hdae_ics_raw  = (q = q₀, p = p₀, λ = λ₀)
idae_ics_raw  = (q = q₀, p = p₀, v = v₀, λ = λ₀)
ldae_ics_raw  = (q = q₀, p = p₀, v = v₀, λ = λ₀)

dae_ics_full   = (q = StateVariable(x₀), λ = AlgebraicVariable(λ₀), μ = AlgebraicVariable(μ₀))
pdae_ics_full  = (q = StateVariable(q₀), p = StateVariable(p₀), λ = AlgebraicVariable(λ₀), μ = AlgebraicVariable(μ₀))
hdae_ics_full  = (q = StateVariable(q₀), p = StateVariable(p₀), λ = AlgebraicVariable(λ₀), μ = AlgebraicVariable(μ₀))
idae_ics_full  = (q = StateVariable(q₀), p = StateVariable(p₀), v = AlgebraicVariable(v₀), λ = AlgebraicVariable(λ₀), μ = AlgebraicVariable(μ₀))
ldae_ics_full  = (q = StateVariable(q₀), p = StateVariable(p₀), v = AlgebraicVariable(v₀), λ = AlgebraicVariable(λ₀), μ = AlgebraicVariable(μ₀))

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


dele_ics = (q̄ = StateVariable(x₀), q = StateVariable(x₁))
