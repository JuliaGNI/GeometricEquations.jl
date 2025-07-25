
using Parameters: @unpack

import GeometricBase: AbstractStochasticProcess
import GeometricEquations: _iode_default_v̄, _lode_default_v̄
import GeometricEquations: _iode_default_g, _lode_default_g
import GeometricEquations: _idae_default_v̄, _ldae_default_v̄

struct TestNoise <: AbstractStochasticProcess end

function ode_v(ẋ, t, x, params)
    ẋ[1] = x[2]
    ẋ[2] = 2x[1]
end

const ode_eqs = (ode_v,)
const ode_igs = (ode_v,)

function sode_v1(v, t, x, params)
    v[1] = x[2]
end

function sode_v2(v, t, x, params)
    v[2] = 2x[1]
end

function sode_q1(x, t, x̄, t̄, params)
    x[1] = x̄[1]
    x[2] = x̄[2]
end

function sode_q2(x, t, x̄, t̄, params)
    x[1] = x̄[1]
    x[2] = x̄[2]
end

const sode_eqs = (sode_v1, sode_v2)
const sode_sols = (sode_q1, sode_q2)
const sode_igs = (ode_v,)

function pode_v(v, t, q, p, params)
    v[1] = p[1]
end

function pode_f(f, t, q, p, params)
    f[1] = 2q[1]
end

function pode_h(t, q, p, params)
    p[1]^2 / 2 + cos(q[1])
end

const pode_eqs = (pode_v, pode_f)
const pode_igs = (pode_v, pode_f)

function hode_h(t, q, p, params)
    p[1]^2 / 2 + cos(q[1])
end

function hode_ω(ω, t, q, p, params)
    ω[1, 1] = sin(q[1])
    ω[1, 2] = 0
    ω[2, 1] = 0
    ω[2, 2] = p[1]
end

const hode_eqs = (pode_v, pode_f, hode_h)
const hode_igs = (pode_v, pode_f)

function iode_ϑ(p, t, q, v, params)
    p[1] = v[1]
end

function iode_f(f, t, q, v, params)
    f[1] = sin(q[1])
end

function iode_u(u, t, q, v, λ, params)
    u[1] = λ[1]
end

function iode_g(g, t, q, v, λ, params)
    g[1] = λ[1]
end

function iode_v(v, t, q, p, params)
    v[1] = p[1]
end

function iode_h(t, q, v, params)
    v[1]^2 / 2 + cos(q[1])
end

const iode_eqs_without_g = (iode_ϑ, iode_f)
const iode_eqs_default_g = (iode_ϑ, iode_f, _iode_default_g)
const iode_eqs = (iode_ϑ, iode_f, iode_g)
const iode_igs = (iode_v, iode_f)

function lode_l(t, q, v, params)
    v[1]^2 / 2 - cos(q[1])
end

function lode_ω(ω, t, q, v, params)
    ω[1, 1] = sin(q[1])
    ω[1, 2] = 0
    ω[2, 1] = 0
    ω[2, 2] = v[1]
end

const lode_v = iode_v
const lode_f = iode_f

const lode_eqs_without_g = (iode_ϑ, lode_f, lode_ω, lode_l)
const lode_eqs_default_g = (iode_ϑ, lode_f, _lode_default_g, lode_ω, lode_l)
const lode_eqs = (iode_ϑ, lode_f, iode_g, lode_ω, lode_l)
const lode_igs = (lode_v, lode_f)

function dae_v(v, t, x, params)
    v[1] = x[1]
    v[2] = x[2]
end

function dae_u(u, t, x, λ, params)
    u[1] = +λ[1]
    u[2] = -λ[1]
end

function dae_ū(u, t, x, λ, params)
    u[1] = +λ[1]
    u[2] = -λ[1]
end

function dae_ϕ(ϕ, t, x, params)
    ϕ[1] = x[2] - x[1]
end

function dae_ψ(ψ, t, x, v, params)
    ψ[1] = v[2] - v[1]
end

const dae_eqs = (dae_v, dae_u, dae_ϕ)
const dae_eqs_full = (dae_v, dae_u, dae_ϕ, dae_ū, dae_ψ)
const dae_igs = (dae_v,)

function pdae_v(v, t, q, p, params)
    v[1] = p[1]
end

function pdae_f(f, t, q, p, params)
    f[1] = q[1]
end

function pdae_p(p, t, q, v, params)
    p[1] = v[1]
end

function pdae_u(u, t, q, p, λ, params)
    u[1] = λ[1]
end

function pdae_g(g, t, q, p, λ, params)
    g[1] = λ[1]
end

function pdae_ϕ(ϕ, t, q, p, params)
    ϕ[1] = p[1] - q[1]
end

function pdae_ψ(ψ, t, q, p, q̇, ṗ, params)
    ψ[1] = μ[1] - λ[1]
end

function pdae_h(t, q, p, params)
    p[1]^2 / 2 + q[1]^2 / 2
end

const pdae_eqs = (pdae_v, pdae_f, pdae_u, pdae_g, pdae_ϕ)
const pdae_eqs_full = (pdae_v, pdae_f, pdae_u, pdae_g, pdae_ϕ, pdae_u, pdae_g, pdae_ψ)
const pdae_igs = (pdae_v, pdae_f)

const hdae_ω = hode_ω

const hdae_eqs = (pdae_eqs..., pdae_h)
const hdae_eqs_full = (pdae_eqs_full..., pdae_h)
const hdae_eqs_main = (pdae_eqs_full..., pdae_v, pdae_f, pdae_h)
const hdae_igs = (pdae_v, pdae_f)

const idae_ϑ = iode_ϑ
const idae_f = iode_f
const idae_v = _idae_default_v̄

idae_u(u, t, q, v, p, λ, params) = pdae_u(u, t, q, p, λ, params)
idae_g(g, t, q, v, p, λ, params) = pdae_g(g, t, q, p, λ, params)
idae_ϕ(ϕ, t, q, v, p, params) = pdae_ϕ(ϕ, t, q, p, params)
idae_ψ(ψ, t, q, v, p, q̇, ṗ, params) = pdae_ψ(ψ, t, q, p, q̇, ṗ, params)

const idae_eqs = (idae_ϑ, idae_f, idae_u, idae_g, idae_ϕ)
const idae_eqs_full = (idae_ϑ, idae_f, idae_u, idae_g, idae_ϕ, idae_u, idae_g, idae_ψ)
const idae_igs = (idae_v, idae_f)

const ldae_ω = lode_ω
const ldae_l = lode_l
const ldae_f = idae_f
const ldae_v = _ldae_default_v̄

const ldae_eqs = (idae_eqs..., ldae_ω, ldae_l)
const ldae_eqs_full = (idae_eqs_full..., ldae_ω, ldae_l)
const ldae_igs = (ldae_v, idae_f)

function sde_v(v, t, q, params)
    @unpack λ = params
    v[1] = λ * q[1]
    v[2] = λ * q[2]
end

function sde_B(B::AbstractVector, t, q, params)
    @unpack μ = params
    B[1] = μ * q[1]
    B[2] = μ * q[2]
end

function sde_B(B::AbstractMatrix, t, q, params)
    @unpack μ = params
    B[1, :] = μ * q[1]
    B[2, :] = μ * q[2]
end

const sde_eqs = (sde_v, sde_B)
const sde_igs = (sde_v,)

function psde_v(v, t, q, p, params)
    v[1] = p[1]
end

function psde_f(f, t, q, p, params)
    f[1] = -q[1]
end

function psde_B(B, t, q, p, params)
    @unpack noise_intensity = params
    B[1, 1] = noise_intensity * p[1]
end

function psde_G(G, t, q, p, params)
    @unpack noise_intensity = params
    G[1, 1] = -noise_intensity * q[1]
end

const psde_eqs = (psde_v, psde_f, psde_B, psde_G)
const psde_igs = (psde_v, psde_f)

function spsde_v(v, t, q, p, params)
    v[1] = p[1]
end

function spsde_f1(f, t, q, p, params)
    f[1] = -q[1]
end

function spsde_f2(f, t, q, p, params)
    f[1] = 0
end

function spsde_B(B, t, q, p, params)
    @unpack noise_intensity = params
    B[1, 1] = noise_intensity * p[1]
end

function spsde_G1(G, t, q, p, params)
    @unpack noise_intensity = params
    G[1, 1] = -noise_intensity * q[1]
end

function spsde_G2(G, t, q, p, params)
    G[1, 1] = 0
end

const spsde_eqs = (spsde_v, spsde_f1, spsde_f2, spsde_B, spsde_G1, spsde_G2)
const spsde_igs = (spsde_v, spsde_f1)

function dele_ld(t₀, t₁, q₀, q₁, params)
    h = (t₁ - t₀)
    q = (q₀[1] + q₁[1]) / 2
    v = (q₁[1] - q₀[1]) / h
    return h * (v^2 / 2 - cos(q))
end

function dele_d1ld(d, t₀, t₁, q₀, q₁, params)
    h = (t₁ - t₀)
    q = (q₀[1] + q₁[1]) / 2
    v = (q₁[1] - q₀[1]) / h
    d[1] = -v + h * sin(q) / 2
    return nothing
end

function dele_d2ld(d, t₀, t₁, q₀, q₁, params)
    h = (t₁ - t₀)
    q = (q₀[1] + q₁[1]) / 2
    v = (q₁[1] - q₀[1]) / h
    d[1] = v + h * sin(q) / 2
    return nothing
end
const

dele_eqs = (dele_ld, dele_d1ld, dele_d2ld)
const
dele_igs = ()
