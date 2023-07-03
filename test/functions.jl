
using Parameters: @unpack

using GeometricEquations: _idae_default_v̄, _ldae_default_v̄


function ode_v(ẋ, t, x, params)
    ẋ[1] = x[2]
    ẋ[2] = 2x[1]
end

ode_eqs = (ode_v,)


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

sode_eqs = (sode_v1, sode_v2)
sode_sols = (sode_q1, sode_q2)


function pode_v(v, t, q, p, params)
    v[1] = p[1]
end

function pode_f(f, t, q, p, params)
    f[1] = 2q[1]
end

function pode_h(t, q, p, params)
    p[1]^2/2 + cos(q[1])
end

pode_eqs = (pode_v, pode_f)


function hode_h(t, q, p, params)
    p[1]^2/2 + cos(q[1])
end

function hode_ω(ω, t, q, p, params)
    ω[1,1] = sin(q[1])
    ω[1,2] = 0
    ω[2,1] = 0
    ω[2,2] = p[1]
end

hode_eqs = (pode_v, pode_f, hode_h)


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
    v[1]^2/2 + cos(q[1])
end

iode_eqs = (iode_ϑ, iode_f, iode_g)


function lode_l(t, q, v, params)
    v[1]^2/2 - cos(q[1])
end

function lode_ω(ω, t, q, v, params)
    ω[1,1] = sin(q[1])
    ω[1,2] = 0
    ω[2,1] = 0
    ω[2,2] = v[1]
end

lode_eqs = (iode_ϑ, iode_f, iode_g, lode_ω, lode_l)


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
    p[1]^2/2 + q[1]^2/2
end

hdae_ω = hode_ω

idae_ϑ = iode_ϑ
idae_f = iode_f
idae_v = _idae_default_v̄

ldae_l = lode_l
ldae_ω = lode_ω
ldae_v = _ldae_default_v̄

idae_u(u, t, q, v, p, λ, params) = pdae_u(u, t, q, p, λ, params)
idae_g(g, t, q, v, p, λ, params) = pdae_g(g, t, q, p, λ, params)
idae_ϕ(ϕ, t, q, v, p, params) = pdae_ϕ(ϕ, t, q, p, params)
idae_ψ(ψ, t, q, v, p, q̇, ṗ, params) = pdae_ψ(ψ, t, q, p, q̇, ṗ, params)


function sde_v(v, t, q, params)
    @unpack λ = params
    v[1] = λ*q[1]
    v[2] = λ*q[2]
end

function sde_B(B::AbstractVector, t, q, params)
    @unpack μ = params
    B[1] = μ*q[1]
    B[2] = μ*q[2]
end

function sde_B(B::AbstractMatrix, t, q, params)
    @unpack μ = params
    B[1,:] = μ*q[1]
    B[2,:] = μ*q[2]
end


function psde_v(v, t, q, p, params)
    v[1] = p[1]
end

function psde_f(f, t, q, p, params)
    f[1] = - q[1]
end

function psde_B(B, t, q, p, params)
    @unpack noise_intensity = params
    B[1,1] = noise_intensity * p[1]
end

function psde_G(G, t, q, p, params)
    @unpack noise_intensity = params
    G[1,1] = - noise_intensity * q[1]
end


function spsde_v(v, t, q, p, params)
    v[1] = p[1]
end

function spsde_f1(f, t, q, p, params)
    f[1] = - q[1]
end

function spsde_f2(f, t, q, p, params)
    f[1] = 0
end

function spsde_B(B, t, q, p, params)
    @unpack noise_intensity = params
    B[1,1] = noise_intensity * p[1]
end

function spsde_G1(G, t, q, p, params)
    @unpack noise_intensity = params
    G[1,1] = - noise_intensity * q[1]
end

function spsde_G2(G, t, q, p, params)
    G[1,1] = 0
end
