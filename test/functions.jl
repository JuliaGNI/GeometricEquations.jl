
using Parameters: @unpack

using GeometricEquations: _idae_default_v̄, _ldae_default_v̄

function ode_v(t, x, ẋ, params)
    ẋ[1] = x[2]
    ẋ[2] = 2x[1]
end


function sode_v1(t, x, v, params)
    v[1] = x[2]
end

function sode_v2(t, x, v, params)
    v[2] = 2x[1]
end

function sode_q1(t, x̄, x, h, params)
    x[1] = x̄[1]
    x[2] = x̄[2]
end

function sode_q2(t, x̄, x, h, params)
    x[1] = x̄[1]
    x[2] = x̄[2]
end

sode_v = (sode_v1, sode_v2)
sode_q = (sode_q1, sode_q2)


function pode_v(t, q, p, v, params)
    v[1] = p[1]
end

function pode_f(t, q, p, f, params)
    f[1] = 2q[1]
end

function pode_h(t, q, p, params)
    p[1]^2/2 + cos(q[1])
end

pode_eqs = (pode_v, pode_f)


function hode_h(t, q, p, params)
    p[1]^2/2 + cos(q[1])
end

function hode_ω(t, q, p, ω, params)
    ω[1,1] = sin(q[1])
    ω[1,2] = 0
    ω[2,1] = 0
    ω[2,2] = p[1]
end

hode_eqs = (pode_v, pode_f, hode_h)


function iode_ϑ(t, q, v, p, params)
    p[1] = v[1]
end

function iode_f(t, q, v, f, params)
    f[1] = sin(q[1])
end

function iode_u(t, q, λ, u, params)
    u[1] = λ[1]
end

function iode_g(t, q, λ, g, params)
    g[1] = λ[1]
end

function iode_v(t, q, p, v, params)
    v[1] = p[1]
end

function iode_h(t, q, v, params)
    v[1]^2/2 + cos(q[1])
end

iode_eqs = (iode_ϑ, iode_f, iode_g)


function lode_l(t, q, v, params)
    v[1]^2/2 - cos(q[1])
end

function lode_ω(t, q, v, ω, params)
    ω[1,1] = sin(q[1])
    ω[1,2] = 0
    ω[2,1] = 0
    ω[2,2] = v[1]
end

lode_eqs = (iode_ϑ, iode_f, iode_g, lode_ω, lode_l)


function dae_v(t, x, v, params)
    v[1] = x[1]
    v[2] = x[2]
end

function dae_u(t, x, λ, u, params)
    u[1] = +λ[1]
    u[2] = -λ[1]
end

function dae_ū(t, x, λ, u, params)
    u[1] = +λ[1]
    u[2] = -λ[1]
end

function dae_ϕ(t, x, ϕ, params)
    ϕ[1] = x[2] - x[1]
end

function dae_ψ(t, x, v, ψ, params)
    ψ[1] = v[2] - v[1]
end


function pdae_v(t, q, p, v, params)
    v[1] = p[1]
end

function pdae_f(t, q, p, f, params)
    f[1] = q[1]
end

function pdae_p(t, q, v, p, params)
    p[1] = v[1]
end

function pdae_u(t, q, p, λ, u, params)
    u[1] = λ[1]
end

function pdae_g(t, q, p, λ, g, params)
    g[1] = λ[1]
end

function pdae_ϕ(t, q, p, ϕ, params)
    ϕ[1] = p[1] - q[1]
end

function pdae_ψ(t, q, p, λ, μ, ψ, params)
    ψ[1] = μ[1] - λ[1]
end

function pdae_h(t, q, p)
    p[1]^2/2 + q[1]^2/2
end

hdae_ω = hode_ω

idae_ϑ = iode_ϑ
idae_f = iode_f
idae_v = _idae_default_v̄

ldae_l = lode_l
ldae_ω = lode_ω
ldae_v = _ldae_default_v̄


function sde_v(t, q, v, params)
    @unpack λ = params
    v[1] = λ*q[1]
    v[2] = λ*q[2]
end

function sde_B(t, q, B::AbstractVector, params)
    @unpack μ = params
    B[1] = μ*q[1]
    B[2] = μ*q[2]
end

function sde_B(t, q, B::AbstractMatrix, params)
    @unpack μ = params
    B[1,:] = μ*q[1]
    B[2,:] = μ*q[2]
end


function psde_v(t, q, p, v, params)
    v[1] = p[1]
end

function psde_f(t, q, p, f, params)
    f[1] = - q[1]
end

function psde_B(t, q, p, B, params)
    @unpack noise_intensity = params
    B[1,1] = noise_intensity * p[1]
end

function psde_G(t, q, p, G, params)
    @unpack noise_intensity = params
    G[1,1] = - noise_intensity * q[1]
end


function spsde_v(t, q, p, v, params)
    v[1] = p[1]
end

function spsde_f1(t, q, p, f, params)
    f[1] = - q[1]
end

function spsde_f2(t, q, p, f, params)
    f[1] = 0
end

function spsde_B(t, q, p, B, params)
    @unpack noise_intensity = params
    B[1,1] = noise_intensity * p[1]
end

function spsde_G1(t, q, p, G, params)
    @unpack noise_intensity = params
    G[1,1] = - noise_intensity * q[1]
end

function spsde_G2(t, q, p, G, params)
    G[1,1] = 0
end
