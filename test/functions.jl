

function ode_v(t, x, ẋ)
    ẋ[1] = x[2]
    ẋ[2] = 2x[1]
end


function v_sode_1(t, x, v)
    v[1] = x[2]
end

function v_sode_2(t, x, v)
    v[2] = 2x[1]
end

function q_sode_1(t, x̄, x)
    x[1] = x̄[1]
    x[2] = x̄[2]
end

function q_sode_2(t, x̄, x)
    x[1] = x̄[1]
    x[2] = x̄[2]
end


function pode_v(t, q, p, v)
    v[1] = p[1]
end

function pode_f(t, q, p, f)
    f[1] = 2q[1]
end

function pode_h(t, q, p)
    p[1]^2/2 + cos(q[1])
end

function iode_ϑ(t, q, v, p)
    p[1] = v[1]
end

function iode_f(t, q, v, f)
    f[1] = sin(q[1])
end

function iode_u(t, q, λ, u)
    u[1] = λ[1]
end

function iode_g(t, q, λ, g)
    g[1] = λ[1]
end

function iode_v(t, q, p, v)
    v[1] = p[1]
end

function iode_h(t, q, v)
    v[1]^2/2 + cos(q[1])
end

function lode_l(t, q, v)
    v[1]^2/2 - cos(q[1])
end

function lode_ω(t, q, v, ω)
    ω[1,1] = sin(q[1])
    ω[1,2] = 0
    ω[2,1] = 0
    ω[2,2] = v[1]
end


function dae_v(t, x, v)
    v[1] = x[1]
    v[2] = x[2]
end

function dae_u(t, x, λ, u)
    u[1] = +λ[1]
    u[2] = -λ[1]
end

function dae_ū(t, x, λ, u)
    u[1] = +λ[1]
    u[2] = -λ[1]
end

function dae_ϕ(t, x, λ, ϕ)
    ϕ[1] = x[2] - x[1]
end

function dae_ψ(t, x, v, λ, ψ)
    ψ[1] = v[2] - v[1]
end


function pdae_v(t, q, p, v)
    v[1] = p[1]
end

function pdae_f(t, q, p, f)
    f[1] = q[1]
end

function pdae_p(t, q, v, p)
    p[1] = v[1]
end

function pdae_u(t, q, p, λ, u)
    u[1] = λ[1]
end

function pdae_g(t, q, p, λ, g)
    g[1] = λ[1]
end

function pdae_ϕ(t, q, p, λ, ϕ)
    ϕ[1] = p[1] - q[1]
end

function pdae_ψ(t, q, p, λ, μ, ψ)
    ψ[1] = μ[1] - λ[1]
end

function pdae_h(t, q, p)
    p[1]^2/2 + q[1]^2/2
end


function sde_v(λ, t, q, v)
    v[1] = λ*q[1]
    v[2] = λ*q[2]
end

function sde_B(μ, t, q, u)
    u[1] = μ*q[1]
    u[2] = μ*q[2]
end


function psde_v(t, q, p, v)
    v[1] = p[1]
end

function psde_f(t, q, p, f)
    f[1] = - q[1]
end

function psde_B(t, q, p, B)
    B[1,1] = noise_intensity * p[1]
end

function psde_G(t, q, p, G)
    G[1,1] = - noise_intensity * q[1]
end


function spsde_v(t, q, p, v)
    v[1] = p[1]
end

function spsde_f1(t, q, p, f)
    f[1] = - q[1]
end

function spsde_f2(t, q, p, f)
    f[1] = 0
end

function spsde_B(t, q, p, B)
    B[1,1] = noise_intensity * p[1]
end

function spsde_G1(t, q, p, G)
    G[1,1] = - noise_intensity * q[1]
end

function spsde_G2(t, q, p, G)
    G[1,1] = 0
end

    