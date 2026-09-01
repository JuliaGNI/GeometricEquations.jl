```@meta
CurrentModule = GeometricEquations
```

# Stochastic Equations

## Theory

### The Stratonovich formulation

The stochastic equation types in this package describe systems driven by an ``m``-dimensional
Wiener process ``W = (W^1, \ldots, W^m)`` in the **Stratonovich** sense, written with the symbol
``\circ``:

```math
dq(t) = v(t, q(t)) \, dt + \sum_{r=1}^{m} B^{\cdot r}(t, q(t)) \circ dW^r(t) ,
\qquad q(t_0) = q_0 .
```

``v`` is the *drift* vector field and ``B`` the ``d \times m`` *diffusion matrix*, one column per
independent Wiener process.

The choice of Stratonovich over Itô is not cosmetic. The Stratonovich differential obeys the
ordinary chain rule, which is what makes it possible to speak of a stochastic Hamiltonian system
at all: a coordinate change, a symplectic form, or a conserved quantity transforms as it does in
the deterministic case. In the Itô calculus each of those acquires a second-order correction term,
and the geometric structure the integrators in
[StochasticIntegrators.jl](https://github.com/JuliaGNI/StochasticIntegrators.jl) preserve would
not be there to preserve. The two formulations describe the same processes and are related by a
shift of the drift,

```math
v^{\text{Itô}} = v^{\text{Strat}} + \tfrac{1}{2} \sum_{r} \frac{\partial B^{\cdot r}}{\partial q} B^{\cdot r} ,
```

so a model written in Itô form must be converted before it is handed to an [`SDE`](@ref).

### Partitioned and split systems

A stochastic Hamiltonian system carries a Hamiltonian ``H`` and one *stochastic Hamiltonian*
``h_r`` per noise dimension:

```math
\begin{aligned}
dq &= \frac{\partial H}{\partial p} \, dt
    + \sum_{r} \frac{\partial h_r}{\partial p} \circ dW^r , \\
dp &= -\frac{\partial H}{\partial q} \, dt
    - \sum_{r} \frac{\partial h_r}{\partial q} \circ dW^r .
\end{aligned}
```

This is a [`PSDE`](@ref), with ``v = \partial H / \partial p``, ``f = -\partial H / \partial q``,
``B = \partial h_r / \partial p`` and ``G = -\partial h_r / \partial q``. Keeping ``q`` and ``p``
separate lets an integrator treat them differently, which is what a symplectic partitioned
Runge-Kutta method needs.

Adding dissipation or an external force gives a *forced* Hamiltonian system,

```math
\begin{aligned}
dq &= \frac{\partial H}{\partial p} \, dt
    + \sum_{r} \frac{\partial h_r}{\partial p} \circ dW^r , \\
dp &= \left[ -\frac{\partial H}{\partial q} + F(q,p) \right] dt
    + \sum_{r} \left[ -\frac{\partial h_r}{\partial q} + f_r(q,p) \right] \circ dW^r ,
\end{aligned}
```

and this is an [`SPSDE`](@ref) — a *split* partitioned system, in which the ``p`` equation is
presented as two separate pieces rather than one sum: ``f_1 = -\partial H / \partial q`` against
``f_2 = F``, and ``G_1 = -\partial h_r / \partial q`` against ``G_2 = f_r``.

The split is not bookkeeping. A Lagrange-d'Alembert variational integrator applies **different**
Runge-Kutta coefficients to the Hamiltonian and the forcing terms, and that extra freedom is
exactly what allows the symplecticity conditions to be satisfied in the presence of forcing. An
integrator handed the sum ``f_1 + f_2`` could not do this. Conversely, a system with no forcing is
naturally a [`PSDE`](@ref); writing it as an [`SPSDE`](@ref) with ``f_2 \equiv 0`` is legal but
tells the integrator nothing it can use.

## Implementation

### Equation and problem types

Each of the three equation types has a matching problem type — an
[`EquationProblem`](@ref) specialised to it — and an ensemble type:

| system | equation | problem | ensemble |
|:--|:--|:--|:--|
| ``dq = v \, dt + B \circ dW`` | [`SDE`](@ref) | `SDEProblem` | `SDEEnsemble` |
| partitioned | [`PSDE`](@ref) | `PSDEProblem` | `PSDEEnsemble` |
| split partitioned | [`SPSDE`](@ref) | `SPSDEProblem` | `SPSDEEnsemble` |

An *equation* holds the vector fields and nothing else; a *problem* adds the time span, the time
step, the initial conditions and the parameters. An *ensemble* is a problem with several initial
conditions or several parameter sets.

### Function signatures

Every function writes its result into its first argument and takes the parameters last. For an
[`SDE`](@ref) in ``d`` dimensions driven by ``m`` Wiener processes:

```julia
function v(v, t, q, params)      # drift, length d
    v[1] = ...
end

function B(B, t, q, params)      # diffusion, size d × m
    B[1,1] = ...
end
```

For a [`PSDE`](@ref) the four functions `v`, `f`, `B`, `G` take both `q` and `p`:

```julia
function v(v, t, q, p, params) end
function f(f, t, q, p, params) end
function B(B, t, q, p, params) end
function G(G, t, q, p, params) end
```

and an [`SPSDE`](@ref) takes the six functions `v`, `f1`, `f2`, `B`, `G1`, `G2` with the same
signature.

### The driving noise

Every stochastic equation carries a noise object in its `noise` field, a subtype of
`GeometricBase.AbstractStochasticProcess`. This says *which* process drives the equation. It does
not hold a realisation: increments are drawn by the integrator, because only the integrator knows
whether the scheme it implements needs strong or weak increments.

```@docs
WienerProcess
GridProcess
```

The dimension of the noise lives with the process, and hence with the problem. Two accessors,
declared in `GeometricBase` and extended here, reach it:

| function | returns |
|:--|:--|
| `noise(equation)`, `noise(problem)` | the driving process |
| `noisedims(process)`, `noisedims(equation)`, `noisedims(problem)` | the number of independent Wiener processes |

`noisedims(problem)` is how an integrator sizes its increment vectors, and it is what fixes the
number of columns of `B` and `G`.

Use a [`GridProcess`](@ref) whenever the driving path must be fixed rather than drawn:

* to make a run reproducible;
* to drive two different integrators along one and the same sample path, which is what a
  comparison of methods requires;
* to measure strong convergence order against a reference solution evaluated on the same path;
* with zero increments, to reduce a stochastic problem to its deterministic drift — a sharp check
  that a stochastic method reproduces the deterministic method it is built from.

## Usage

A one-dimensional Kubo oscillator — a harmonic oscillator whose frequency is perturbed by noise —
as an [`SDE`](@ref):

```julia
using GeometricEquations

function kubo_v(v, t, q, params)
    v[1] =  q[2]
    v[2] = -q[1]
end

function kubo_B(B, t, q, params)
    B[1,1] =  params.ν * q[2]
    B[2,1] = -params.ν * q[1]
end

prob = SDEProblem(kubo_v, kubo_B, WienerProcess(1),
                  (0.0, 1.0), 0.01, [0.5, 0.0]; parameters = (ν = 0.1,))

noisedims(prob)   # 1
```

The same system in partitioned form, which is what a symplectic method needs:

```julia
kubo_psde_v(v, t, q, p, params) = (v[1] =  p[1])
kubo_psde_f(f, t, q, p, params) = (f[1] = -q[1])
kubo_psde_B(B, t, q, p, params) = (B[1,1] =  params.ν * p[1])
kubo_psde_G(G, t, q, p, params) = (G[1,1] = -params.ν * q[1])

pprob = PSDEProblem(kubo_psde_v, kubo_psde_f, kubo_psde_B, kubo_psde_G,
                    WienerProcess(1), (0.0, 1.0), 0.01, [0.5], [0.0];
                    parameters = (ν = 0.1,))
```

Prescribing the increments instead of drawing them:

```julia
nt = 100
Δt = 0.01
ΔW = randn(1, nt) .* sqrt(Δt)          # one sample path, fixed in advance
fixed = SDEProblem(kubo_v, kubo_B, GridProcess(ΔW),
                   (0.0, Δt*nt), Δt, [0.5, 0.0]; parameters = (ν = 0.1,))
```

and, with zero increments, the deterministic limit of the same problem:

```julia
deterministic = SDEProblem(kubo_v, kubo_B, GridProcess(zeros(1, nt)),
                           (0.0, Δt*nt), Δt, [0.5, 0.0]; parameters = (ν = 0.1,))
```

Ready-made stochastic problems live in
[GeometricProblems.jl](https://github.com/JuliaGNI/GeometricProblems.jl) — see its
`KuboOscillator` module, which provides all three formulations, their ensembles, a damped variant
and closed-form solutions. Integrators for them are in
[StochasticIntegrators.jl](https://github.com/JuliaGNI/StochasticIntegrators.jl).

!!! note "Ensemble constructors"
    Unlike `ODEEnsemble` and its relatives, `SDEEnsemble`, `PSDEEnsemble` and `SPSDEEnsemble` are
    currently type aliases without a convenience constructor. An ensemble of stochastic problems
    has to be assembled by building the equation and then calling `EnsembleProblem` directly.
