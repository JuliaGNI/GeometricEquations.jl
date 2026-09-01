
@doc raw"""
`WienerProcess`: an ``m``-dimensional standard Wiener process.

The type records only *which* noise drives an equation — its dimension — not any realisation of
it. Increments are drawn by the integrator, which is what knows whether it needs strong
(Gaussian) or weak (discrete) increments, and with which random number generator. Two runs of the
same problem therefore follow different sample paths; use [`GridProcess`](@ref) to prescribe a
path.

The dimension is a type parameter so that an integrator can specialise its caches on it.

### Parameters

* `m`: the number of independent Wiener processes

### Constructors

```julia
WienerProcess(m::Int)
```
"""
struct WienerProcess{m} <: AbstractStochasticProcess
    function WienerProcess(m::Int)
        @assert m > 0 "A Wiener process must have at least one dimension, got m = $m."
        new{m}()
    end
end

@doc raw"""
`GridProcess`: noise increments prescribed on the time grid.

Holds the increments for every step of a run, so the driving path is fixed in advance rather than
drawn as the integration proceeds. That makes a run reproducible, lets two different integrators
be compared on one and the same sample path, and — with `ΔW` and `ΔZ` set to zero — reduces a
stochastic problem to its deterministic drift, which is how a stochastic scheme is checked
against its deterministic counterpart.

Both arrays are indexed `[noise dimension, step]`.

### Fields

* `ΔW`: increments ``\Delta W^r_n`` of the Wiener process, size ``m \times n_t``
* `ΔZ`: the second family of increments, size ``m \times n_t`` — for a strong scheme the time
  integrals ``\Delta Z^r_n``, for a weak one the second discrete random variable

### Constructors

```julia
GridProcess(ΔW, ΔZ)
GridProcess(ΔW)             # ΔZ set to zero
```
"""
struct GridProcess{DT, AT <: AbstractMatrix{DT}} <: AbstractStochasticProcess
    ΔW::AT
    ΔZ::AT

    function GridProcess(ΔW::AT, ΔZ::AT) where {DT, AT <: AbstractMatrix{DT}}
        @assert axes(ΔW)==axes(ΔZ) "ΔW and ΔZ must have the same axes, got $(axes(ΔW)) and $(axes(ΔZ))."
        new{DT, AT}(ΔW, ΔZ)
    end
end

GridProcess(ΔW::AbstractMatrix) = GridProcess(ΔW, zero(ΔW))

GeometricBase.noisedims(::WienerProcess{m}) where {m} = m
GeometricBase.noisedims(process::GridProcess) = size(process.ΔW, 1)

# The stochastic equation types each hold the driving process in a `noise` field and define
# `noise(equation)` themselves; everything else forwards. `AbstractEquationSDE` covers `SDE`,
# `AbstractEquationPSDE` covers both `PSDE` and `SPSDE`, and `GeometricProblem` covers the
# single-problem and ensemble wrappers alike.
const StochasticEquation = Union{AbstractEquationSDE, AbstractEquationPSDE}

GeometricBase.noisedims(equation::StochasticEquation) = noisedims(noise(equation))
GeometricBase.noise(prob::GeometricProblem{<:StochasticEquation}) = noise(equation(prob))
function GeometricBase.noisedims(prob::GeometricProblem{<:StochasticEquation})
    noisedims(noise(prob))
end

"Number of time steps a `GridProcess` prescribes increments for."
nsteps(process::GridProcess) = size(process.ΔW, 2)

function Base.:(==)(p1::GridProcess, p2::GridProcess)
    p1.ΔW == p2.ΔW && p1.ΔZ == p2.ΔZ
end

function Base.show(io::IO, ::WienerProcess{m}) where {m}
    print(io, "Wiener process of dimension ", m)
end

function Base.show(io::IO, process::GridProcess)
    print(io, "Prescribed noise increments of dimension ", noisedims(process),
        " for ", nsteps(process), " time steps")
end
