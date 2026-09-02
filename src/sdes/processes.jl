
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
WienerProcess{m}()
```
"""
struct WienerProcess{m} <: AbstractStochasticProcess
    function WienerProcess{m}() where {m}
        @assert m isa Int "The dimension of a Wiener process must be an `Int`, got m = $m of type $(typeof(m))."
        @assert m > 0 "A Wiener process must have at least one dimension, got m = $m."
        new{m}()
    end
end

WienerProcess(m::Int) = WienerProcess{m}()

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
struct GridProcess{DT, WT <: AbstractMatrix{DT}, ZT <: AbstractMatrix{DT}} <:
       AbstractStochasticProcess
    # The two arrays carry separate type parameters, tied only by their element type, because
    # they routinely arrive as different array types even when they describe the same path:
    # `zero` of a view or of an adjoint is a plain `Matrix`, so requiring one shared type would
    # reject the one-argument constructor's own output for any ΔW that is not already a `Matrix`.
    ΔW::WT
    ΔZ::ZT

    function GridProcess(ΔW::AbstractMatrix{DT}, ΔZ::AbstractMatrix{DT}) where {DT}
        @assert axes(ΔW)==axes(ΔZ) "ΔW and ΔZ must have the same axes, got $(axes(ΔW)) and $(axes(ΔZ))."
        new{DT, typeof(ΔW), typeof(ΔZ)}(ΔW, ΔZ)
    end
end

# Mismatched element types are promoted rather than rejected; the inner constructor above is the
# more specific method and takes over once both arrays agree.
function GridProcess(ΔW::AbstractMatrix, ΔZ::AbstractMatrix)
    DT = promote_type(eltype(ΔW), eltype(ΔZ))
    GridProcess(convert(AbstractMatrix{DT}, ΔW), convert(AbstractMatrix{DT}, ΔZ))
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

"""
    ntime(process::GridProcess)

The number of time steps a `GridProcess` prescribes increments for.

This is `ntime` rather than `nsteps` for the same reason it is on a problem: in this package
`nsteps` counts the substeps of an [`SODE`](@ref) splitting, and the number of steps along the
time grid is `ntime`.
"""
GeometricBase.ntime(process::GridProcess) = size(process.ΔW, 2)

"""
    check_noise(equation, timespan, timestep)

Verify that the noise driving `equation` can supply the whole run, and throw if it cannot.

Only a [`GridProcess`](@ref) can come up short: it carries a fixed number of increments, and an
integrator handed too few would read past the end of `ΔW` partway through the run. A
[`WienerProcess`](@ref) draws its increments as it goes and so is never short.
"""
function check_noise(equation::StochasticEquation, timespan, timestep)
    check_noise(noise(equation), timespan, timestep)
end

check_noise(::AbstractStochasticProcess, timespan, timestep) = true

function check_noise(process::GridProcess, timespan, timestep)
    nt = ntimesteps(timespan, timestep)
    ntime(process) ≥ nt || throw(ArgumentError(
        "GridProcess prescribes increments for $(ntime(process)) time steps, " *
        "but the problem takes $nt steps over $timespan with a time step of $timestep."))
    true
end

function Base.:(==)(p1::GridProcess, p2::GridProcess)
    p1.ΔW == p2.ΔW && p1.ΔZ == p2.ΔZ
end

Base.hash(process::GridProcess, h::UInt) = hash(process.ΔW, hash(process.ΔZ, h))

function Base.show(io::IO, ::WienerProcess{m}) where {m}
    print(io, "Wiener process of dimension ", m)
end

function Base.show(io::IO, process::GridProcess)
    print(io, "Prescribed noise increments of dimension ", noisedims(process),
        " for ", ntime(process), " time steps")
end
