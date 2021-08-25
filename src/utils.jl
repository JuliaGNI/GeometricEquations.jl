
function_dummy_v(t, q, v) = nothing
function_dummy_v(t, q, p, v) = nothing


"""
    promote_tspan(tspan)
Convert the `tspan` field of a `GeometricProblem` to a `(tmin, tmax)` tuple, where both
elements are of the same type.
"""
promote_tspan((t1,t2)::Tuple{T,S}) where {T,S} = promote(t1, t2)
promote_tspan(tspan::Number) = (zero(tspan),tspan)
promote_tspan(tspan::Nothing) = (nothing,nothing)
promote_tspan(tspan::AbstractArray) = length(tspan) == 2 ? promote_tspan((first(tspan),last(tspan))) : throw(error("The length of tspan must be two (and preferably, tspan should be a tuple, i.e. (0.0,1.0))."))

function promote_tspan_and_tstep((t1,t2)::Tuple, Δt::Number)
    t = promote(t1, t2, Δt)
    return (t[1], t[2]), t[3]
end


function initial_multiplier(q₀::AbstractVector{DT}, λ₀::AbstractVector{DT}) where {DT <: Number}
    zero(λ₀)
end

function initial_multiplier(q₀::AbstractVector{DT}, λ₀::AbstractVector{AT}) where {DT <: Number, AT <: AbstractArray{DT}}
    zero(λ₀[begin])
end

function initial_multiplier(q₀::AbstractVector{AT}, λ₀::AbstractVector{DT}) where {DT <: Number, AT <: AbstractArray{DT}}
    [zero(λ₀) for i in eachindex(q₀)]
end

function initial_multiplier(q₀::AbstractVector{AT}, λ₀::AbstractVector{AT}) where {DT <: Number, AT <: AbstractArray{DT}}
    [zero(λ₀[begin]) for i in eachindex(q₀)]
end


function symplectic_matrix(::Type{DT}, d::Int) where {DT}
    ω = zeros(DT, 2d, 2d)

    for i in 1:d
        ω[i, d+i] = +1
        ω[d+i, i] = -1
    end

    return ω
end

function symplectic_matrix(::Number, q::AbstractVector{DT}, p::AbstractVector{DT}) where {DT}
    @assert length(q) == length(p)
    symplectic_matrix(DT, length(q))
end
