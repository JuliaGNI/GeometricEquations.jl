
function_v_dummy(t, q, p, v) = nothing


function get_λ₀(q₀::AbstractVector{DT}, λ₀::AbstractVector{DT}) where {DT <: Number}
    zero(λ₀)
end

function get_λ₀(q₀::AbstractVector{DT}, λ₀::AbstractVector{AT}) where {DT <: Number, AT <: AbstractArray{DT}}
    zero(λ₀[begin])
end

function get_λ₀(q₀::AbstractVector{AT}, λ₀::AbstractVector{DT}) where {DT <: Number, AT <: AbstractArray{DT}}
    [zero(λ₀) for i in eachindex(q₀)]
end

function get_λ₀(q₀::AbstractVector{AT}, λ₀::AbstractVector{AT}) where {DT <: Number, AT <: AbstractArray{DT}}
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
