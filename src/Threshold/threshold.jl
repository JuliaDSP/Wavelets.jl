using ..Threshold


# THRESHOLD TYPES

abstract type THType end
struct HardTH <: THType end
struct SoftTH <: THType end
struct SemiSoftTH <: THType end
struct SteinTH <: THType end
struct BiggestTH <: THType end
struct PosTH <: THType end
struct NegTH <: THType end



# biggest m-term approximation (best m-term approximation for orthogonal transforms)
# result is m-sparse
function threshold!(x::AbstractArray{<:Number}, TH::BiggestTH, m::Int)
    @assert m >= 0
    n = length(x)
    m > n && (m = n)
    ind = sortperm(x, alg=QuickSort, by=abs)
    @inbounds begin
        for i = 1:n-m
            x[ind[i]] = 0
        end
    end
    return x
end

# hard
function threshold!(x::AbstractArray{<:Number}, TH::HardTH, t::Real)
    @assert t >= 0
    @inbounds for i in eachindex(x)
        if abs(x[i]) <= t
            x[i] = 0
        end
    end
    return x
end

# soft
function threshold!(x::AbstractArray{<:Number}, TH::SoftTH, t::Real)
    @assert t >= 0
    @inbounds for i in eachindex(x)
        sh = abs(x[i]) - t
        if sh < 0
            x[i] = 0
        else
            x[i] = sign(x[i]) * sh
        end
    end
    return x
end

# semisoft
function threshold!(x::AbstractArray{<:Number}, TH::SemiSoftTH, t::Real)
    @assert t >= 0
    @inbounds for i in eachindex(x)
        if x[i] <= 2 * t
            sh = abs(x[i]) - t
            if sh < 0
                x[i] = 0
            elseif sh - t < 0
                x[i] = sign(x[i]) * sh * 2
            end
        end
    end
    return x
end

# stein
function threshold!(x::AbstractArray{<:Number}, TH::SteinTH, t::Real)
    @assert t >= 0
    @inbounds for i in eachindex(x)
        sh = 1 - t * t / (x[i] * x[i])
        if sh < 0
            x[i] = 0
        else
            x[i] = x[i] * sh
        end
    end
    return x
end

# shrink negative elements to 0
function threshold!(x::AbstractArray{<:Number}, TH::NegTH)
    @inbounds for i in eachindex(x)
        if x[i] < 0
            x[i] = 0
        end
    end
    return x
end

# shrink positive elements to 0
function threshold!(x::AbstractArray{<:Number}, TH::PosTH)
    @inbounds for i in eachindex(x)
        if x[i] > 0
            x[i] = 0
        end
    end
    return x
end

# the non inplace functions
function threshold(x::AbstractArray{T}, TH::THType, t::Real) where {T<:Number}
    y = Vector{T}(undef, size(x))
    return threshold!(copyto!(y, x), TH, t)
end


function threshold(x::AbstractArray{T}, TH::THType) where {T<:Number}
    y = Vector{T}(undef, size(x))
    return threshold!(copyto!(y, x), TH)
end
