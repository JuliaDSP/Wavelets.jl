using ..Util
using ..WT
using ..Transforms
using LinearAlgebra: norm

# PLOTTING UTILITIES


"""
Return levels and detail coefficient centers on the interval [0,r) above (>=) threshold t
as tuple (d,l).
Parameters:
- x: AbstractVector - input vector

- t: Real - threshold value

- r: Real - range value
"""
function wplotdots(x::AbstractVector, t::Real=0, r::Real=1)
    if !isdyadic(x)
        throw(ArgumentError("array must be of dyadic size"))
    end
    n = length(x)                       # length of input vector
    c = wcount(x, t, level=0)           # number of detail coefficients above threshold
    d = Vector{Float64}(undef, c)       # detail coefficient centers
    l = Vector{Int}(undef, c)           # detail coefficient levels
    p = LinRange(0, r * (n - 1) / n, n) # positions
    J = ndyadicscales(n)                # number of scales
    k = 1

    @inbounds for j in 0:J-1
        Δ = 1 << (J - j)                 # Δ=2^(J-j) : step size at level j
        start = Δ >> 1                   # start= Δ÷2 = 2^(J-j-1)
        rind = start:Δ:n                 # indices for detail coefficients
        for i in 1:dyadicdetailn(j)
            if x[dyadicdetailindex(j, i)] |> abs ≥ t
                d[k] = p[rind[i]]
                l[k] = j
                k += 1
            end
        end
    end
    return (d, l)
end

"""
Return a matrix of detail coefficient values where row j+1 is level j.
"""
function wplotim(x::AbstractVector)
    if !isdyadic(x)
        throw(ArgumentError("vector must be of dyadic size"))
    end
    n = length(x)
    J = ndyadicscales(n)
    A = zeros(Float64, J, n)

    @inbounds for j in 0:J-1
        dr = dyadicdetailrange(j)      # indices in x for level j details
        m = 1 << (J - j)
        for i in eachindex(dr)
            A[j+1, 1+(i-1)*m:i*m] .= x[dr[i]]
        end
    end

    return A
end

"""
Return an array of scaled detail coefficients and unscaled scaling coefficients
ready to be plotted as an image.
Parameters:
    - x::AbstractArray: The input array.

    - L::Integer: The number of levels.

    - wt::Union{DiscreteWavelet,Nothing}=nothing: The wavelet type.

    - wabs::Bool=true: Whether to take the absolute value of the coefficients.

    - power::Real=0.7: The power to raise the coefficients to.

    - pnorm::Real=1.0: The p-norm to use.
"""
function wplotim(
    x::AbstractArray,
    L::Integer,
    wt::Union{DiscreteWavelet,Nothing}=nothing;
    wabs::Bool=true,
    power::Real=0.7,
    pnorm::Real=1.0
)

    dim = ndims(x)
    n = size(x, 1)
    cn = size(x, 3)  # color dimension

    # We perform various checks :
    # - Array must be of dyadic size
    # - Dimension must be 2 or 3
    # - Array must be square
    # - Third dimension (color dimension) must be 1 or 3

    if !isdyadic(x)
        throw(ArgumentError("array must be of dyadic size"))
    end

    if dim ∉ (2, 3)
        throw(ArgumentError("dimension $(dim) not supported, only dimensions 2 and 3 are supported"))
    end

    if size(x, 2) != n
        throw(ArgumentError("array must be square"))
    end

    if cn ∉ (1, 3)
        throw(ArgumentError("third dimension $(cn) not supported, expected 1 or 3"))
    end

    J = ndyadicscales(n)
    nsc = 2^(J - L)

    # perform wavelet transform
    if !isnothing(wt)
        x = size(x, 3) > 1 ? dwtc(x, wt, L) : dwt(x, wt, L)
    end

    # scaling coefficients
    scs = x[1:nsc, 1:nsc, :]
    scale01!(scs)

    # detail coefficients
    xts = wabs ? abs.(copy(x)) : copy(x)
    xts[1:nsc, 1:nsc, :] .= 0
    scale01!(xts)

    for i ∈ 1:n, j ∈ 1:n
        @inbounds xts[i, j, :] .= norm(xts[i, j, :], pnorm) .^ (power)
    end

    # merge and reshape the final image
    scale01!(xts)
    xts[1:nsc, 1:nsc, :] = scs
    return xts
end


# scale elements of z to the interval [0,1]
function scale01!(z)
    # we have to perform an extra check if z is uniform
    # and in that case we return the unscaled z
    mi, ma = extrema(z)
    ma - mi ≈ 0 ? z : (z .- mi) ./ (ma - mi)
end
