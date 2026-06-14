
"""
Perform one level of the maximal overlap discrete wavelet transform (MODWT) on
the the vector `v`, which is either the `j`th level scaling coefficients or,
if `j==1`, the raw signal.  The vectors `h` and `g` are the MODWT detail and
scaling filters.

Returns a tuple `(v, w)` of the scaling and detail coefficients at level `j+1`.
"""
function modwt_step(v::AbstractVector{T}, j::Integer, h::Array{S,1},
        g::Array{S,1}) where {T <: Number, S <: Number}
    N = length(v)
    L = length(h)
    v1 = zeros(T, N)
    w1 = zeros(T, N)
    for t in 1:N
        k = t
        w1[t] = h[1] * v[k]
        v1[t] = g[1] * v[k]
        for n in 2:L
            k -= 2^(j-1)
            if k <= 0
                k = mod1(k, N)
            end
            w1[t] += h[n] * v[k]
            v1[t] += g[n] * v[k]
        end
    end
    return v1, w1
end


"""
Perform a maximal overlap discrete wavelet transform (MODWT) of the array
`x`.  The wavelet type `wt` determines the transform type and the wavelet
class, see `wavelet`. (Note that the wavelet filter coefficients are scaled
by 1/√2 for the MODWT so that the transform maintains unit energy).

The number of transform levels `L` can be any number <= `maxmodwttransformlevels(x)`
(the default value).

Returns an `n × L+1` matrix (where `n` is the length of `x`) with the wavelet
coefficients for level j in column j.  The scaling coefficients are in the
last (L+1th) column.
"""
function modwt(x::AbstractVector{T}, wt::OrthoFilter,
        L::Integer=maxmodwttransformlevels(x)) where T <: Number
    L <= maxmodwttransformlevels(x) ||
        throw(ArgumentError("Too many transform levels (length(x) < 2^L)"))
    L >= 1 || throw(ArgumentError("L must be >= 1"))
    g, h = WT.makereverseqmfpair(wt)
    g /= sqrt(2)
    h /= sqrt(2)
    N = length(x)
    W = zeros(T, N, L)
    V = deepcopy(x)
    for j in 1:L
        V[:], W[:, j] = modwt_step(V, j, h, g)
    end
    return [W V]
end


"""
Perform one level of the inverse maximal overlap discrete wavelet transform
(MODWT). Takes the `j`th-level scaling and detail coefficients `v` and `w`
and returns a vector of the `j-1`th level scaling coefficients. The vectors
`h` and `g` are the MODWT detail and scaling filters.
"""
function imodwt_step(v::AbstractVector{T}, w::AbstractVector{T}, j::Integer,
        h::Array{S,1}, g::Array{S,1}) where {T<:Number, S<:Number}
    length(v) == length(w) ||
        throw(DimensionMismatch("Input array sizes must match"))
    length(h) == length(g) ||
        throw(DimensionMismatch("Filter sizes must match"))
    N = length(v)
    L = length(h)
    v0 = zeros(T, N)
    for t in 1:N
        k = t
        v0[t] = h[1] * w[k] + g[1] * v[k]
        for n in 2:L
            k += 2^(j-1)
            if k > N
                k = mod1(k, N)
            end
            v0[t] += h[n] * w[k] + g[n] * v[k]
        end
    end
    return v0
end

"""
Perform an inverse maximal overlap discrete wavelet transform (MODWT) of `xw`,
the inverse of `modwt(x, wt, L)`.
"""
function imodwt(xw::Array{T, 2}, wt::OrthoFilter) where T <: Number
    g, h = WT.makereverseqmfpair(wt)
    g /= sqrt(2)
    h /= sqrt(2)
    N, L = size(xw)
    x = deepcopy(xw[:, L])
    for j in L-1:-1:1
        x[:] = imodwt_step(x, xw[:, j], j, h, g)
    end
    return x
end
