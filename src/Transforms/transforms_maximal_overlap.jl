
"""
Perform one level of the maximal overlap discrete wavelet transform (MODWT) on
the the vector `v`, which is either the `j`th level scaling coefficients or,
if `j==1`, the raw signal. The vectors `h` and `g` are the MODWT detail and
scaling filters.

Returns a tuple `(v, w)` of the scaling and detail coefficients at level `j+1`.
"""
function modwt_step!(
    vtmp::Vector{T}, v::Vector{T},
    W::Matrix{T},
    j::Integer,
    h::Vector{T}, g::Vector{T}
) where {T<:Number}
    N = length(v)
    L = length(h)
    mod_delta = mod(2^(j - 1), N)
    for t in 1:N
        k = t
        W[t, j] = h[1] * v[k]
        vtmp[t] = g[1] * v[k]
        for n in 2:L
            k -= mod_delta
            if k <= 0
                k += N
            end
            W[t, j] += h[n] * v[k]
            vtmp[t] += g[n] * v[k]
        end
    end
    return nothing
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
function modwt(
    x::AbstractVector{T}, wt::OrthoFilter,
    L::Integer=maxmodwttransformlevels(x)
) where T<:Number
    L <= maxmodwttransformlevels(x) ||
        throw(ArgumentError("Too many transform levels (length(x) < 2^L)"))
    L >= 1 || throw(ArgumentError("L must be >= 1"))
    g, h = WT.makereverseqmfpair(wt, true, T)
    g ./= sqrt(2)
    h ./= sqrt(2)
    N = length(x)
    W = similar(x, T, (N, L + 1))
    V = copyto!(similar(x, T), x)
    vtmp = similar(V)
    for j in 1:L
        modwt_step!(vtmp, V, W, j, h, g)
        vtmp, V = V, vtmp
    end
    W[:, end] = V
    return W
end


"""
Perform one level of the inverse maximal overlap discrete wavelet transform
(MODWT). Takes the `j`th-level scaling and detail coefficients `v` and `w`
and returns a vector of the `j-1`th level scaling coefficients. The vectors
`h` and `g` are the MODWT detail and scaling filters.
"""
function imodwt_step!(
    vn::Vector{T}, v::Vector{T},
    w::Matrix{T},
    j::Integer,
    h::Vector{T}, g::Vector{T}
) where {T<:Number}
    N = length(v)
    L = length(h)
    N == size(w, 1) ||
        throw(DimensionMismatch("Column sizes must match"))
    L == length(g) ||
        throw(DimensionMismatch("Filter sizes must match"))
    mod_delta = mod(2^(j - 1), N)
    for t in 1:N
        k = t
        vn[t] = h[1] * w[k, j] + g[1] * v[k]
        for n in 2:L
            k += mod_delta
            if k > N
                k -= N
            end
            vn[t] += h[n] * w[k, j] + g[n] * v[k]
        end
    end
end

"""
Perform an inverse maximal overlap discrete wavelet transform (MODWT) of `xw`,
the inverse of `modwt(x, wt, L)`.
"""
function imodwt(xw::Matrix{T}, wt::OrthoFilter) where T<:Number
    g, h = WT.makereverseqmfpair(wt, true, T)
    g ./= sqrt(2)
    h ./= sqrt(2)
    L = size(xw, 2)
    x = xw[:, L]
    xn = similar(x)
    for j in L-1:-1:1
        imodwt_step!(xn, x, xw, j, h, g)
        x, xn = xn, x
    end
    return x
end

#####################################
### Only for benchmarks and tests ###
#####################################

function modwt_step(
    v::AbstractVector{T},
    j::Integer,
    h::Vector{T}, g::Vector{T}
) where {T<:Number}
    N = length(v)
    L = length(h)

    v1 = similar(v)
    w1 = fill!(similar(v, (N, j)), zero(T))

    h_new = copyto!(similar(v, L), h)
    g_new = copyto!(similar(v, L), g)

    modwt_step!(v1, v, w1, j, h_new, g_new)
    return v1, w1
end

function imodwt_step(
    v::AbstractVector{T}, w::AbstractMatrix{T},
    j::Integer,
    h::Vector{T}, g::Vector{T}
) where {T<:Number}
    N = length(v)
    L = length(h)
    N == size(w, 1) || throw(DimensionMismatch("Column sizes must match"))
    L == length(g)  || throw(DimensionMismatch("Filter sizes must match"))

    vn = zero(v)

    h_new = copyto!(similar(v, L), h)
    g_new = copyto!(similar(v, L), g)

    imodwt_step!(vn, v, w, j, h_new, g_new)
    return vn
end
