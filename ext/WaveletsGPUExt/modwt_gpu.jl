# MODWT for GPU arrays
# Uses KernelAbstractions kernels for the modwt_step! and imodwt_step! inner loops

# ============================================================================
# modwt_step! kernel — each thread computes one output element
# ============================================================================

import Wavelets.Transforms: modwt_step!
import Wavelets.Transforms: imodwt_step!

# One launch computes both scaling and detail coefficients for a single MODWT level.
@kernel function _modwt_step_kernel!(
    v1::AbstractVector, w1::AbstractMatrix,
    @Const(v), @Const(h), @Const(g),
    j::Int32, mod_jpow::Int32,
    N::Int32, L::Int32
)
    I = @index(Global)
    t = Int32(I)
    k = t
    wt = h[1] * v[k]
    vt = g[1] * v[k]
    for n in 2:L
        k -= mod_jpow
        if k <= 0
            k += N
        end
        wt += @inbounds h[n] * v[k]
        vt += @inbounds g[n] * v[k]
    end
    @inbounds v1[t] = vt
    @inbounds w1[t, j] = wt
end

function modwt_step!(
    v1::AbstractGPUVector, v::AbstractGPUVector,
    w1::AbstractGPUMatrix,
    j::Int,
    h::AbstractGPUVector, g::AbstractGPUVector,
)
    N = length(v)
    L = length(h)
    mod_jpow = mod(2^(j - 1), N)

    backend = KernelAbstractions.get_backend(v)
    kernel! = _modwt_step_kernel!(backend, 256)

    kernel!(v1, w1, v, h, g, Int32(j), Int32(mod_jpow), Int32(N), Int32(L); ndrange = N)
    return nothing
end

# ============================================================================
# imodwt_step! kernel
# ============================================================================

# Each thread reconstructs one sample from the wrapped scaling/detail inputs.
@kernel function _imodwt_step_kernel!(
    v0::AbstractVector,
    @Const(v), @Const(w),
    @Const(h), @Const(g),
    j::Int32, mod_jpow::Int32,
    N::Int32, L::Int32
)
    I = @index(Global)
    t = Int32(I)
    k = t
    val = h[1] * w[k, j] + g[1] * v[k]
    for n in 2:L
        k += mod_jpow
        if k > N
            k -= N
        end
        val += @inbounds h[n] * w[k, j] + g[n] * v[k]
    end
    @inbounds v0[t] = val
end

function imodwt_step!(
    vn::AbstractGPUVector{T}, v::AbstractGPUVector{T},
    w::AbstractGPUMatrix{T},
    j::Integer,
    h::AbstractGPUVector{T}, g::AbstractGPUVector{T}
) where T<:Number
    N = length(v)
    L = length(h)
    mod_jpow = mod(2^(j - 1), N)

    N == size(w, 1) || throw(DimensionMismatch("Column sizes must match"))
    L == length(g)  || throw(DimensionMismatch("Filter sizes must match"))

    backend = KernelAbstractions.get_backend(v)
    kernel! = _imodwt_step_kernel!(backend, 256)

    kernel!(vn, v, w, h, g, Int32(j), Int32(mod_jpow), Int32(N), Int32(L); ndrange = N)
    return nothing
end

# ============================================================================
# modwt — full MODWT for GPU
# ============================================================================

function modwt(
    x::AbstractGPUVector{T}, wt::OrthoFilter,
    L::Integer=maxmodwttransformlevels(x)
) where {T<:Number}
    L <= maxmodwttransformlevels(x) ||
        throw(ArgumentError("Too many transform levels (length(x) < 2^L)"))
    L >= 1 || throw(ArgumentError("L must be >= 1"))
    g, h = WT.makereverseqmfpair(wt, true, T)
    g ./= sqrt(2)
    h ./= sqrt(2)
    N = length(x)
    backend = KernelAbstractions.get_backend(x)

    W = KernelAbstractions.zeros(backend, T, (N, L + 1))
    v = copy(x)
    next = similar(v)

    # Move filters to GPU once
    h_gpu = to_device(backend, h)
    g_gpu = to_device(backend, g)

    for j in 1:L
        # Ping-pong the scaling buffer while writing details directly into the result matrix.
        modwt_step!(next, v, W, j, h_gpu, g_gpu)
        v, next = next, v
    end
    W[:, end] = v
    return W
end

# ============================================================================
# imodwt — full inverse MODWT for GPU
# ============================================================================

function imodwt(xw::AbstractGPUMatrix{T}, wt::OrthoFilter) where {T <: Number}
    g, h = WT.makereverseqmfpair(wt, true, T)
    g ./= sqrt(2)
    h ./= sqrt(2)
    Lp1 = size(xw, 2)
    L = Lp1 - 1
    backend = KernelAbstractions.get_backend(xw)

    h_gpu = to_device(backend, h)
    g_gpu = to_device(backend, g)

    current = xw[:, Lp1]
    next = similar(current)
    for j in L:-1:1
        # Inverse levels also ping-pong buffers to avoid per-level allocations.
        imodwt_step!(next, current, xw, j, h_gpu, g_gpu)
        current, next = next, current
    end
    return current
end
