# MODWT for GPU arrays
# Uses KernelAbstractions kernels for the modwt_step and imodwt_step inner loops

# ============================================================================
# modwt_step kernel — each thread computes one output element
# ============================================================================

# One launch computes both scaling and detail coefficients for a single MODWT level.
@kernel function _modwt_step_kernel!(
        v1, w1, @Const(v), @Const(h), @Const(g),
        N::Int32, L::Int32, j_pow::Int32
    )
    I = @index(Global, Cartesian)
    t = Int32(I[1])
    k = t
    wt = h[1] * v[k]
    vt = g[1] * v[k]
    for n in 2:L
        k -= j_pow
        if k <= 0
            k = mod1(k, N)
        end
        @inbounds wt += h[n] * v[k]
        @inbounds vt += g[n] * v[k]
    end
    @inbounds v1[t] = vt
    @inbounds w1[t] = wt
end

function modwt_step(
        v::AbstractGPUVector{T}, j::Integer, h::AbstractVector{S},
        g::AbstractVector{S}
    ) where {T <: Number, S <: Number}
    N = length(v)
    L = length(h)
    backend = KernelAbstractions.get_backend(v)

    v1 = KernelAbstractions.zeros(backend, T, N)
    w1 = KernelAbstractions.zeros(backend, T, N)

    h_gpu = h isa AbstractGPUArray ? h : to_device(backend, h)
    g_gpu = g isa AbstractGPUArray ? g : to_device(backend, g)

    j_pow = 2^(j - 1)
    _modwt_step!(v1, w1, v, h_gpu, g_gpu, j_pow)
    return v1, w1
end

function _modwt_step!(
        v1::AbstractGPUVector, w1::AbstractGPUVector, v::AbstractGPUVector,
        h::AbstractGPUVector, g::AbstractGPUVector, j_pow::Integer
    )
    N = length(v)
    L = length(h)
    backend = KernelAbstractions.get_backend(v)
    kernel = _modwt_step_kernel!(backend, 256)
    kernel(v1, w1, v, h, g, Int32(N), Int32(L), Int32(j_pow); ndrange = N)
    return nothing
end

# ============================================================================
# imodwt_step kernel
# ============================================================================

# Each thread reconstructs one sample from the wrapped scaling/detail inputs.
@kernel function _imodwt_step_kernel!(
        v0, @Const(v), @Const(w), @Const(h), @Const(g),
        N::Int32, L::Int32, j_pow::Int32
    )
    I = @index(Global, Cartesian)
    t = Int32(I[1])
    k = t
    val = h[1] * w[k] + g[1] * v[k]
    for n in 2:L
        k += j_pow
        if k > N
            k = mod1(k, N)
        end
        @inbounds val += h[n] * w[k] + g[n] * v[k]
    end
    @inbounds v0[t] = val
end

function imodwt_step(
        v::AbstractGPUVector{T}, w::AbstractGPUVector{T}, j::Integer,
        h::AbstractVector{S}, g::AbstractVector{S}
    ) where {T <: Number, S <: Number}
    length(v) == length(w) ||
        throw(DimensionMismatch("Input array sizes must match"))
    length(h) == length(g) ||
        throw(DimensionMismatch("Filter sizes must match"))
    N = length(v)
    L = length(h)
    backend = KernelAbstractions.get_backend(v)

    v0 = KernelAbstractions.zeros(backend, T, N)

    h_gpu = h isa AbstractGPUArray ? h : to_device(backend, h)
    g_gpu = g isa AbstractGPUArray ? g : to_device(backend, g)

    j_pow = 2^(j - 1)
    _imodwt_step!(v0, v, w, h_gpu, g_gpu, j_pow)
    return v0
end

function _imodwt_step!(
        v0::AbstractGPUVector, v::AbstractGPUVector, w::AbstractGPUVector,
        h::AbstractGPUVector, g::AbstractGPUVector, j_pow::Integer
    )
    length(v) == length(w) ||
        throw(DimensionMismatch("Input array sizes must match"))
    length(h) == length(g) ||
        throw(DimensionMismatch("Filter sizes must match"))
    N = length(v)
    L = length(h)
    backend = KernelAbstractions.get_backend(v)
    kernel = _imodwt_step_kernel!(backend, 256)
    kernel(v0, v, w, h, g, Int32(N), Int32(L), Int32(j_pow); ndrange = N)
    return nothing
end

# ============================================================================
# modwt — full MODWT for GPU
# ============================================================================

function modwt(
        x::AbstractGPUVector{T}, wt::OrthoFilter,
        L::Integer = maxmodwttransformlevels(x)
    ) where {T <: Number}
    L <= maxmodwttransformlevels(x) ||
        throw(ArgumentError("Too many transform levels (length(x) < 2^L)"))
    L >= 1 || throw(ArgumentError("L must be >= 1"))
    g, h = WT.makereverseqmfpair(wt)
    g /= sqrt(2)
    h /= sqrt(2)
    N = length(x)
    backend = KernelAbstractions.get_backend(x)

    W = KernelAbstractions.zeros(backend, T, N, L + 1)
    current = copy(x)
    next = similar(current)
    detail = similar(current)

    # Move filters to GPU once
    h_gpu = to_device(backend, h)
    g_gpu = to_device(backend, g)

    for j in 1:L
        # Ping-pong the scaling buffer while writing details directly into the result matrix.
        _modwt_step!(next, detail, current, h_gpu, g_gpu, 2^(j - 1))
        copyto!(@view(W[:, j]), detail)
        current, next = next, current
    end
    copyto!(@view(W[:, L + 1]), current)
    return W
end

# ============================================================================
# imodwt — full inverse MODWT for GPU
# ============================================================================

function imodwt(xw::AbstractGPUMatrix{T}, wt::OrthoFilter) where {T <: Number}
    g, h = WT.makereverseqmfpair(wt)
    g /= sqrt(2)
    h /= sqrt(2)
    N, Lp1 = size(xw)
    L = Lp1 - 1
    backend = KernelAbstractions.get_backend(xw)

    h_gpu = to_device(backend, h)
    g_gpu = to_device(backend, g)

    current = copy(@view(xw[:, L + 1]))
    next = similar(current)
    for j in L:-1:1
        # Inverse levels also ping-pong buffers to avoid per-level allocations.
        _imodwt_step!(next, current, @view(xw[:, j]), h_gpu, g_gpu, 2^(j - 1))
        current, next = next, current
    end
    return current
end
