# Batched filter-bank DWT/WPT kernels for GPU arrays

# Fused forward pass: approximation and detail reuse the same input traversal.
@kernel function _filtdown_pair_lines_kernel!(
        out, @Const(out_bases), out_stride::Int, out_offset1::Int, out_offset2::Int,
    @Const(x), @Const(x_bases), x_stride::Int,
        nout::Int,
        @Const(f1), shift1::Int, ss1::Bool,
        @Const(f2), shift2::Int, ss2::Bool,
        flen::Int
    )
    idx = @index(Global)
    line = ((idx - 1) ÷ nout) + 1
    k = idx - (line - 1) * nout
    nx = nout << 1
    phase1 = flen - 1 + (ss1 ? 1 : 0)
    phase2 = flen - 1 + (ss2 ? 1 : 0)
    dsidx1 = phase1 + 2 * (k - 1)
    dsidx2 = phase2 + 2 * (k - 1)
    in_base = get_base(x_bases, line)
    out_base = get_base(out_bases, line)

    val1 = zero(eltype(out))
    val2 = zero(eltype(out))
    for j in 0:(flen - 1)
        xidx1 = mod(dsidx1 - j + shift1, nx)
        xidx2 = mod(dsidx2 - j + shift2, nx)
        val1 += @inbounds f1[j + 1] * x[in_base + xidx1 * x_stride]
        val2 += @inbounds f2[j + 1] * x[in_base + xidx2 * x_stride]
    end
    @inbounds out[out_base + out_offset1 + (k - 1) * out_stride] = val1
    @inbounds out[out_base + out_offset2 + (k - 1) * out_stride] = val2
end

# Each thread reconstructs one output sample by traversing the wrapped half-rate input.
@kernel function _filtup_lines_kernel!(
        out, @Const(out_bases), out_stride::Int, out_offset::Int,
    @Const(x), @Const(x_bases), x_stride::Int, x_offset::Int,
        nout::Int, @Const(f), flen::Int,
        shift::Int, ss::Bool, add2out::Bool
    )
    idx = @index(Global)
    line = ((idx - 1) ÷ nout) + 1
    i = idx - (line - 1) * nout
    nx = nout >> 1
    istart = flen - (shift % 2)
    dsshift = ss ? 1 : 0
    shift_half = shift >> 1
    n = istart + (i - 1)
    in_base = get_base(x_bases, line) + x_offset
    out_base = get_base(out_bases, line) + out_offset

    val = zero(eltype(out))
    for j in 0:(flen - 1)
        pos = mod1(n - j, nout)
        if isodd(pos + dsshift)
            xk = mod(((pos - 1) >> 1) + shift_half, nx)
            val += @inbounds f[j + 1] * x[in_base + xk * x_stride]
        end
    end

    if add2out
        @inbounds out[out_base + (i - 1) * out_stride] += val
    else
        @inbounds out[out_base + (i - 1) * out_stride] = val
    end
end

# Fused inverse pass: approximation and detail contributions are accumulated in one launch.
@kernel function _filtup_pair_lines_kernel!(
        out, @Const(out_bases), out_stride::Int, out_offset::Int,
    @Const(x1), @Const(x1_bases), x1_stride::Int, x_offset1::Int,
    @Const(x2), @Const(x2_bases), x2_stride::Int, x_offset2::Int,
        nout::Int,
        @Const(f1), shift1::Int, ss1::Bool,
        @Const(f2), shift2::Int, ss2::Bool,
        flen::Int
    )
    idx = @index(Global)
    line = ((idx - 1) ÷ nout) + 1
    i = idx - (line - 1) * nout
    nx = nout >> 1
    in_base1 = get_base(x1_bases, line)
    in_base2 = get_base(x2_bases, line)
    out_base = get_base(out_bases, line) + out_offset

    istart1 = flen - (shift1 % 2)
    dsshift1 = ss1 ? 1 : 0
    shift_half1 = shift1 >> 1
    n1 = istart1 + (i - 1)

    istart2 = flen - (shift2 % 2)
    dsshift2 = ss2 ? 1 : 0
    shift_half2 = shift2 >> 1
    n2 = istart2 + (i - 1)

    val = zero(eltype(out))
    for j in 0:(flen - 1)
        pos1 = mod1(n1 - j, nout)
        if isodd(pos1 + dsshift1)
            xk1 = mod(((pos1 - 1) >> 1) + shift_half1, nx)
            val += @inbounds f1[j + 1] * x1[in_base1 + x_offset1 + xk1 * x1_stride]
        end

        pos2 = mod1(n2 - j, nout)
        if isodd(pos2 + dsshift2)
            xk2 = mod(((pos2 - 1) >> 1) + shift_half2, nx)
            val += @inbounds f2[j + 1] * x2[in_base2 + x_offset2 + xk2 * x2_stride]
        end
    end

    @inbounds out[out_base + (i - 1) * out_stride] = val
end

function batched_filtdown_pair!(
    out, out_bases, out_stride::Int, out_offset1::Int, out_offset2::Int,
    x, x_bases, x_stride::Int, nout::Int,
    f1, shift1::Int, ss1::Bool,
    f2, shift2::Int, ss2::Bool
)
    nout == 0 && return out
    nlines = length(out_bases)
    nlines == 0 && return out
    length(f1) == length(f2) || throw(DimensionMismatch("Filter sizes must match"))
    kernel = _filtdown_pair_lines_kernel!(KernelAbstractions.get_backend(out))
    kernel(
        out, out_bases, Int(out_stride), Int(out_offset1), Int(out_offset2), x, x_bases, Int(x_stride),
        Int(nout), f1, Int(shift1), ss1, f2, Int(shift2), ss2, Int(length(f1)); ndrange = nlines * nout
    )
    return out
end

function batched_filtup!(
    out, out_bases, out_stride::Int, out_offset::Int,
    x, x_bases, x_stride::Int, x_offset::Int, nout::Int,
    f, shift::Int, ss::Bool, add2out::Bool
)
    nout == 0 && return out
    nlines = length(out_bases)
    nlines == 0 && return out
    kernel = _filtup_lines_kernel!(KernelAbstractions.get_backend(out))
    kernel(
        out, out_bases, Int(out_stride), Int(out_offset), x, x_bases, Int(x_stride), Int(x_offset),
        Int(nout), f, Int(length(f)), Int(shift), ss, add2out; ndrange = nlines * nout
    )
    return out
end

function batched_filtup_pair!(
    out, out_bases, out_stride::Int, out_offset::Int,
    x1, x1_bases, x1_stride::Int, x_offset1::Int,
    x2, x2_bases, x2_stride::Int, x_offset2::Int,
    nout::Int,
    f1, shift1::Int, ss1::Bool,
    f2, shift2::Int, ss2::Bool
)
    nout == 0 && return out
    nlines = length(out_bases)
    nlines == 0 && return out
    length(f1) == length(f2) || throw(DimensionMismatch("Filter sizes must match"))
    kernel = _filtup_pair_lines_kernel!(KernelAbstractions.get_backend(out))
    kernel(
        out, out_bases, Int(out_stride), Int(out_offset),
        x1, x1_bases, Int(x1_stride), Int(x_offset1),
        x2, x2_bases, Int(x2_stride), Int(x_offset2),
        Int(nout), f1, Int(shift1), ss1, f2, Int(shift2), ss2, Int(length(f1)); ndrange = nlines * nout
    )
    return out
end

function _dwt!(
    y::AbstractGPUVector{Ty}, x::AbstractGPUVector{Tx},
    filter::OrthoFilter, L::Integer,
    fw::Bool
) where {Tx<:Number,Ty<:Number}
    T = promote_type(Tx, Ty)
    n = length(x)
    size(x) == size(y) || throw(DimensionMismatch("in and out array size must match"))
    L >= 0 || throw(ArgumentError("L must be non-negative"))
    sufficientpoweroftwo(y, L) || throw(ArgumentError("size must have a sufficient power of 2 factor"))
    y === x && throw(ArgumentError("in array is out array"))
    L == 0 && return copyto!(y, x)

    backend = KernelAbstractions.get_backend(x)
    scfilter_cpu, dcfilter_cpu = WT.makereverseqmfpair(filter, fw, T)
    scfilter = to_device(backend, scfilter_cpu)
    dcfilter = to_device(backend, dcfilter_cpu)
    filtlen = length(filter)
    s = x
    snew = L > 1 ? similar(x, T, n >> 1) : similar(x, T, 0)
    lrange = fw ? (1:1:L) : (L:-1:1)
    bases = LineBases(1)

    for l in lrange
        if fw
            half = detailn(n, l)
            batched_filtdown_pair!(
                y, bases, 1, detailindex(n, l, 1) - 1, 0,
                s, bases, 1, half,
                dcfilter, -filtlen + 1, true,
                scfilter, 0, false
            )
        else
            nout = detailn(n, l - 1)
            half = nout >> 1
            batched_filtup_pair!(
                y, bases, 1, 0,
                s, bases, 1, 0,
                x, bases, 1, detailindex(n, l, 1) - 1,
                nout,
                scfilter, -filtlen + 1, false,
                dcfilter, 0, true
            )
        end
        if l != lrange[end]
            copyto!(snew, 1, y, 1, detailn(n, fw ? l : l - 1))
            s = snew
        end
    end

    return y
end

function _dwt!(
    y::AbstractGPUMatrix{Ty}, x::AbstractGPUMatrix{Tx},
    filter::OrthoFilter, L::Integer,
    fw::Bool
) where {Tx<:Number,Ty<:Number}
    m, n = size(x)
    T = promote_type(Tx, Ty)
    size(x) == size(y) || throw(DimensionMismatch("in and out array size must match"))
    L >= 0 || throw(ArgumentError("L must be non-negative"))
    sufficientpoweroftwo(y, L) || throw(ArgumentError("size must have a sufficient power of 2 factor"))
    y === x && throw(ArgumentError("in array is out array"))
    L == 0 && return copyto!(y, x)

    backend = KernelAbstractions.get_backend(x)
    scfilter_cpu, dcfilter_cpu = WT.makereverseqmfpair(filter, fw, T)
    scfilter = to_device(backend, scfilter_cpu)
    dcfilter = to_device(backend, dcfilter_cpu)
    tmp = similar(y, T, size(y))
    filtlen = length(filter)

    if fw
        msub = m
        nsub = n
    else
        msub = div(m, 2^(L - 1))
        nsub = div(n, 2^(L - 1))
        copyto!(y, x)
    end
    current = x
    for _ in 1:L
        col_bases = LineBases(nsub, nsub, m, 0, 1)
        row_bases = LineBases(msub, msub, 1, 0, 1)
        half_cols = msub >> 1
        half_rows = nsub >> 1
        if fw
            batched_filtdown_pair!(tmp, row_bases, m, half_rows * m, 0, current, row_bases, m, half_rows, dcfilter, -filtlen + 1, true, scfilter, 0, false)
            batched_filtdown_pair!(y, col_bases, 1, half_cols, 0, tmp, col_bases, 1, half_cols, dcfilter, -filtlen + 1, true, scfilter, 0, false)
            msub >>= 1
            nsub >>= 1
        else
            batched_filtup_pair!(tmp, col_bases, 1, 0, current, col_bases, 1, 0, current, col_bases, 1, half_cols, msub, scfilter, -filtlen + 1, false, dcfilter, 0, true)
            batched_filtup_pair!(y, row_bases, m, 0, tmp, row_bases, m, 0, tmp, row_bases, m, half_rows * m, nsub, scfilter, -filtlen + 1, false, dcfilter, 0, true)
            msub <<= 1
            nsub <<= 1
        end
        current = y
    end

    return y
end

function _dwt!(
    y::AbstractGPUArray{Ty,3}, x::AbstractGPUArray{Tx,3},
    filter::OrthoFilter, L::Integer,
    fw::Bool
) where {Tx<:Number,Ty<:Number}
    m, n, d = size(x)
    T = promote_type(Tx, Ty)
    size(x) == size(y) || throw(DimensionMismatch("in and out array size must match"))
    L >= 0 || throw(ArgumentError("L must be non-negative"))
    sufficientpoweroftwo(y, L) || throw(ArgumentError("size must have a sufficient power of 2 factor"))
    y === x && throw(ArgumentError("in array is out array"))
    L == 0 && return copyto!(y, x)

    backend = KernelAbstractions.get_backend(x)
    scfilter_cpu, dcfilter_cpu = WT.makereverseqmfpair(filter, fw, T)
    scfilter = to_device(backend, scfilter_cpu)
    dcfilter = to_device(backend, dcfilter_cpu)
    tmp1 = similar(y, T, size(y))
    tmp2 = similar(y, T, size(y))
    filtlen = length(filter)
    row_stride = m
    plane_stride = m * n

    if fw
        msub, nsub, dsub = m, n, d
    else
        msub = div(m, 2^(L - 1))
        nsub = div(n, 2^(L - 1))
        dsub = div(d, 2^(L - 1))
        copyto!(y, x)
    end
    current = x
    for _ in 1:L
        col_bases = LineBases(nsub * dsub, nsub, m, plane_stride, 1)
        row_bases = LineBases(msub * dsub, msub, 1, plane_stride, 1)
        plane_bases = LineBases(msub * nsub, msub, 1, m, 1)
        half_m = msub >> 1
        half_n = nsub >> 1
        half_d = dsub >> 1
        if fw
            batched_filtdown_pair!(tmp1, plane_bases, plane_stride, half_d * plane_stride, 0, current, plane_bases, plane_stride, half_d, dcfilter, -filtlen + 1, true, scfilter, 0, false)
            batched_filtdown_pair!(tmp2, row_bases, row_stride, half_n * row_stride, 0, tmp1, row_bases, row_stride, half_n, dcfilter, -filtlen + 1, true, scfilter, 0, false)
            batched_filtdown_pair!(y, col_bases, 1, half_m, 0, tmp2, col_bases, 1, half_m, dcfilter, -filtlen + 1, true, scfilter, 0, false)
            msub >>= 1
            nsub >>= 1
            dsub >>= 1
        else
            batched_filtup_pair!(tmp1, col_bases, 1, 0, current, col_bases, 1, 0, current, col_bases, 1, half_m, msub, scfilter, -filtlen + 1, false, dcfilter, 0, true)
            batched_filtup_pair!(tmp2, row_bases, row_stride, 0, tmp1, row_bases, row_stride, 0, tmp1, row_bases, row_stride, half_n * row_stride, nsub, scfilter, -filtlen + 1, false, dcfilter, 0, true)
            batched_filtup_pair!(y, plane_bases, plane_stride, 0, tmp2, plane_bases, plane_stride, 0, tmp2, plane_bases, plane_stride, half_d * plane_stride, dsub, scfilter, -filtlen + 1, false, dcfilter, 0, true)
            msub <<= 1
            nsub <<= 1
            dsub <<= 1
        end
        current = y
    end

    return y
end

function _wpt!(
    y::AbstractGPUVector{Ty}, x::AbstractGPUVector{Tx},
    filter::OrthoFilter, tree::BitVector,
    fw::Bool
) where {Tx<:Number,Ty<:Number}
    T = promote_type(Tx, Ty)
    size(x) == size(y) || throw(DimensionMismatch("in and out array size must match"))
    y === x && throw(ArgumentError("in array is out array"))
    isvalidtree(y, tree) || throw(ArgumentError("invalid tree"))
    !tree[1] && return copyto!(y, x)

    backend = KernelAbstractions.get_backend(x)
    scfilter_cpu, dcfilter_cpu = WT.makereverseqmfpair(filter, fw, T)
    scfilter = to_device(backend, scfilter_cpu)
    dcfilter = to_device(backend, dcfilter_cpu)
    tmp = copyto!(similar(y, T, length(y)), x)
    current = tmp
    output = y
    filtlen = length(filter)
    n = length(x)
    Lmax = maxtransformlevels(n)

    for L in 1:Lmax
        Lfw = fw ? (L - 1) : (Lmax - L)
        seglen = detailn(n, Lfw)
        bases_cpu = segment_bases(n, seglen)
        bases = LineBases(n ÷ seglen, n ÷ seglen, seglen, 0, 1)
        copy_lines!(output, bases, 1, current, bases, 1, seglen)
        treeind = 2^Lfw - 1
        active_idx = tree[(treeind+1):(treeind+length(bases_cpu))]
        if !isempty(active_idx)
            active_bases = to_device(backend, bases_cpu[active_idx])
            half = seglen >> 1
            if fw
                batched_filtdown_pair!(output, active_bases, 1, half, 0, current, active_bases, 1, half, dcfilter, -filtlen + 1, true, scfilter, 0, false)
            else
                batched_filtup!(output, active_bases, 1, 0, current, active_bases, 1, 0, seglen, scfilter, -filtlen + 1, false, false)
                batched_filtup!(output, active_bases, 1, 0, current, active_bases, 1, half, seglen, dcfilter, 0, true, true)
            end
        end
        current, output = output, current
    end

    current === y || copyto!(y, current)
    return y
end
