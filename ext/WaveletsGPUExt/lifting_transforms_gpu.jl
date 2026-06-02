struct GPULiftStep{A}
    coef::A
    shift::Int
    predict::Bool
end

function gpu_makescheme(backend, ::Type{T}, scheme::GLS, fw::Bool) where {T <: Number}
    nsteps = length(scheme.step)
    steps = Vector{GPULiftStep}(undef, nsteps)
    for i in 1:nsteps
        j = fw ? i : (nsteps + 1 - i)
        step = scheme.step[j]
        coef = convert(Vector{T}, step.param.coef .* (fw ? -1 : 1))
        steps[i] = GPULiftStep(to_device(backend, coef), step.param.shift, step.steptype isa WT.PredictStep)
    end
    norm1 = convert(T, fw ? scheme.norm1 : inv(scheme.norm1))
    norm2 = convert(T, fw ? scheme.norm2 : inv(scheme.norm2))
    return steps, norm1, norm2
end

@kernel function _split_lines_kernel!(
        out, out_bases, out_stride::Int,
        @Const(x), x_bases, x_stride::Int, ns::Int
    )
    idx = @index(Global)
    line = ((idx - 1) ÷ ns) + 1
    pos = idx - (line - 1) * ns
    half = ns >> 1
    srcpos = pos <= half ? (2 * pos - 1) : (2 * (pos - half))
    out_idx = get_base(out_bases, line) + (pos - 1) * out_stride
    x_idx = get_base(x_bases, line) + (srcpos - 1) * x_stride
    @inbounds out[out_idx] = x[x_idx]
end

@kernel function _merge_lines_kernel!(
        out, out_bases, out_stride::Int,
        @Const(x), x_bases, x_stride::Int, ns::Int
    )
    idx = @index(Global)
    line = ((idx - 1) ÷ ns) + 1
    pos = idx - (line - 1) * ns
    half = ns >> 1
    srcpos = isodd(pos) ? ((pos + 1) >> 1) : (half + (pos >> 1))
    out_idx = get_base(out_bases, line) + (pos - 1) * out_stride
    x_idx = get_base(x_bases, line) + (srcpos - 1) * x_stride
    @inbounds out[out_idx] = x[x_idx]
end

@kernel function _normalize_lines_kernel!(x, bases, stride::Int, ns::Int, norm1, norm2)
    idx = @index(Global)
    line = ((idx - 1) ÷ ns) + 1
    pos = idx - (line - 1) * ns
    half = ns >> 1
    factor = pos <= half ? norm1 : norm2
    @inbounds x[get_base(bases, line) + (pos - 1) * stride] *= factor
end

@kernel function _normalize_copy_lines_kernel!(
        out, out_bases, out_stride::Int,
        @Const(x), x_bases, x_stride::Int,
        ns::Int, norm1, norm2
    )
    idx = @index(Global)
    line = ((idx - 1) ÷ ns) + 1
    pos = idx - (line - 1) * ns
    half = ns >> 1
    factor = pos <= half ? norm1 : norm2
    out_idx = get_base(out_bases, line) + (pos - 1) * out_stride
    x_idx = get_base(x_bases, line) + (pos - 1) * x_stride
    @inbounds out[out_idx] = factor * x[x_idx]
end

@kernel function _lift_step_lines_kernel!(
        x, bases, stride::Int,
        half::Int, @Const(coef), nc::Int,
        shift::Int, predict::Bool
    )
    idx = @index(Global)
    line = ((idx - 1) ÷ half) + 1
    i = idx - (line - 1) * half
    base = get_base(bases, line)

    if predict
        dst = base + (i - 1) * stride
        @inbounds val = x[dst]
        for k in 1:nc
            srcpos = half + mod1(i + k - 1 - shift, half)
            @inbounds val += coef[k] * x[base + (srcpos - 1) * stride]
        end
        @inbounds x[dst] = val
    else
        dstpos = half + i
        dst = base + (dstpos - 1) * stride
        @inbounds val = x[dst]
        for k in 1:nc
            srcpos = mod1(i + k - 1 - shift, half)
            @inbounds val += coef[k] * x[base + (srcpos - 1) * stride]
        end
        @inbounds x[dst] = val
    end
end

function split_lines!(out, out_bases, out_stride::Int, x, x_bases, x_stride::Int, ns::Int)
    ns == 0 && return out
    nlines = length(out_bases)
    nlines == 0 && return out
    kernel = _split_lines_kernel!(KernelAbstractions.get_backend(out), 256)
    kernel(out, out_bases, Int(out_stride), x, x_bases, Int(x_stride), Int(ns); ndrange = nlines * ns)
    return out
end

function merge_lines!(out, out_bases, out_stride::Int, x, x_bases, x_stride::Int, ns::Int)
    ns == 0 && return out
    nlines = length(out_bases)
    nlines == 0 && return out
    kernel = _merge_lines_kernel!(KernelAbstractions.get_backend(out), 256)
    kernel(out, out_bases, Int(out_stride), x, x_bases, Int(x_stride), Int(ns); ndrange = nlines * ns)
    return out
end

function normalize_lines!(x, bases, stride::Int, ns::Int, norm1, norm2)
    ns == 0 && return x
    nlines = length(bases)
    nlines == 0 && return x
    kernel = _normalize_lines_kernel!(KernelAbstractions.get_backend(x), 256)
    kernel(x, bases, Int(stride), Int(ns), norm1, norm2; ndrange = nlines * ns)
    return x
end

function normalize_copy_lines!(out, out_bases, out_stride::Int, x, x_bases, x_stride::Int, ns::Int, norm1, norm2)
    ns == 0 && return out
    nlines = length(out_bases)
    nlines == 0 && return out
    kernel = _normalize_copy_lines_kernel!(KernelAbstractions.get_backend(out), 256)
    kernel(out, out_bases, Int(out_stride), x, x_bases, Int(x_stride), Int(ns), norm1, norm2; ndrange = nlines * ns)
    return out
end

function lift_step_lines!(x, bases, stride::Int, half::Int, step::GPULiftStep)
    half == 0 && return x
    nlines = length(bases)
    nlines == 0 && return x
    kernel = _lift_step_lines_kernel!(KernelAbstractions.get_backend(x), 256)
    kernel(
        x, bases, Int(stride), Int(half), step.coef, Int(length(step.coef)), Int(step.shift), step.predict;
        ndrange = nlines * half
    )
    return x
end

function lifting_forward_lines!(out, x, bases, stride::Int, ns::Int, steps, norm1, norm2)
    split_lines!(out, bases, stride, x, bases, stride, ns)
    half = ns >> 1
    for step in steps
        lift_step_lines!(out, bases, stride, half, step)
    end
    normalize_lines!(out, bases, stride, ns, norm1, norm2)
    return out
end

function lifting_inverse_lines!(out, x, work, bases, stride::Int, ns::Int, steps, norm1, norm2)
    normalize_copy_lines!(work, bases, stride, x, bases, stride, ns, norm1, norm2)
    half = ns >> 1
    for step in steps
        lift_step_lines!(work, bases, stride, half, step)
    end
    merge_lines!(out, bases, stride, work, bases, stride, ns)
    return out
end

function _dwt!(y::AbstractGPUVector{T}, scheme::GLS, L::Integer, fw::Bool) where {T <: Number}
    n = length(y)
    0 <= L || throw(ArgumentError("L must be positive"))
    sufficientpoweroftwo(y, L) || throw(ArgumentError("size must have a sufficient power of 2 factor"))
    L == 0 && return y

    backend = KernelAbstractions.get_backend(y)
    steps, norm1, norm2 = gpu_makescheme(backend, T, scheme, fw)
    bases = LineBases(1)
    temp = similar(y, T, size(y))
    work = similar(y, T, size(y))

    if fw
        current = y
        output = temp
        ns = n
        for _ in 1:L
            current === output || copyto!(output, current)
            lifting_forward_lines!(output, current, bases, 1, ns, steps, norm1, norm2)
            current, output = output, current
            ns >>= 1
        end
        current === y || copyto!(y, current)
    else
        current = y
        output = temp
        ns = div(n, 2^(L - 1))
        for _ in L:-1:1
            current === output || copyto!(output, current)
            lifting_inverse_lines!(output, current, work, bases, 1, ns, steps, norm1, norm2)
            current, output = output, current
            ns <<= 1
        end
        current === y || copyto!(y, current)
    end

    return y
end

function _dwt!(y::AbstractGPUMatrix{T}, scheme::GLS, L::Integer, fw::Bool) where {T <: Number}
    n = size(y, 1)
    iscube(y) || throw(ArgumentError("array must be square/cube"))
    0 <= L || throw(ArgumentError("L must be positive"))
    sufficientpoweroftwo(y, L) || throw(ArgumentError("size must have a sufficient power of 2 factor"))
    L == 0 && return y

    backend = KernelAbstractions.get_backend(y)
    steps, norm1, norm2 = gpu_makescheme(backend, T, scheme, fw)
    temp = similar(y, T, size(y))
    work = similar(y, T, size(y))

    if fw
        current = y
        nsub = n
        for _ in 1:L
            row_bases = LineBases(nsub, nsub, 1, 0, 1)
            col_bases = LineBases(nsub, nsub, n, 0, 1)
            lifting_forward_lines!(temp, current, row_bases, n, nsub, steps, norm1, norm2)
            lifting_forward_lines!(y, temp, col_bases, 1, nsub, steps, norm1, norm2)
            current = y
            nsub >>= 1
        end
    else
        current = y
        nsub = div(n, 2^(L - 1))
        for _ in L:-1:1
            col_bases = LineBases(nsub, nsub, n, 0, 1)
            row_bases = LineBases(nsub, nsub, 1, 0, 1)
            lifting_inverse_lines!(temp, current, work, col_bases, 1, nsub, steps, norm1, norm2)
            lifting_inverse_lines!(y, temp, work, row_bases, n, nsub, steps, norm1, norm2)
            current = y
            nsub <<= 1
        end
    end

    return y
end

function _dwt!(y::AbstractGPUArray{T, 3}, scheme::GLS, L::Integer, fw::Bool) where {T <: Number}
    n = size(y, 1)
    iscube(y) || throw(ArgumentError("array must be square/cube"))
    0 <= L || throw(ArgumentError("L must be positive"))
    sufficientpoweroftwo(y, L) || throw(ArgumentError("size must have a sufficient power of 2 factor"))
    L == 0 && return y

    backend = KernelAbstractions.get_backend(y)
    steps, norm1, norm2 = gpu_makescheme(backend, T, scheme, fw)
    temp1 = similar(y, T, size(y))
    temp2 = similar(y, T, size(y))
    work = similar(y, T, size(y))
    row_stride = n
    plane_stride = n * n

    if fw
        current = y
        nsub = n
        for _ in 1:L
            plane_bases = LineBases(nsub * nsub, nsub, 1, n, 1)
            row_bases = LineBases(nsub * nsub, nsub, 1, n * n, 1)
            col_bases = LineBases(nsub * nsub, nsub, n, n * n, 1)
            lifting_forward_lines!(temp1, current, plane_bases, plane_stride, nsub, steps, norm1, norm2)
            lifting_forward_lines!(temp2, temp1, row_bases, row_stride, nsub, steps, norm1, norm2)
            lifting_forward_lines!(y, temp2, col_bases, 1, nsub, steps, norm1, norm2)
            current = y
            nsub >>= 1
        end
    else
        current = y
        nsub = div(n, 2^(L - 1))
        for _ in L:-1:1
            col_bases = LineBases(nsub * nsub, nsub, n, n * n, 1)
            row_bases = LineBases(nsub * nsub, nsub, 1, n * n, 1)
            plane_bases = LineBases(nsub * nsub, nsub, 1, n, 1)
            lifting_inverse_lines!(temp1, current, work, col_bases, 1, nsub, steps, norm1, norm2)
            lifting_inverse_lines!(temp2, temp1, work, row_bases, row_stride, nsub, steps, norm1, norm2)
            lifting_inverse_lines!(y, temp2, work, plane_bases, plane_stride, nsub, steps, norm1, norm2)
            current = y
            nsub <<= 1
        end
    end

    return y
end

function _wpt!(
        y::AbstractGPUVector{T}, scheme::GLS, tree::BitVector, fw::Bool,
        tmp::AbstractVector{T} = similar(y, T, 0)
    ) where {T <: Number}
    n = length(y)
    isvalidtree(y, tree) || throw(ArgumentError("invalid tree"))
    !tree[1] && return y

    backend = KernelAbstractions.get_backend(y)
    steps, norm1, norm2 = gpu_makescheme(backend, T, scheme, fw)
    temp = similar(y, T, size(y))
    work = similar(y, T, size(y))
    current = y
    output = temp
    Lmax = maxtransformlevels(n)

    for L in Lmax:-1:1
        Lfw = fw ? (Lmax - L) : (L - 1)
        seglen = detailn(n, Lfw)
        bases_cpu = segment_bases(n, seglen)
        bases = LineBases(n ÷ seglen, n ÷ seglen, seglen, 0, 1)
        copy_lines!(output, bases, 1, current, bases, 1, seglen)
        treeind = 2^Lfw - 1
        active_idx = findall(k -> tree[treeind + k], 1:length(bases_cpu))
        if !isempty(active_idx)
            active_bases = to_device(backend, Int.(bases_cpu[active_idx]))
            if fw
                lifting_forward_lines!(output, current, active_bases, 1, seglen, steps, norm1, norm2)
            else
                lifting_inverse_lines!(output, current, work, active_bases, 1, seglen, steps, norm1, norm2)
            end
        end
        current, output = output, output === y ? temp : y
    end

    current === y || copyto!(y, current)
    return y
end
