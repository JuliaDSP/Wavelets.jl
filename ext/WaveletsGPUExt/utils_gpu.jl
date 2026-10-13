# GPU utility kernels and helpers

function to_device(backend, x::AbstractArray)
    out = KernelAbstractions.allocate(backend, eltype(x), size(x))
    copyto!(out, x)
    return out
end

"""
    LineBases(nlines, inner_count, inner_stride, outer_stride, offset)

Lightweight descriptor that computes the starting index for each "line" via
integer arithmetic instead of materializing a GPU index array.

Every line-based kernel indexes data as `base + (pos - 1) * stride`, where
the base for line `t` (1-based) is:

    offset + ((t-1) % inner_count) * inner_stride
           + ((t-1) ÷ inner_count) * outer_stride
"""
struct LineBases
    nlines::Int
    inner_count::Int
    inner_stride::Int
    outer_stride::Int
    offset::Int
end

# Single line starting at `offset`
LineBases(offset::Int) = LineBases(1, 1, 0, 0, offset)

Base.length(b::LineBases) = b.nlines
Adapt.adapt_structure(to, lb::LineBases) = lb

@inline function get_base(b::LineBases, line::Int)
    t = line - 1
    q, r = divrem(t, b.inner_count)
    return b.offset + r * b.inner_stride + q * b.outer_stride
end
@inline get_base(b::AbstractArray, line::Int) = @inbounds b[line]

@kernel function _copy_lines_kernel!(
        out, out_bases::LineBases, out_stride::Int,
    @Const(x), x_bases::LineBases,   x_stride::Int,
        line_len::Int
    )
    idx = @index(Global)
    line = ((idx - 1) ÷ line_len) + 1
    pos = idx - (line - 1) * line_len
    @inbounds out[get_base(out_bases, line) + (pos - 1) * out_stride] =
        x[get_base(x_bases, line) + (pos - 1) * x_stride]
end

function copy_lines!(
    out::AbstractGPUArray{<:Number}, out_bases::LineBases, out_stride::Int,
    x::AbstractGPUArray{<:Number},     x_bases::LineBases,   x_stride::Int,
    line_len::Int
)
    line_len == 0 && return
    nlines = length(out_bases)
    nlines == 0 && return
    kernel! = _copy_lines_kernel!(KernelAbstractions.get_backend(out))
    kernel!(out, out_bases, Int(out_stride), x, x_bases, Int(x_stride), Int(line_len); ndrange = nlines * line_len)
    return
end

# segment_bases is kept for WPT where irregular subsets of lines are needed
segment_bases(n::Int, seglen::Int) = Int[1 + (k - 1) * seglen for k in 1:(n ÷ seglen)]
