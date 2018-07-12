module Transforms
export  dwt, idwt, dwt!, idwt!,
        wpt, iwpt, wpt!, iwpt!,
        modwt, imodwt
using ..Util, ..WT
using Compat: AbstractRange, copyto!, undef

# TODO Use StridedArray instead of AbstractArray where writing to array.
const DWTArray = AbstractArray
const WPTArray = AbstractVector
const ValueType = Union{AbstractFloat, Complex}

# DWT

"""
    dwt(x, wt[, L=maxtransformlevels(x)])

Perform a discrete wavelet transform of the array `x`.
The wavelet type `wt` determines the transform type
(filter or lifting) and the wavelet class, see `wavelet`.

The number of transformation levels `L` can be any non-negative
integer such that the size of `x` is divisible by `L`.
Performs the identity transformation if `L==0`.

# Examples
```julia
dwt(x, wavelet(WT.coif6))
```

**See also:** `idwt`, `dwt!`, `wavelet`
"""
function dwt end

"""
    idwt(x, wt[, L=maxtransformlevels(x)])

Perform an inverse discrete wavelet transform of the array `x`,
the inverse of `dwt(x, wt, L)`.

**See also:** `dwt`, `idwt!`
"""
function idwt end

"""
    dwt!(y, x, wt::OrthoFilter[, L=maxtransformlevels(x)])

    dwt!(y, wt::GLS[, L=maxtransformlevels(x)])

Same as `dwt` but without array allocation.
Perform "out of place" transform with a filter, or
a inplace transform with a lifting scheme. The difference
between the filter and lifting methods is due to the
structure of the transform algorithms.

**See also:** `idwt!`
"""
function dwt! end

"""
    idwt!(y, x, wt::OrthoFilter[, L=maxtransformlevels(x)])

    idwt!(y, wt::GLS[, L=maxtransformlevels(x)])

The inverse of `dwt!`.

**See also:** `dwt!`
"""
function idwt! end


# WPT

"""
    wpt

Perform a discrete wavelet packet transform of the array `x`.
**See also:** `dwt`, `wavelet`
"""
function wpt end

"""
    iwpt

Perform an inverse discrete wavelet packet transform of the array `x`.
**See also:** `idwt`, `wavelet`
"""
function iwpt end

"""
    wpt!

Same as `wpt` but without array allocation.
"""
function wpt! end

"""
    iwpt!

Same as `iwpt` but without array allocation.
"""
function iwpt! end


# DWT (discrete wavelet transform)

for (Xwt, Xwt!, _Xwt!, fw) in ((:dwt, :dwt!, :_dwt!, true),
                                (:idwt, :idwt!, :_dwt!, false))
@eval begin
    # filter
    function ($Xwt)(x::DWTArray{T}, filter::OrthoFilter,
                    L::Integer=maxtransformlevels(x)) where T<:ValueType
        y = Array{T}(undef, size(x))
        return ($_Xwt!)(y, x, filter, L, $fw)
    end
    function ($Xwt!)(y::DWTArray{T}, x::DWTArray{T}, filter::OrthoFilter,
                    L::Integer=maxtransformlevels(x)) where T<:ValueType
        return ($_Xwt!)(y, x, filter, L, $fw)
    end
    # lifting
    function ($Xwt)(x::DWTArray{T}, scheme::GLS,
                    L::Integer=maxtransformlevels(x)) where T<:ValueType
        y = Array{T}(undef, size(x))
        copyto!(y, x)
        return ($_Xwt!)(y, scheme, L, $fw)
    end
    function ($Xwt!)(y::DWTArray{T}, scheme::GLS,
                    L::Integer=maxtransformlevels(y)) where T<:ValueType
        return ($_Xwt!)(y, scheme, L, $fw)
    end
end # begin
end # for


# WPT (wavelet packet transform)

for (Xwt, Xwt!, _Xwt!, fw) in ((:wpt, :wpt!, :_wpt!, true),
                                (:iwpt, :iwpt!, :_wpt!, false))
@eval begin
    function ($Xwt)(x::WPTArray{T}, wt::DiscreteWavelet,
                    L::Integer=maxtransformlevels(x)) where T<:ValueType
        return ($Xwt)(x, wt, maketree(length(x), L, :full))
    end
    # filter
    function ($Xwt)(x::WPTArray{T}, filter::OrthoFilter,
                    tree::BitVector=maketree(x, :full)) where T<:ValueType
        y = Array{T}(undef, size(x))
        return ($_Xwt!)(y, x, filter, tree, $fw)
    end
    function ($Xwt!)(y::WPTArray{T}, x::WPTArray{T},
                    filter::OrthoFilter) where T<:ValueType
        return ($Xwt!)(y, x, filter, maketree(x, :full))
    end
    function ($Xwt!)(y::WPTArray{T}, x::WPTArray{T},
                    filter::OrthoFilter, L::Integer) where T<:ValueType
        return ($Xwt!)(y, x, filter, maketree(length(x), L, :full))
    end
    function ($Xwt!)(y::WPTArray{T}, x::WPTArray{T},
                    filter::OrthoFilter, tree::BitVector) where T<:ValueType
        return ($_Xwt!)(y, x, filter, tree, $fw)
    end
    # lifting
    function ($Xwt)(x::WPTArray{T}, scheme::GLS,
                    tree::BitVector=maketree(x, :full)) where T<:ValueType
        y = Array{T}(undef, size(x))
        copyto!(y, x)
        return ($_Xwt!)(y, scheme, tree, $fw)
    end
    function ($Xwt!)(y::WPTArray{T}, scheme::GLS) where T<:ValueType
        return ($Xwt!)(y, scheme, maketree(y, :full))
    end
    function ($Xwt!)(y::WPTArray{T}, scheme::GLS, L::Integer) where T<:ValueType
        return ($Xwt!)(y, scheme, maketree(length(x), L, :full))
    end
    function ($Xwt!)(y::WPTArray{T}, scheme::GLS, tree::BitVector) where T<:ValueType
        return ($_Xwt!)(y, scheme, tree, $fw)
    end
end # begin
end # for


# DWTC (column-wise discrete wavelet transform)
#dwtc(::AbstractArray, ::DiscreteWavelet)
#idwtc(::AbstractArray, ::DiscreteWavelet)

# SWT (stationary wavelet transform)
#swt(::AbstractVector, ::DiscreteWavelet)
#iswt(::AbstractVector, ::DiscreteWavelet)

# CWT (continuous wavelet transform)
#cwt(::AbstractVector, ::ContinuousWavelet)
#icwt(::AbstractVector, ::ContinuousWavelet)

# CWTFT (continuous wavelet transform via FFT)
#cwtft(::AbstractVector, ::ContinuousWavelet)
#icwtft(::AbstractVector, ::ContinuousWavelet)

# Int -> Float
for Xwt in (:dwt, :idwt, :dwtc, :idwtc, :wpt, :iwpt)
    @eval $Xwt(x::AbstractArray{<:Integer}, args...) = $Xwt(float(x), args...)
end

# non-exported "out of place" functions
for (Xwt_oop!, Xwt!) in ((:dwt_oop!, :dwt!), (:idwt_oop!, :idwt!))
@eval begin
    # filter
    function ($Xwt_oop!)(y::DWTArray{T}, x::DWTArray{T}, filter::OrthoFilter,
                        L::Integer=maxtransformlevels(x)) where T<:ValueType
        return ($Xwt!)(y, x, filter, L)
    end
    # lifting
    function ($Xwt_oop!)(y::DWTArray{T}, x::DWTArray{T}, scheme::GLS,
                        L::Integer=maxtransformlevels(x)) where T<:ValueType
        copyto!(y, x)
        return ($Xwt!)(y, scheme, L)
    end
end # begin
end # for

# Array with shared memory
function unsafe_vectorslice(A::Array{T}, i::Int, n::Int) where T
    return unsafe_wrap(Array, pointer(A, i), n)::Vector{T}
end

# linear indices of start of rows/cols/planes
# 2-D, size(A) = (m,n)
row_idx(i, m) = i
col_idx(i, m) = 1 + (i-1)*m
# 3-D, size(A) = (m,n,d)
row_idx(i, j, m, n=m) = row_idx(i, n) + (j-1)*n*m
col_idx(i, j, m, n=m) = col_idx(i, m) + (j-1)*n*m
plane_idx(i, j, m) = i + (j-1)*m

# filter transforms
include("transforms_filter.jl")

# lifting transforms
include("transforms_lifting.jl")

include("transforms_maximal_overlap.jl")

end # module
