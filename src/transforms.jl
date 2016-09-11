module Transforms
using ..Util, ..WT
using Compat
export  dwt, idwt, dwt!, idwt!,
        wpt, iwpt, wpt!, iwpt!

# TODO Use StridedArray instead of AbstractArray where writing to array.
typealias DWTArray AbstractArray
typealias WPTArray AbstractVector
typealias ValueType Union{AbstractFloat, Complex}

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
    function ($Xwt){T<:ValueType}(x::DWTArray{T},
                                    filter::OrthoFilter,
                                    L::Integer=maxtransformlevels(x))
        y = Array(T, size(x))
        return ($_Xwt!)(y, x, filter, L, $fw)
    end
    function ($Xwt!){T<:ValueType}(y::DWTArray{T}, x::DWTArray{T},
                                    filter::OrthoFilter,
                                    L::Integer=maxtransformlevels(x))
        return ($_Xwt!)(y, x, filter, L, $fw)
    end
    # lifting
    function ($Xwt){T<:ValueType}(x::DWTArray{T},
                                    scheme::GLS,
                                    L::Integer=maxtransformlevels(x))
        y = Array(T, size(x))
        copy!(y, x)
        return ($_Xwt!)(y, scheme, L, $fw)
    end
    function ($Xwt!){T<:ValueType}(y::DWTArray{T},
                                    scheme::GLS,
                                    L::Integer=maxtransformlevels(y))
        return ($_Xwt!)(y, scheme, L, $fw)
    end
end # begin
end # for


# WPT (wavelet packet transform)

for (Xwt, Xwt!, _Xwt!, fw) in ((:wpt, :wpt!, :_wpt!, true),
                                (:iwpt, :iwpt!, :_wpt!, false))
@eval begin
    function ($Xwt){T<:ValueType}(x::WPTArray{T},
                                    wt::DiscreteWavelet,
                                    L::Integer=maxtransformlevels(x))
        return ($Xwt)(x, wt, maketree(length(x), L, :full))
    end
    # filter
    function ($Xwt){T<:ValueType}(x::WPTArray{T},
                                    filter::OrthoFilter,
                                    tree::BitVector=maketree(x, :full))
        y = Array(T, size(x))
        return ($_Xwt!)(y, x, filter, tree, $fw)
    end
    function ($Xwt!){T<:ValueType}(y::WPTArray{T}, x::WPTArray{T},
                                    filter::OrthoFilter)
        return ($Xwt!)(y, x, filter, maketree(x, :full))
    end
    function ($Xwt!){T<:ValueType}(y::WPTArray{T}, x::WPTArray{T},
                                    filter::OrthoFilter,
                                    L::Integer)
        return ($Xwt!)(y, x, filter, maketree(length(x), L, :full))
    end
    function ($Xwt!){T<:ValueType}(y::WPTArray{T}, x::WPTArray{T},
                                    filter::OrthoFilter,
                                    tree::BitVector)
        return ($_Xwt!)(y, x, filter, tree, $fw)
    end
    # lifting
    function ($Xwt){T<:ValueType}(x::WPTArray{T},
                                    scheme::GLS,
                                    tree::BitVector=maketree(x, :full))
        y = Array(T, size(x))
        copy!(y, x)
        return ($_Xwt!)(y, scheme, tree, $fw)
    end
    function ($Xwt!){T<:ValueType}(y::WPTArray{T},
                                    scheme::GLS)
        return ($Xwt!)(y, scheme, maketree(y, :full))
    end
    function ($Xwt!){T<:ValueType}(y::WPTArray{T},
                                    scheme::GLS,
                                    L::Integer)
        return ($Xwt!)(y, scheme, maketree(length(x), L, :full))
    end
    function ($Xwt!){T<:ValueType}(y::WPTArray{T},
                                    scheme::GLS,
                                    tree::BitVector)
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
    @eval $Xwt{T<:Integer}(x::AbstractArray{T}, args...) = $Xwt(float(x), args...)
end

# non-exported "out of place" functions
for (Xwt_oop!, Xwt!) in ((:dwt_oop!, :dwt!), (:idwt_oop!, :idwt!))
@eval begin
    # filter
    function ($Xwt_oop!){T<:ValueType}(y::DWTArray{T}, x::DWTArray{T},
                                    filter::OrthoFilter,
                                    L::Integer=maxtransformlevels(x))
        return ($Xwt!)(y, x, filter, L)
    end
    # lifting
    function ($Xwt_oop!){T<:ValueType}(y::DWTArray{T}, x::DWTArray{T},
                                    scheme::GLS,
                                    L::Integer=maxtransformlevels(x))
        copy!(y, x)
        return ($Xwt!)(y, scheme, L)
    end
end # begin
end # for

# Array with shared memory
function unsafe_vectorslice{T}(A::Array{T}, i::Int, n::Int)#::Vector{T}
    return unsafe_wrap(Array, pointer(A, i), n, false)::Vector{T}
end

# 2-D
row_idx(i, n) = i
col_idx(i, n) = 1 + (i-1)*n
# 3-D
row_idx(i, j, n) = row_idx(i, n) + (j-1)*n*n
col_idx(i, j, n) = col_idx(i, n) + (j-1)*n*n
hei_idx(i, j, n) = i + (j-1)*n

# filter transforms
include("transforms_filter.jl")

# lifting transforms
include("transforms_lifting.jl")


end # module
