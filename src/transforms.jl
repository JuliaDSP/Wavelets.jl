module Transforms
using ..Util, ..WT
using Compat
export  dwt, idwt, dwt!, idwt!,
        wpt, iwpt, wpt!, iwpt!,
        dwtc, idwtc

# TODO Use StridedArray instead of AbstractArray where writing to array.
typealias DWTArray AbstractArray
typealias WPTArray AbstractVector


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
    function ($Xwt){T<:AbstractFloat}(x::DWTArray{T},
                                    filter::OrthoFilter,
                                    L::Integer=maxtransformlevels(x))
        y = Array(T, size(x))
        return ($_Xwt!)(y, x, filter, L, $fw)
    end
    function ($Xwt!){T<:AbstractFloat}(y::DWTArray{T}, x::DWTArray{T},
                                    filter::OrthoFilter,
                                    L::Integer=maxtransformlevels(x))
        return ($_Xwt!)(y, x, filter, L, $fw)
    end
    # lifting
    function ($Xwt){T<:AbstractFloat}(x::DWTArray{T},
                                    scheme::GLS,
                                    L::Integer=maxtransformlevels(x))
        y = Array(T, size(x))
        copy!(y, x)
        return ($_Xwt!)(y, scheme, L, $fw)
    end
    function ($Xwt!){T<:AbstractFloat}(y::DWTArray{T},
                                    scheme::GLS,
                                    L::Integer=maxtransformlevels(x))
        return ($_Xwt!)(y, scheme, L, $fw)
    end
end # begin
end # for


# WPT (wavelet packet transform)

for (Xwt, Xwt!, _Xwt!, fw) in ((:wpt, :wpt!, :_wpt!, true),
                                (:iwpt, :iwpt!, :_wpt!, false))
@eval begin
    function ($Xwt){T<:AbstractFloat}(x::WPTArray{T},
                                    wt::DiscreteWavelet,
                                    L::Integer=maxtransformlevels(x))
        return ($Xwt)(x, wt, maketree(length(x), L, :full))
    end
    # filter
    function ($Xwt){T<:AbstractFloat}(x::WPTArray{T},
                                    filter::OrthoFilter,
                                    tree::BitVector=maketree(x, :full))
        y = Array(T, size(x))
        return ($_Xwt!)(y, x, filter, tree, $fw)
    end
    function ($Xwt!){T<:AbstractFloat}(y::WPTArray{T}, x::WPTArray{T},
                                    filter::OrthoFilter,
                                    tree::BitVector=maketree(x, :full))
        return ($_Xwt!)(y, x, filter, tree, $fw)
    end
    function ($Xwt!){T<:AbstractFloat}(y::WPTArray{T}, x::WPTArray{T},
                                    filter::OrthoFilter,
                                    L::Integer=maxtransformlevels(x))
        return ($Xwt!)(y, x, filter, maketree(length(x), L, :full))
    end
    # lifting
    function ($Xwt){T<:AbstractFloat}(x::WPTArray{T},
                                    scheme::GLS,
                                    tree::BitVector=maketree(x, :full))
        y = Array(T, size(x))
        copy!(y, x)
        return ($_Xwt!)(y, scheme, tree, $fw)
    end
    function ($Xwt!){T<:AbstractFloat}(y::WPTArray{T},
                                    scheme::GLS,
                                    tree::BitVector=maketree(y, :full))
        return ($_Xwt!)(y, scheme, tree, $fw)
    end
    function ($Xwt!){T<:AbstractFloat}(y::WPTArray{T},
                                    scheme::GLS,
                                    L::Integer=maxtransformlevels(x))
        return ($Xwt!)(y, scheme, maketree(length(x), L, :full))
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

@deprecate dwt!(y, x, filter::OrthoFilter, L, fw) (fw ? dwt!(y,x,filter,L) : idwt!(y,x,filter,L))
@deprecate dwt!(y, scheme::GLS, L, fw) (fw ? dwt!(y,scheme,L) : idwt!(y,scheme,L))
@deprecate wpt!(y, x, filter::OrthoFilter, L, fw) (fw ? wpt!(y,x,filter,L) : iwpt!(y,x,filter,L))
@deprecate wpt!(y, scheme::GLS, L, fw) (fw ? wpt!(y,scheme,L) : wpt!(y,scheme,L))

# Int -> Float
for Xwt in (:dwt, :idwt, :dwtc, :idwtc, :wpt, :iwpt)
@eval begin
    ($Xwt){T<:Integer}(x::AbstractArray{T}, args...) = ($Xwt)(float(x), args...)
end
end

# non-exported "out of place" functions
for (Xwt_oop!, Xwt!) in ((:dwt_oop!, :dwt!), (:idwt_oop!, :idwt!))
@eval begin
    # filter
    function ($Xwt_oop!){T<:AbstractFloat}(y::DWTArray{T}, x::DWTArray{T},
                                    filter::OrthoFilter,
                                    L::Integer=maxtransformlevels(x))
        return ($Xwt!)(y, x, filter, L)
    end
    # lifting
    function ($Xwt_oop!){T<:AbstractFloat}(y::DWTArray{T}, x::DWTArray{T},
                                    scheme::GLS,
                                    L::Integer=maxtransformlevels(x))
        copy!(y, x)
        return ($Xwt!)(y, scheme, L)
    end
end # begin
end # for



# column-wise transforms or color images, transform each x[:,...,:,i] separately (default)
# or transform each x[:,...,i,:,...] separately at dim td
for (Xwtc, Xwt) in ((:dwtc, :dwt!), (:idwtc, :idwt!))
@eval begin
    function $Xwtc{T<:AbstractFloat}(x::AbstractArray{T}, wt::DiscreteWavelet, L::Integer, td::Integer=ndims(x))
        dim = ndims(x)
        (1 <= td <= dim) || throw(BoundsError())
        sizex = size(x)
        sizexc = maketfsize(sizex, td)
        y = Array(T, sizex)
        xc = Array(T, sizexc)
        yc = Array(T, sizexc)

        ind = Array(Any, dim)
        for i = 1:dim
            ind[i] = 1:sizex[i]
        end

        for d = 1:sizex[td]
            ind[td] = d
            if dim > 2
                xc[:] = x[ind...]
                ($Xwt)(yc, xc, wt, L)
                y[ind...] = yc
            else  # fast copy methods in 2-D
                Util.copygeneral2!(xc, x, ind...)
                ($Xwt)(yc, xc, wt, L)
                Util.copygeneral1!(y, ind..., yc)
            end
        end
        return y
    end
end
end


# utils

# for dwtc
function maketfsize(t::NTuple, td::Integer)
    s = Array(eltype(t[1]), length(t)-1)
    k = 1
    for i in 1:length(t) # TODO eachindex for tuples
        if i != td
            s[k] = t[i]
            k += 1
        end
    end
    return tuple(s...)
end

# Array with shared memory
function unsafe_vectorslice(A::Array, i::Int, n::Int)
    return pointer_to_array(pointer(A, i), n, false)::Vector
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
