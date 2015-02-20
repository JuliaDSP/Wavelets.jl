module Transforms
using ..Util, ..WaveletTypes
export  dwt, idwt, dwt!, 
        dwtc, idwtc, wpt, 
        iwpt, wpt!

# TODO Use StridedArray instead of AbstractArray where writing to array.

## API

# DWT (discrete wavelet transform)
#dwt(::AbstractVector, ::DiscreteWavelet)
#idwt(::AbstractVector, ::DiscreteWavelet)

# DWTC (column-wise discrete wavelet transform)
#dwtc(::AbstractArray, ::DiscreteWavelet)
#idwtc(::AbstractArray, ::DiscreteWavelet)

# SWT (stationary wavelet transform)
#swt(::AbstractVector, ::DiscreteWavelet)
#iswt(::AbstractVector, ::DiscreteWavelet)

# WPT (wavelet packet transform)
#wpt(::AbstractVector, ::DiscreteWavelet)
#iwpt(::AbstractVector, ::DiscreteWavelet)

# CWT (continuous wavelet transform)
#cwt(::AbstractVector, ::ContinuousWavelet)
#icwt(::AbstractVector, ::ContinuousWavelet)

# CWTFT (continuous wavelet transform via FFT)
#cwtft(::AbstractVector, ::ContinuousWavelet)
#icwtft(::AbstractVector, ::ContinuousWavelet)


## GENERAL TRANSFORM FUNCTIONS

# general functions for all types and dimensions
for Xwt in (:dwt, :idwt, :dwtc, :idwtc, :wpt, :iwpt)
@eval begin
    # assume full transform
    $Xwt(x::AbstractArray, wt::DiscreteWavelet) = $Xwt(x, wt, maxtranformlevels(x))
    # int -> float
    $Xwt{T<:Integer}(x::AbstractArray{T}, wt::DiscreteWavelet, L::Integer) = $Xwt(float(x), wt, L)
end
end

# DWT, WPT methods with Array allocation
for (Xwt, fw, Xwtip) in ((:dwt, true, :dwt!), 
                         (:idwt, false, :dwt!), 
                         (:wpt, true, :wpt!), 
                         (:iwpt, false, :wpt!))
@eval begin
    function $Xwt{T<:FloatingPoint}(x::AbstractArray{T}, wt::DiscreteWavelet, L::Integer)
        y = Array(T, size(x))
        $Xwtip(y, x, wt, L, $fw)
        return y
    end
end
if Xwt == :wpt || Xwt == :iwpt
@eval begin
    function $Xwt{T<:FloatingPoint}(x::AbstractArray{T}, wt::DiscreteWavelet, tree::BitVector)
        y = Array(T, size(x))
        $Xwtip(y, x, wt, tree, $fw)
        return y
    end
end
end
end

# column-wise transforms or color images, transform each x[:,...,:,i] separately (default)
# or transform each x[:,...,i,:,...] separately at dim td
for (Xwtc, Xwt, dir) in ((:dwtc, :dwt!, :true), (:idwtc, :dwt!, :false))
@eval begin
    function $Xwtc{T<:FloatingPoint}(x::AbstractArray{T}, wt::DiscreteWavelet, L::Integer, td::Integer=ndims(x))
        dim = ndims(x)
        @assert 1 <= td <= dim
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
                ($Xwt)(yc, xc, wt, L, $dir)
                y[ind...] = yc
            else  # fast copy methods in 2-D
                Util.copygeneral2!(xc, x, ind...)
                ($Xwt)(yc, xc, wt, L, $dir)
                Util.copygeneral1!(y, ind..., yc)
            end
        end
        return y
    end
end
end

# for dwtc
function maketfsize(t::NTuple, td::Integer)
    s = Array(eltype(t[1]), length(t)-1)
    k = 1
    for i=1:length(t)
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

# filter transforms
include("transforms_filter.jl")

# lifting transforms
include("transforms_lifting.jl")


end # module

