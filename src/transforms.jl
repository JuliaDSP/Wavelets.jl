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
    $Xwt(x::AbstractArray, wt::DiscreteWavelet) = $Xwt(x, wt, nscales(size(x,1)))
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

# column-wise transforms or color images, transform each x[:,...,:,i] separately
for (Xwtc, Xwt) in ((:dwtc, :dwt), (:idwtc, :idwt))
@eval begin
    function $Xwtc{T<:FloatingPoint}(x::AbstractArray{T}, wt::DiscreteWavelet, L::Integer)
        dim = ndims(x)
        cn = size(x, dim)
        y = Array(eltype(x), size(x))
        
        ind = Array(Any, dim)
        for i = 1:dim
            ind[i] = 1:size(x, i)
        end
        
        for d = 1:cn
            ind[dim] = d
            xc = reshape(x[ind...], size(x)[1:end-1]...)
            y[ind...] = $Xwt(xc, wt, L)
        end
        return y
    end
end
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

