module Transforms
using ..Util, ..WT
using Compat
export  dwt, idwt, dwt!, idwt!,
        wpt, iwpt, wpt!, iwpt!,
        dwtc, idwtc, cwtft
VERSION < v"0.4-" && using Docile

# TODO Use StridedArray instead of AbstractArray where writing to array.
typealias DWTArray AbstractArray
typealias WPTArray AbstractVector


# DWT (discrete wavelet transform)

for (Xwt, Xwt!, _Xwt!, fw) in ((:dwt, :dwt!, :_dwt!, true),
                                (:idwt, :idwt!, :_dwt!, false))
@eval begin
    # filter
    function ($Xwt){T<:FloatingPoint}(x::DWTArray{T},
                                    filter::OrthoFilter,
                                    L::Integer=maxtransformlevels(x))
        y = Array(T, size(x))
        return ($_Xwt!)(y, x, filter, L, $fw)
    end
    function ($Xwt!){T<:FloatingPoint}(y::DWTArray{T}, x::DWTArray{T},
                                    filter::OrthoFilter,
                                    L::Integer=maxtransformlevels(x))
        return ($_Xwt!)(y, x, filter, L, $fw)
    end
    # lifting
    function ($Xwt){T<:FloatingPoint}(x::DWTArray{T},
                                    scheme::GLS,
                                    L::Integer=maxtransformlevels(x))
        y = Array(T, size(x))
        copy!(y, x)
        return ($_Xwt!)(y, scheme, L, $fw)
    end
    function ($Xwt!){T<:FloatingPoint}(y::DWTArray{T},
                                    scheme::GLS,
                                    L::Integer=maxtransformlevels(x))
        return ($_Xwt!)(y, scheme, L, $fw)
    end
end # begin
end # for

@doc """
`dwt(x, wt[, L=maxtransformlevels(x)])`

Perform a discrete wavelet transform of the array `x`.
The wavelet type `wt` determines the transform type
(filter or lifting) and the wavelet class, see `wavelet`.

The number of transformation levels `L` can be any non-negative
integer such that the size of `x` is divisible by `L`.
Performs the identity transformation if `L==0`.

**Example:**
```julia
dwt(x, wavelet(WT.coif6))
```

**See also:** `idwt`, `dwt!`, `wavelet`
""" -> dwt

@doc """
`idwt(x, wt[, L=maxtransformlevels(x)])`

Perform an inverse discrete wavelet transform of the array `x`,
the inverse of `dwt(x, wt, L)`.

**See also:** `dwt`, `idwt!`
""" -> idwt

@doc """
`dwt!(y, x, wt::OrthoFilter[, L=maxtransformlevels(x)])`

`dwt!(y, wt::GLS[, L=maxtransformlevels(x)])`

Same as `dwt` but without array allocation.
Perform "out of place" transform with a filter, or
a inplace transform with a lifting scheme. The difference
between the filter and lifting methods is due to the
structure of the transform algorithms.

**See also:** `idwt!`
""" -> dwt!

@doc """
`idwt!(y, x, wt::OrthoFilter[, L=maxtransformlevels(x)])`

`idwt!(y, wt::GLS[, L=maxtransformlevels(x)])`

The inverse of `dwt!`.

**See also:** `dwt!`
""" -> idwt!


# WPT (wavelet packet transform)

for (Xwt, Xwt!, _Xwt!, fw) in ((:wpt, :wpt!, :_wpt!, true),
                                (:iwpt, :iwpt!, :_wpt!, false))
@eval begin
    function ($Xwt){T<:FloatingPoint}(x::WPTArray{T},
                                    wt::DiscreteWavelet,
                                    L::Integer=maxtransformlevels(x))
        return ($Xwt)(x, wt, maketree(length(x), L, :full))
    end
    # filter
    function ($Xwt){T<:FloatingPoint}(x::WPTArray{T},
                                    filter::OrthoFilter,
                                    tree::BitVector=maketree(x, :full))
        y = Array(T, size(x))
        return ($_Xwt!)(y, x, filter, tree, $fw)
    end
    function ($Xwt!){T<:FloatingPoint}(y::WPTArray{T}, x::WPTArray{T},
                                    filter::OrthoFilter,
                                    tree::BitVector=maketree(x, :full))
        return ($_Xwt!)(y, x, filter, tree, $fw)
    end
    function ($Xwt!){T<:FloatingPoint}(y::WPTArray{T}, x::WPTArray{T},
                                    filter::OrthoFilter,
                                    L::Integer=maxtransformlevels(x))
        return ($Xwt!)(y, x, filter, maketree(length(x), L, :full))
    end
    # lifting
    function ($Xwt){T<:FloatingPoint}(x::WPTArray{T},
                                    scheme::GLS,
                                    tree::BitVector=maketree(x, :full))
        y = Array(T, size(x))
        copy!(y, x)
        return ($_Xwt!)(y, scheme, tree, $fw)
    end
    function ($Xwt!){T<:FloatingPoint}(y::WPTArray{T},
                                    scheme::GLS,
                                    tree::BitVector=maketree(y, :full))
        return ($_Xwt!)(y, scheme, tree, $fw)
    end
    function ($Xwt!){T<:FloatingPoint}(y::WPTArray{T},
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

## After Torrence & Compo (1998):

function cwtft{T<:Real}(Y::AbstractArray{T},dt::Number; pad::Bool=false,dj::Number=0.25,s0::Number=2.*dt,J1::Number=-1,mother::WT.ContinuousWavelet=WT.morlet,param::Number=WT.sparam(mother))

#Y=Y[:];
n1 = length(Y);

if J1 == -1
        J1=floor(Int,(log(n1*dt/s0)/log(2.))/dj);
end
#....construct time series to analyze, pad if necessary
x = Y - mean(Y);
if pad
        base2 = floor(Int,log(n1)/log(2) + 0.4999);   # power of 2 nearest to N
        x = [x, zeros(2^(floor(Int,base2)+1)-n1)];
end
n = length(x);

#....construct wavenumber array used in transform [Eqn(5)]
k = 1:floor(Int,n/2);
k = [0.;  k;  -k[floor(Int,(n-1)/2):-1:1]]*((2*pi)/(n*dt));
#....compute FFT of the (padded) time series
f = fft(x);    # [Eqn(3)]
#....construct SCALE array & empty PERIOD & WAVE arrays
scale = s0*2.^((0:J1)*dj);

wave = zeros(Complex,J1+1,n);  # define the wavelet array


# loop through all scales and compute transform
for a1 in 1:J1+1
        daughter=WT.Daughter(mother,scale[a1],k,param,n)
        wave[a1,:] = ifft(f.*daughter)  # wavelet transform[Eqn(4)]
end
fourier_factor=WT.FourierFactor(mother,param);
period = fourier_factor*scale;
coi = WT.COI(mother,fourier_factor).*dt*[1E-5; 1:((n1+1)/2-1); flipdim((1:(n1/2-1)),1); 1E-5];  # COI [Sec.3g]
wave = wave[:,1:n1];  # get rid of padding before returning


return wave,period,scale,coi
end

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
    function ($Xwt_oop!){T<:FloatingPoint}(y::DWTArray{T}, x::DWTArray{T},
                                    filter::OrthoFilter,
                                    L::Integer=maxtransformlevels(x))
        return ($Xwt!)(y, x, filter, L)
    end
    # lifting
    function ($Xwt_oop!){T<:FloatingPoint}(y::DWTArray{T}, x::DWTArray{T},
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
    function $Xwtc{T<:FloatingPoint}(x::AbstractArray{T}, wt::DiscreteWavelet, L::Integer, td::Integer=ndims(x))
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
