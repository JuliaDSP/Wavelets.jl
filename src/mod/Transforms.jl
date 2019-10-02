module Transforms
export  dwt, idwt, dwt!, idwt!,
        wpt, iwpt, wpt!, iwpt!,
        modwt, imodwt, cwt
using ..Util, ..WT
using FFTW

# TODO Use StridedArray instead of AbstractArray where writing to array.
# TODO change integer dependent wavelets to parametric types (see "Value types", https://docs.julialang.org/en/v1/manual/types/index.html#%22Value-types%22-1)
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

# CWT (continuous wavelet transform)
# cwt(Y::AbstractVector, ::ContinuousWavelet)

@doc """
     wave = cwt(Y::AbstractArray{T,N}, c::CFW{W, S, WT}, daughters, rfftPlan =
             plan_rfft([1]), fftPlan = plan_fft([1])) where {N, T<:Real,
                                                             S<:Real,
                                                             U<:Number,
                                                             W<:WT.WaveletBoundary,
                                                             WT<:Union{<:WT.Morlet,
                                                                       <:WT.Paul}}

  return the continuous wavelet transform along the first axis with averaging.
  `wave`, is (signalLength)×(nscales+1)×(previous dimensions), of type T of
  Y. averagingLength defines the number of octaves (powers of 2) that are
  replaced by an averaging function. This has form averagingType, which can be
  one of `Mother()` or `Dirac()`- in the `Mother()` case, it uses the same form
  as for the wavelets, while the `Dirac` uses a constant window. J1 is the
  total number of scales; default (when J1=NaN, or is negative) is just under
  the maximum possible number, i.e. the log base 2 of the length of the signal,
  times the number of wavelets per octave. If you have sampling information,
  you will need to scale wave by δt^(1/2).

  """
function cwt(Y::AbstractArray{T,N}, c::CFW{W, S, WaTy}, daughters, rfftPlan::AbstractFFTs.Plan =
             plan_rfft([1]), fftPlan = plan_fft([1])) where {N, T<:Real,
                                                             S<:Real,
                                                             W<:WT.WaveletBoundary,
                                                             WaTy<:Union{<:WT.Morlet,
                                                                         <:WT.Paul}}
    # This is for analytic wavelets, so we need to treat the positive and
    # negative frequencies differently, even for real data

    # TODO: complex input version of this
    @assert typeof(N)<:Integer
    # vectors behave a bit strangely, so we reshape them
    if N==1
        Y= reshape(Y,(length(Y), 1))
    end
    n1 = size(Y, 1);
    
    nScales = getNScales(n1, c)
    #....construct time series to analyze, pad if necessary
    x = reflect(Y, boundaryType(c)()) #this function is defined below

    # check if the plans we were given are dummies or not
    if size(rfftPlan)==(1,)
        rfftPlan = plan_rfft(x, 1)
    end
    if size(fftPlan)==(1,)
        fftPlan = plan_fft(x, 1)
    end
    n = size(x, 1)
    
    x̂ = rfftPlan * x
    # If the vector isn't long enough to actually have any other scales, just
    # return the averaging
    if nScales <= 0 || size(daughters,2) == 1
        daughters = daughters[:,1:1]
        nScales = 0
    end

    isAve = (c.averagingLength > 0 && !(typeof(c.averagingType) <: WT.NoAve)) ? 1 : 0

    wave = zeros(Complex{T}, size(x)..., nScales + isAve);  # result array;
    # faster if we put the example index on the outside loop through all scales
    # and compute transform
    actuallyTransform!(wave, daughters,x̂, fftPlan, c.waveType, c.averagingType)
    wave = permutedims(wave, [1, ndims(wave), (2:(ndims(wave)-1))...])
    ax = axes(wave)
    wave = wave[1:n1, ax[2:end]...] 
    if N==1
        wave = dropdims(wave, dims=3)
    end

    return wave
end


function cwt(Y::AbstractArray{T,N}, c::CFW{W, S, WaTy}, daughters, rfftPlan =
             plan_rfft([1])) where {N, T<:Real, S<:Real, U<:Number,
                                    W<:WT.WaveletBoundary, WaTy<:WT.Dog}
    # Dog doesn't need a fft because it is strictly real
    # TODO: complex input version of this
    @assert typeof(N)<:Integer
    # vectors behave a bit strangely, so we reshape them
    if N==1
        Y= reshape(Y,(length(Y), 1))
    end

    n1 = size(Y, 1);
    
    nScales = getNScales(n1, c)
    #....construct time series to analyze, pad if necessary
    x = reflect(Y, boundaryType(c)())

    # check if the plans we were given are dummies or not
    if size(rfftPlan)==(1,)
        rfftPlan = plan_rfft(x, 1)
    end
    n = size(x, 1)

    # If the vector isn't long enough to actually have any other scales, just
    # return the averaging. Or if there's only averaging
    if nScales <= 0 || size(daughters,2) == 1
        daughters = daughters[:,1:1]
        nScales = 0
    end

    x̂ = rfftPlan * x
    
    isAve = (c.averagingLength > 0 && !(typeof(c.averagingType) <: WT.NoAve)) ? 1 : 0

    wave = zeros(Complex{T}, size(x)..., nScales + isAve);  # result array;
    # faster if we put the example index on the outside loop through all scales
    # and compute transform


    actuallyTransform!(wave, daughters,x̂, rfftPlan, c.waveType)
    wave = permutedims(wave, [1, ndims(wave), (2:(ndims(wave)-1))...])
    ax = axes(wave)
    wave = wave[1:n1, ax[2:end]...] 

    if N==1
        wave = dropdims(wave, dims=3)
    end

    return real.(wave)
end



function getNScales(n1, c)
    nOctaves = log2(max(n1, 2)) - c.averagingLength
    nWaveletsInOctave = reverse([max(1, round(Int,
                                              c.scalingFactor/x^(c.decreasing)))
                                 for x = 1:round(Int, nOctaves)])
    nScales = max(sum(nWaveletsInOctave), 0)
end

function reflect(Y, bt)
    n1 = size(Y, 1)
    if bt == WT.padded
        base2 = round(Int,log(n1)/log(2));   # power of 2 nearest to N
        x = cat(Y, zeros(2^(base2+1)-n1, size(Y)[2:end]...), dims=1)
    elseif bt == WT.DEFAULT_BOUNDARY
        x = cat(Y, reverse(Y,dims = 1), dims = 1)
    else
        x = Y
    end
    return x
end

function actuallyTransform!(wave, daughters, x̂, fftPlan, analytic::Union{<:WT.Morlet,
                                                                         <:WT.Paul},
                            averagingType::Union{WT.Father, WT.Dirac})
    outer = axes(x̂)[2:end]
    n1 = size(x̂, 1)
    isSourceOdd = mod(size(wave,1)+1,2)
    # the averaging function isn't analytic, so we need to do both positive and
    # negative frequencies
    tmpWave = x̂ .* daughters[:,1]
    wave[(n1+1):end, outer..., 1] = reverse(conj.(tmpWave[2:end-isSourceOdd, outer...]),dims=1)
    wave[1:n1, outer..., 1] = tmpWave
    wave[:, outer..., 1] = fftPlan \ (wave[:, outer..., 1])  # wavelet transform
    for j in 2:size(daughters,2)
        wave[1:n1, outer..., j] = x̂ .* daughters[:,j]
        wave[:, outer..., j] = fftPlan \ (wave[:, outer..., j])  # wavelet transform
    end
end
function actuallyTransform!(wave, daughters, x̂, fftPlan, analytic::Union{<:WT.Morlet,
                                                                         <:WT.Paul},
                            averagingType::WT.NoAve)
    outer = axes(x̂)[2:end]
    n1 = size(x̂, 1)
    # the averaging function isn't analytic, so we need to do both positive and
    # negative frequencies
    for j in 1:size(daughters,2)
        wave[1:n1, outer..., j] = x̂ .* daughters[:,j]
        wave[:, outer..., j] = fftPlan \ (wave[:, outer..., j])  # wavelet transform
    end
end

function actuallyTransform!(wave, daughters, x̂, rfftPlan, analytic::Union{<:WT.Dog})
    outer = axes(x̂)[2:end]
    n1 = size(x̂, 1)
    for j in 1:size(daughters,2)
        wave[1:n1, outer..., j] = x̂ .* daughters[:,j]
        wave[:, outer..., j] = rfftPlan \ (wave[1:n1, outer..., j])  # wavelet transform
    end
end



function cwt(Y::AbstractArray{T}, c::CFW{W}; J1::Int64=-1, dt::S=NaN, s0::V=NaN) where {T<:Number, S<:Real, V<: Real,
                                                                                        W<:WT.WaveletBoundary}
    daughters,ω = computeWavelets(size(Y, 1), c; J1=J1, dt=dt, s0=s0) 
    return cwt(Y, c, daughters)
end


"""
period,scale, coi = caveats(Y::AbstractArray{T}, c::CFW{W}; J1::S=NaN) where {T<:Real, S<:Real, W<:WT.WaveletBoundary}

returns the period, the scales, and the cone of influence for the given wavelet transform. If you have sampling information, you will need to scale the vector scale appropriately by 1/δt, and the actual transform by δt^(1/2).
"""
function caveats(n1, c::CFW{W}; J1::Int64=-1, dt::S=NaN, s0::V=NaN) where {S<:Real, W<:WT.WaveletBoundary, V <: Real}
    # don't alter scaling with sampling information if it doesn't exists
    fλ = (4*π) / (c.σ[1] + sqrt(2 + c.σ[1]^2))
    if isnan(dt) || (dt<0)
        dt = 1
    end

    # smallest resolvable scale
    if isnan(s0) || (s0<0)
        s0 = 2 * dt / fλ
    end
    sj = s0 * 2.0.^(collect(0:J1)./c.scalingFactor)
    # Fourier equivalent frequencies
    freqs = 1 ./ (fλ .* sj)

    # Determines the cone-of-influence. Note that it is returned as a function
    # of time in Fourier periods. Uses triangualr Bartlett window with
    # non-zero end-points.
    coi = (n1 / 2 .- abs.(collect(0:n1-1) .- (n1 - 1) ./ 2))
    coi = (fλ * dt / sqrt(2)).*coi


    n1 = length(Y);
    # J1 is the total number of elements
    if isnan(J1) || (J1<0)
        J1=floor(Int,(log2(n1))*c.scalingFactor);
    end
    #....construct time series to analyze, pad if necessary
    if boundaryType(c) == WT.ZPBoundary
        base2 = round(Int,log(n1)/log(2));   # power of 2 nearest to N
        n = length(Y)+2^(base2+1)-n1
    elseif boundaryType(c) == WT.PerBoundary
        n = length(Y)*2
    end
    ω = [0:floor(Int, n/2); -floor(Int,n/2)+1:-1]*2π
    period = c.fourierFactor*2 .^((0:J1)/c.scalingFactor)
    scale = [1E-5; 1:((n1+1)/2-1); reverse((1:(n1/2-1)),dims=1); 1E-5]
    coi = c.coi*scale  # COI [Sec.3g]
    return sj, freqs, period, scale, coi
end
cwt(Y::AbstractArray{T}, w::WT.ContinuousWaveletClass; J1::Int64=-1, dt::S=NaN, s0::V=NaN) where {T<:Real, S<:Real, V<:Real} = cwt(Y,CFW(w),J1=J1,dt=dt,s0=s0)
caveats(Y::AbstractArray{T}, w::WT.ContinuousWaveletClass; J1::S=NaN) where {T<: Real, S<: Real} = caveats(Y,CFW(w),J1=J1)
cwt(Y::AbstractArray{T}) where T<:Real = cwt(Y,WT.Morlet())
caveats(Y::AbstractArray{T}) where T<:Real = caveats(Y,WT.Morlet())

"""
icwt(WT::AbstractArray{T}, c::CFW{W}, sj::AbstractArray; dt::S=NaN, dj::V=1/12) where {T<:Complex{Real}, S<:Real, V<:Real, W<:WT.WaveletBoundary}

return the inverse continuous wavelet transform
"""
function icwt(WT::AbstractArray, c::CFW{W}, sj::AbstractArray; dt::S=NaN, dj::V=1/12) where {S<:Real, V<:Real, W<:WT.WaveletBoundary}
    # Torrence and Compo (1998), eq. (11)
    iW = (dj * sqrt(dt) / 0.776 * psi(c, 0)) .* sum((real.(WT) ./ sqrt.(sj)), dims=1)

    return iW
end
cwt(Y::AbstractArray{T}, w::WT.ContinuousWaveletClass, sj::AbstractArray; dt::S=NaN, dj::V=NaN) where {T<:Real, S<:Real, V<:Real} = cwt(Y,CFW(w),sj,dt=dt,dj=dj)

function psi(c::CFW{W}, t::Int64) where W<:WT.WaveletBoundary
    return π^(-0.25) * exp.(im*c.σ[1]*t - t^2 / 2)
end


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

# CWT (continuous wavelet transform directly) TODO: direct if sufficiently small

# TODO: continuous inverse, when defined
#icwt(::AbstractVector, ::ContinuousWavelet)

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
