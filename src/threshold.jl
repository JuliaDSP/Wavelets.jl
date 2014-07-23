module Threshold
using ..Util, ..WaveletTypes, ..Transforms
export 
    # denoising types
    DNFT,
    VisuShrink,
    # denoising functions
    denoise,
    noisest,
    
    thf,
    # threshold with parameter m
    biggestterms!,
    biggestterms,
    # threshold with parameter t
    thresholdhard!,
    thresholdhard,
    thresholdsoft!,
    thresholdsoft,
    thresholdsemisoft!,
    thresholdsemisoft,
    thresholdsemistein!,
    thresholdsemistein,
    # threshold without parameters
    thresholdneg!,
    thresholdneg,
    thresholdpos!,
    thresholdpos

# thresholding and denoising utilities

abstract DNFT

type VisuShrink <: DNFT
    f::Function     # thresholding function (inplace)
    t::Real         # threshold for noise level sigma=1, use sigma*t in application
end
# define type for signal length n
function VisuShrink(n::Int)
    return VisuShrink(thf("hard"), sqrt(2*log(n)))
end

const out_type = Float64                # for non inplace functions
const def_wavelet = POfilter("sym5")    # default wavelet type

# denoise signal x by thresholding in wavelet space
function denoise{T<:WaveletType,S<:DNFT}(x::AbstractArray;  
                                    wt::Union(T,Nothing)=def_wavelet, 
                                    level::Int=max(nscales(size(x,1))-6,1),
                                    dnt::S=VisuShrink(size(x,1)),
                                    sigma::Real=noisest(x, wt=wt),
                                    TI::Bool=false,
                                    nspin::Union(Int,Tuple)=tuple([8 for i=1:length(size(x))]...) )
    @assert iscube(x)
    
    if TI
        wt == nothing && error("TI not supported with wt=nothing")
        y = zeros(eltype(x), size(x))
        L = level2tl(size(x,1),level)
        pns = prod(nspin)
        
        if ndims(x)==1
            z = Array(eltype(x), size(x))
            T<:GPLS && (tmp=Array(eltype(x),length(x)>>2))
            T<:POfilter && (xt=Array(eltype(x),length(x)))
            for i = 1:pns
                shift = nspin2circ(nspin, i)[1]
                circshift!(z, x, shift)
                
                if T<:GPLS
                    dwt!(z, L, wt, true, tmp)
                    dnt.f(z, sigma*dnt.t)   # threshold
                    dwt!(z, L, wt, false, tmp)
                elseif T<:POfilter
                    dwt!(xt, z, L, wt, true)
                    dnt.f(xt, sigma*dnt.t)   # threshold
                    dwt!(z, xt, L, wt, false)
                else
                    dwt!(z, L, wt, true)
                    dnt.f(z, sigma*dnt.t)   # threshold
                    dwt!(z, L, wt, false)
                end
               shiftadd!(y,z,-shift)
            end
        else # ndims > 1
            for i = 1:pns
                shift = nspin2circ(nspin, i)
                z = circshift(x, shift)
                
                dwt!(z, L, wt, true)
                dnt.f(z, sigma*dnt.t)   # threshold
                dwt!(z, L, wt, false)

                z = circshift(z, -shift)
                broadcast!(+,y,y,z)
            end
        end
        scale!(y,1/pns)
    else
        if wt == nothing
            y = copy(x)
            dnt.f(y, sigma*dnt.t)
        else
            L = level2tl(size(x,1),level)
            y = fwt(x, L, wt)
            dnt.f(y, sigma*dnt.t)   # threshold
            dwt!(y, L, wt, false)
        end
    end
    
    return y
end
# shift z and add to y
function shiftadd!(y,z,shift)
    @assert shift<=0
    @assert length(y)==length(z)
    n = length(y)
    for i = 1:n+shift
        @inbounds y[i] += z[i-shift]
    end
    sh = n+shift
    for i = n+shift+1:n
        @inbounds y[i] += z[i-sh]
    end
    return y
end

# estimate the std. dev. of the signal noise, assuming Gaussian distribution
function noisest{T<:WaveletType}(x::AbstractArray; wt::Union(T,Nothing)=def_wavelet)
    if wt == nothing
        y = x
    else
        y = fwt(x, 1, wt)
    end
    ind = detailrange(maxlevel(size(y,1)))
    dr = y[ind]
    return mad!(dr)/0.6745
end
# Median absolute deviation
function mad!(y::AbstractArray)
    m = median!(y)
    for i in 1:length(y)
        y[i] = abs(y[i]-m)
    end
    return median!(y, checknan=false)
end
function mad(x::AbstractArray)
    y = copy(x)
    mad!(y)
end

# convert index i to a circshift array starting at 0 shift
function nspin2circ(nspin::Union(Int,Tuple), i::Int)
    typeof(nspin) == Int && (nspin = (nspin,))
    c1 = ind2sub(nspin,i)
    c = Array(Int,length(c1))
    for k = 1:length(c1)
        c[k] = c1[k]-1
    end
    return c
end

# return an inplace threshold function
function thf(th::String="hard")
    if th=="hard"
        return thresholdhard!
    elseif th=="soft"
        return thresholdsoft!
    elseif th=="semisoft"
        return thresholdsemisoft!
    elseif th=="stein"
        return thresholdstein!
    end
    error("threshold ", th, " not defined")
end


# WITH 1 PARAMETER t OR m

# biggest m-term approximation (best m-term approximation for orthogonal transforms)
# returns a m-sparse array
function biggestterms!(x::AbstractArray, m::Int)
    @assert m >= 0
    n = length(x)
    m > n && (m = n)
    ind = sortperm(sub(x,1:n), alg=QuickSort, by=abs)
    @inbounds begin
        for i = 1:n-m
            x[ind[i]] = 0
        end
    end
    return x
end

# hard
function thresholdhard!(x::AbstractArray, t::Real)
    @assert t >= 0
    @inbounds begin
        for i = 1:length(x)
            if abs(x[i]) <= t
                x[i] = 0
            end
        end
    end
    return x
end

# soft
function thresholdsoft!(x::AbstractArray, t::Real)
    @assert t >= 0
    @inbounds begin
        for i = 1:length(x)
            sh = abs(x[i]) - t
            if sh < 0
                x[i] = 0
            else
                x[i] = sign(x[i])*sh
            end
        end
    end
    return x
end

# semisoft
function thresholdsemisoft!(x::AbstractArray, t::Real)
    @assert t >= 0
    @inbounds begin
        for i = 1:length(x)
            if x[i] <= 2*t
                sh = abs(x[i]) - t
                if sh < 0
                    x[i] = 0
                elseif sh - t < 0
                    x[i] = sign(x[i])*sh*2
                end
            end
        end
    end
    return x
end

# stein
function thresholdstein!(x::AbstractArray, t::Real)
    @assert t >= 0
    @inbounds begin
        for i = 1:length(x)
            sh = 1 - t*t/(x[i]*x[i])
            if sh < 0
                x[i] = 0
            else
                x[i] = x[i]*sh
            end
        end
    end
    return x
end

# the non inplace functions
for (fn,fn!) in (   (:biggestterms,     :biggestterms!),
                    (:thresholdhard,    :thresholdhard!),
                    (:thresholdsoft,    :thresholdsoft!),
                    (:thresholdsemisoft,:thresholdsemisoft!),
                    (:thresholdstein,   :thresholdstein!)
                 )
@eval begin
function ($fn)(x::AbstractArray, t::Real) 
    y = Array(out_type, size(x))
    return ($fn!)(copy!(y,x), t)
end
end # eval begin
end #for


# WITHOUT PARAMETERS

# shrink negative elements to 0
function thresholdneg!(x::AbstractArray)
    @inbounds begin
        for i = 1:length(x)
            if x[i] < 0
                x[i] = 0
            end
        end
    end
    return x
end

# shrink positive elements to 0
function thresholdpos!(x::AbstractArray)
    @inbounds begin
        for i = 1:length(x)
            if x[i] > 0
                x[i] = 0
            end
        end
    end
    return x
end

# the non inplace functions
for (fn,fn!) in (   (:thresholdneg, :thresholdneg!),
                    (:thresholdpos, :thresholdpos!)
                 )
@eval begin
function ($fn)(x::AbstractArray) 
    y = Array(out_type, size(x))
    return ($fn!)(copy!(y,x))
end
end # eval begin
end #for

end

