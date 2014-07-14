module Threshold
using ..Util, ..POfilters, ..FilterTransforms, ..LiftingSchemes, ..LiftingTransforms
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
    # treshold without parameters
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
function VisuShrink(n::Real)
    return VisuShrink(thf("hard"), sqrt(2*log(n)))
end

const out_type = Float64                # for non inplace functions
const def_wavelet = POfilter("sym5")    # default wavelet type

function denoise{T<:WaveletType,S<:DNFT}(x::AbstractArray;  
                                    wt::Union(T,Nothing)=def_wavelet, 
                                    level::Int=max(nscales(size(x,1))-6,1),
                                    dnt::S=VisuShrink(size(x,1)),
                                    sigma::Real=noisest(x, wt=wt) )
    if wt == nothing
        y = copy(x)
        dnt.f(y, sigma*dnt.t)
    else
        L = nscales(size(x,1)) - level
        
        y = fwt(x, L, wt)
        # threshold
        dnt.f(y, sigma*dnt.t)
        
        if T <: WaveletLS
            dwt!(y, L, wt, false)
        else
            y = iwt(y, L, wt)
        end
    end
    return y
end

# estimate the std. dev. of the signal noise, assuming Gaussian distribution
function noisest{T<:WaveletType}(x::AbstractArray; wt::Union(T,Nothing)=def_wavelet)
    if wt == nothing
        y = copy(x)
    else
        y = fwt(x, 1, wt)
    end
    ind = detailrange(maxlevel(size(y,1)))
    return mad(y[ind])/0.6745
end
# Median absolute deviation
function mad(x::AbstractArray)
    y = copy(x)
    m = median!(y)
    for i in 1:length(y)
        y[i] = abs(y[i]-m)
    end
    return median!(y, checknan=false)
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
    m < 0 && error("m negative")
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
    t < 0 && error("t negative")
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
    t < 0 && error("t negative")
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
    t < 0 && error("t negative")
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
    t < 0 && error("t negative")
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

