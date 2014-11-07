module Threshold
using ..Util, ..WaveletTypes, ..Transforms
export 
    # denoising types
    DNFT,
    VisuShrink,
    # denoising functions
    denoise,
    noisest,
    matchingpursuit,
    
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

const DEF_WAVELET = waveletfilter("sym5")    # default wavelet type

# denoise signal x by thresholding in wavelet space
function denoise{T<:DiscreteWavelet,S<:DNFT}(x::AbstractArray;  
                                    wt::Union(T,Nothing)=DEF_WAVELET, 
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
                    dwt!(z, wt, L, true, tmp)
                    dnt.f(z, sigma*dnt.t)   # threshold
                    dwt!(z, wt, L, false, tmp)
                elseif T<:POfilter
                    dwt!(xt, z, wt, L, true)
                    dnt.f(xt, sigma*dnt.t)   # threshold
                    dwt!(z, xt, wt, L, false)
                else
                    dwt!(z, wt, L, true)
                    dnt.f(z, sigma*dnt.t)   # threshold
                    dwt!(z, wt, L, false)
                end
               shiftadd!(y,z,-shift)
            end
        else # ndims > 1
            for i = 1:pns
                shift = nspin2circ(nspin, i)
                z = circshift(x, shift)
                
                dwt!(z, wt, L, true)
                dnt.f(z, sigma*dnt.t)   # threshold
                dwt!(z, wt, L, false)

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
            y = dwt(x, wt, L)
            dnt.f(y, sigma*dnt.t)   # threshold
            dwt!(y, wt, L, false)
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
function noisest{T<:DiscreteWavelet}(x::AbstractArray; wt::Union(T,Nothing)=DEF_WAVELET)
    if wt == nothing
        y = x
    else
        y = dwt(x, wt, 1)
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


# Matching Pursuit
# see: Mallat (2009) p.642 "A wavelet tour of signal processing"
# find sparse vector y such that ||x - f(y)|| < tol approximately
# f is the operation of a M by N (M<N) dictionary/matrix
# ft is a function defining the transpose of f
function matchingpursuit(x::AbstractVector, f::Function, ft::Function, tol::Real, nmax::Int=-1, oop::Bool=false, N::Int=0)
    @assert nmax >= -1
    @assert tol > 0
    r = x
    n = 1
    
    if !oop
        y = zeros(eltype(x), length(ft(x)))
    else # out of place functions f and ft
        y = zeros(eltype(x), N)
        tmp = Array(eltype(x), N)
        ftr = Array(eltype(x), N)
        aphi = Array(eltype(x), length(x))
    end
    spat = zeros(eltype(x), length(y))  # sparse for atom computation
    nmax == -1 && (nmax = length(y))
    
    while vecnorm(r) > tol && n <= nmax
        # find largest inner product
        !oop && (ftr = ft(r))
        oop  && ft(ftr, r, tmp)
        i = findmaxabs(ftr)
        
        # project on i-th atom
        spat[i] = ftr[i]
        !oop && (aphi = f(spat))
        oop  && f(aphi, spat, tmp)
        spat[i] = 0
        
        # update residual, r = r - aphi
        broadcast!(-, r, r, aphi)
        
        y[i] += ftr[i]
        n += 1
    end
    return y
end
function findmaxabs(x::AbstractVector)
    m = abs(x[1])
    k = 1
    @inbounds for i = 1:length(x)
        if abs(x[i]) > m
            k = i
            m = abs(x[i])
        end
    end
    return k
end




# WITH 1 PARAMETER t OR m

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
function ($fn){T<:Number}(x::AbstractArray{T}, t::Real) 
    y = Array(T, size(x))
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
function ($fn){T<:Number}(x::AbstractArray{T}) 
    y = Array(T, size(x))
    return ($fn!)(copy!(y,x))
end
end # eval begin
end #for

end

