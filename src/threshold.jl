module Threshold
using ..Util, ..WaveletTypes, ..Transforms
export 
    # threshold
    threshold!,
    threshold,
    HardTH,
    SoftTH,
    SemiSoftTH,
    SteinTH,
    BiggestTH,
    PosTH,
    NegTH,
    # denoising
    DNFT,
    VisuShrink,
    denoise,
    noisest,
    # basis functions
    matchingpursuit,
    bestbasistree,
    # entropy
    coefentropy,
    Entropy,
    ShannonEntropy,
    LogEnergyEntropy


# THRESHOLD TYPES AND FUNCTIONS

abstract THType
immutable HardTH     <: THType end
immutable SoftTH     <: THType end
immutable SemiSoftTH <: THType end
immutable SteinTH    <: THType end
immutable BiggestTH  <: THType end
immutable PosTH      <: THType end
immutable NegTH      <: THType end

const DEF_TH = HardTH()

# biggest m-term approximation (best m-term approximation for orthogonal transforms)
# result is m-sparse
function threshold!{T<:Number}(x::AbstractArray{T}, TH::BiggestTH, m::Int)
    @assert m >= 0
    n = length(x)
    m > n && (m = n)
    ind = sortperm(x, alg=QuickSort, by=abs)
    @inbounds begin
        for i = 1:n-m
            x[ind[i]] = 0
        end
    end
    return x
end

# hard
function threshold!{T<:Number}(x::AbstractArray{T}, TH::HardTH, t::Real)
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
function threshold!{T<:Number}(x::AbstractArray{T}, TH::SoftTH, t::Real)
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
function threshold!{T<:Number}(x::AbstractArray{T}, TH::SemiSoftTH, t::Real)
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
function threshold!{T<:Number}(x::AbstractArray{T}, TH::SteinTH, t::Real)
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

# shrink negative elements to 0
function threshold!{T<:Number}(x::AbstractArray{T}, TH::NegTH)
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
function threshold!{T<:Number}(x::AbstractArray{T}, TH::PosTH)
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
function threshold{T<:Number}(x::AbstractArray{T}, TH::THType, t::Real) 
    y = Array(T, size(x))
    return threshold!(copy!(y,x), TH, t)
end
function threshold{T<:Number}(x::AbstractArray{T}, TH::THType) 
    y = Array(T, size(x))
    return threshold!(copy!(y,x), TH)
end


# DENOISING

abstract DNFT

type VisuShrink <: DNFT
    th::THType      # threshold type
    t::Float64      # threshold for noise level sigma=1, use sigma*t in application
end
# define type for signal length n
function VisuShrink(n::Int)
    return VisuShrink(DEF_TH, sqrt(2*log(n)))
end

const DEFAULT_WAVELET = waveletfilter(WT.sym5)    # default wavelet type

# denoise signal x by thresholding in wavelet space
function denoise{T<:DiscreteWavelet,S<:DNFT}(x::AbstractArray;  
                                    wt::Union(T,Nothing)=DEFAULT_WAVELET, 
                                    level::Int=max(nscales(size(x,1))-6,1),
                                    dnt::S=VisuShrink(size(x,1)),
                                    sigma::Real=noisest(x, wt=wt),
                                    TI::Bool=false,
                                    nspin::Union(Int,Tuple)=tuple([8 for i=1:length(size(x))]...) )
    @assert iscube(x)
    
    if TI
        wt == nothing && error("TI not supported with wt=nothing")
        y = zeros(eltype(x), size(x))
        xt = Array(eltype(x), size(x))
        L = level2tl(size(x,1), level)
        pns = prod(nspin)
        
        if ndims(x) == 1
            z = Array(eltype(x), size(x))
            for i = 1:pns
                shift = nspin2circ(nspin, i)[1]
                circshift!(z, x, shift)
                
                dwt!(xt, z, wt, L, true)
                threshold!(xt, dnt.th, sigma*dnt.t)
                dwt!(z, xt, wt, L, false)
                
                circshift!(xt, z, -shift)
                arrayadd!(y, xt)
            end
        else # ndims > 1
            for i = 1:pns
                shift = nspin2circ(nspin, i)
                z = circshift(x, shift)
                
                dwt!(xt, z, wt, L, true)
                threshold!(xt, dnt.th, sigma*dnt.t)
                dwt!(z, xt, wt, L, false)
                
                z = circshift(z, -shift)
                arrayadd!(y, z)
            end
        end
        scale!(y, 1/pns)
    else # !TI
        if wt == nothing
            y = copy(x)
            threshold!(y, dnt.th, sigma*dnt.t)
        else
            L = level2tl(size(x,1),level)
            y = dwt(x, wt, L)
            threshold!(y, dnt.th, sigma*dnt.t)
            dwt!(y, wt, L, false)
        end
    end
    
    return y
end
# add z to y
function arrayadd!(y::AbstractArray, z::AbstractArray)
    @assert length(y) == length(z)
    for i = 1:length(y)
        @inbounds y[i] += z[i]
    end
    return y
end


# estimate the std. dev. of the signal noise, assuming Gaussian distribution
function noisest{T<:DiscreteWavelet}(x::AbstractArray; wt::Union(T,Nothing)=DEFAULT_WAVELET)
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
nspin2circ(nspin::Int, i::Int) = nspin2circ((nspin,), i)
function nspin2circ(nspin::Tuple, i::Int)
    c1 = ind2sub(nspin,i)
    c = Array(Int, length(c1))
    for k = 1:length(c1)
        c[k] = c1[k]-1
    end
    return c
end


# BASIS FUNCTIONS

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


# ENTROPY MEASURES

abstract Entropy
immutable ShannonEntropy <: Entropy end  #Coifman-Wickerhauser
immutable LogEnergyEntropy <: Entropy end

# Entropy measures: Additive with coefentropy(0) = 0
# all coefs assumed to be on [-1,1] after normalization with nrm
# given x and y, where x has "more concentrated energy" than y
# then coefentropy(x, et, norm) <= coefentropy(y, et, norm) should be satisfied.

function coefentropy{T<:FloatingPoint}(x::T, et::ShannonEntropy, nrm::T)
    s = (x/nrm)^2
    if s == 0.0
        return -zero(T)
    else
        return -s*log(s)
    end
end
function coefentropy{T<:FloatingPoint}(x::T, et::LogEnergyEntropy, nrm::T)
    s = (x/nrm)^2
    if s == 0.0
        return -zero(T)
    else
        return -log(s)
    end
end
function coefentropy{T<:FloatingPoint}(x::AbstractArray{T}, et::Entropy, nrm::T=vecnorm(x))
    @assert nrm >= 0
    sum = zero(T)
    nrm == sum && return sum
    for i in 1:length(x)
        @inbounds sum += coefentropy(x[i], et, nrm)
    end
    return sum
end



function bestbasistree{T<:FloatingPoint}(y::AbstractVector{T}, wt::DiscreteWavelet, L::Integer=nscales(y), et::Entropy=ShannonEntropy())
    bestbasistree(y, wt, maketree(length(y), L, :full), et)
end
function bestbasistree{T<:FloatingPoint}(y::AbstractVector{T}, wt::DiscreteWavelet, tree::BitVector, et::Entropy=ShannonEntropy())

    n = length(y)
    J = nscales(n)
    @assert isdyadic(y)
    @assert isvalidtree(y, tree)
    tree[1] || return besttree      # do nothing
    
    besttree = copy(tree)
    x = copy(y)
    tmp = Array(T, size(y))
    entr = Array(T, length(tree)+length(y))
    nrm = vecnorm(y)
    
    L = J
    k = 1
    while L > 0
        ix = 1
        Lfw = J-L
        nj = detailn(tl2level(n, Lfw))
        
        dtmp = Transforms.unsafe_vectorslice(tmp, 1, nj)
        while ix <= n
            dx = Transforms.unsafe_vectorslice(x, ix, nj)
            
            entr[k] = coefentropy(dx, et, nrm)
            
            dwt!(dtmp, dx, wt, 1, true)
            copy!(dx, dtmp)
            
            ix += nj
            k += 1
        end
        L -= 1
    end
    
    for i = 1:length(y)
        entr[k] = coefentropy(y[i], et, nrm)
        k += 1
    end
    
    # make the best tree
    for i=1:length(tree)
        if (i>1 && !besttree[i>>1]) || !tree[i]  # parent is 0 or input tree-node is 0
            besttree[i] = false
        else
            if entr[i] <= bestsubtree(entr, i)
                besttree[i] = false
            else
                besttree[i] = true
            end
        end
    end
    
    @assert isvalidtree(y, besttree)
    return besttree
end

# the entropy of best subtree
function bestsubtree(entr::Array, i::Int)
    n = length(entr)
    
    if i<<1+1 > n  # bottom of tree
        return entr[i]
    else
        if (i<<1+1)<<1+1 > n  # above bottom of tree
            sum = entr[i<<1]
            sum += entr[i<<1+1]
        else
            sum = bestsubtree(entr, i<<1)
            sum += bestsubtree(entr, i<<1+1)
        end
        
        return min(entr[i], sum)
    end
end



end

