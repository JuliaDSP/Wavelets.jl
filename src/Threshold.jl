module Threshold
using ..Util, ..WT, ..Transforms
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

const DEFAULT_TH = HardTH()

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
        for i in eachindex(x)
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
        for i in eachindex(x)
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
        for i in eachindex(x)
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
        for i in eachindex(x)
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
        for i in eachindex(x)
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
        for i in eachindex(x)
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
    VisuShrink(th, t) = new(th, t)
end
# define type for signal length n
function VisuShrink(n::Int)
    return VisuShrink(DEFAULT_TH, sqrt(2*log(n)))
end

const DEFAULT_WAVELET = wavelet(WT.sym5, WT.Filter)    # default wavelet type

# denoise signal x by thresholding in wavelet space
# estnoise is (x::AbstractArray, wt::Union{DiscreteWavelet,Void})
function denoise{S<:DNFT}(x::AbstractArray,
                        wt::Union{DiscreteWavelet,Void}=DEFAULT_WAVELET;
                        L::Int=min(maxtransformlevels(x),6),
                        dnt::S=VisuShrink(size(x,1)),
                        estnoise::Function=noisest,
                        TI::Bool=false,
                        nspin::Union{Int,Tuple}=tuple([8 for i=1:ndims(x)]...) )
    iscube(x) || throw(ArgumentError("array must be square/cube"))
    sigma = estnoise(x, wt)

    if TI
        wt == nothing && error("TI not supported with wt=nothing")
        y = zeros(eltype(x), size(x))
        xt = Array(eltype(x), size(x))
        pns = prod(nspin)

        if ndims(x) == 1
            z = Array(eltype(x), size(x))
            for i = 1:pns
                shift = nspin2circ(nspin, i)[1]
                Util.circshift!(z, x, shift)

                Transforms.dwt_oop!(xt, z, wt, L)
                threshold!(xt, dnt.th, sigma*dnt.t)
                Transforms.idwt_oop!(z, xt, wt, L)

                Util.circshift!(xt, z, -shift)
                arrayadd!(y, xt)
            end
        else # ndims > 1
            for i = 1:pns
                shift = nspin2circ(nspin, i)
                z = circshift(x, shift)

                Transforms.dwt_oop!(xt, z, wt, L)
                threshold!(xt, dnt.th, sigma*dnt.t)
                Transforms.idwt_oop!(z, xt, wt, L)

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
            y = dwt(x, wt, L)
            threshold!(y, dnt.th, sigma*dnt.t)
            if isa(wt, GLS)
                idwt!(y, wt, L)
            else
                y2 = idwt(y, wt, L)
                y = y2
            end
        end
    end

    return y
end
# add z to y
function arrayadd!(y::AbstractArray, z::AbstractArray)
    length(y) == length(z) || throw(DimensionMismatch("lengths must be equal"))
    for i in eachindex(y)
        @inbounds y[i] += z[i]
    end
    return y
end


# estimate the std. dev. of the signal noise, assuming Gaussian distribution
function noisest(x::AbstractArray, wt::Union{DiscreteWavelet,Void}=DEFAULT_WAVELET, L::Integer = 1)
    if wt == nothing
        y = x
    else
        y = dwt(x, wt, L)
    end
    dr = y[detailrange(y, L)]
    return mad!(dr)/0.6745
end
# Median absolute deviation
function mad!(y::AbstractArray)
    m = median!(y)
    for i in eachindex(y)
        y[i] = abs(y[i]-m)
    end
    return median!(y)
end
#function mad(x::AbstractArray)
#    y = copy(x)
#    mad!(y)
#end

# convert index i to a circshift array starting at 0 shift
nspin2circ(nspin::Int, i::Int) = nspin2circ((nspin,), i)
function nspin2circ(nspin::Tuple, i::Int)
    c1 = ind2sub(nspin,i)
    c = Array(Int, length(c1))
    for k in 1:length(c1)
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
    @inbounds for i in eachindex(x)
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

function coefentropy{T<:AbstractFloat}(x::T, et::ShannonEntropy, nrm::T)
    s = (x/nrm)^2
    if s == 0.0
        return -zero(T)
    else
        return -s*log(s)
    end
end
function coefentropy{T<:AbstractFloat}(x::T, et::LogEnergyEntropy, nrm::T)
    s = (x/nrm)^2
    if s == 0.0
        return -zero(T)
    else
        return -log(s)
    end
end
function coefentropy{T<:AbstractFloat}(x::AbstractArray{T}, et::Entropy, nrm::T=vecnorm(x))
    @assert nrm >= 0
    sum = zero(T)
    nrm == sum && return sum
    for i in eachindex(x)
        @inbounds sum += coefentropy(x[i], et, nrm)
    end
    return sum
end


# find the best tree that is a subset of the input tree (use :full to find the best tree)
# for wpt
function bestbasistree{T<:AbstractFloat}(y::AbstractVector{T}, wt::DiscreteWavelet, L::Integer=maxtransformlevels(y), et::Entropy=ShannonEntropy())
    bestbasistree(y, wt, maketree(length(y), L, :full), et)
end
function bestbasistree{T<:AbstractFloat}(y::AbstractVector{T}, wt::DiscreteWavelet, tree::BitVector, et::Entropy=ShannonEntropy())

    isvalidtree(y, tree) || throw(ArgumentError("invalid tree"))

    tree[1] || copy(tree)      # do nothing

    x = copy(y)
    n = length(y)
    tmp = Array(T, n)
    ntree = length(tree)
    entr_bf = Array(T, ntree)
    nrm = vecnorm(y)

    Lmax = maxtransformlevels(n)
    L = Lmax
    k = 1
    while L > 0
        ix = 1
        Lfw = Lmax-L
        nj = detailn(n, Lfw)

        @assert nj <= n
        dtmp = Transforms.unsafe_vectorslice(tmp, 1, nj)
        while ix <= n
            @assert nj+ix-1 <= n
            dx = Transforms.unsafe_vectorslice(x, ix, nj)

            entr_bf[k] = coefentropy(dx, et, nrm)

            dwt!(dtmp, dx, wt, 1)
            copy!(dx, dtmp)

            ix += nj
            k += 1
        end
        L -= 1
    end

    # entropy of fully transformed signal (end nodes)
    n_af = 2^(Lmax-1)
    entr_af = Array(T, n_af)
    n_coef_af = div(n, n_af)
    for i in 1:n_af
        range = (i-1)*n_coef_af+1 : i*n_coef_af
        entr_af[i] = coefentropy(x[range], et, nrm)
    end

    # make the best tree
    besttree = copy(tree)
    for i in 1:ntree
        if (i>1 && !besttree[i>>1]) || !tree[i]  # parent is 0 or input tree-node is 0
            besttree[i] = false
        else
            if entr_bf[i] <= bestsubtree_entropy(entr_bf, entr_af, i)
                besttree[i] = false
            else
                besttree[i] = true
            end
        end
    end

    @assert isvalidtree(y, besttree)
    return besttree::BitVector
end

# the entropy of best subtree
function bestsubtree_entropy(entr_bf::Array, entr_af::Array, i::Int)
    n = length(entr_bf)
    n_af = length(entr_af)
    @assert isdyadic(n+1)
    @assert isdyadic(n_af)
    @assert n + 1 == n_af<<1

    if n < (i<<1)  # bottom of tree
        sum  = entr_af[i - n_af + 1]
    else
        sum  = bestsubtree_entropy(entr_bf, entr_af, i<<1)
        sum += bestsubtree_entropy(entr_bf, entr_af, i<<1+1)
    end

    lowestentropy = min(entr_bf[i], sum)
    return lowestentropy
end


end # module
