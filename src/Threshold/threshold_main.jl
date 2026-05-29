using ..Util, ..WT, ..Transforms

# THRESHOLD TYPES AND FUNCTIONS

abstract type THType end
struct HardTH     <: THType end
struct SoftTH     <: THType end
struct SemiSoftTH <: THType end
struct SteinTH    <: THType end
struct BiggestTH  <: THType end
struct PosTH      <: THType end
struct NegTH      <: THType end

const DEFAULT_TH = HardTH()

# biggest m-term approximation (best m-term approximation for orthogonal transforms)
# result is m-sparse
function threshold!(x::AbstractArray{<:Number}, TH::BiggestTH, m::Int)
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
function threshold!(x::AbstractArray{<:Number}, TH::HardTH, t::Real)
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
function threshold!(x::AbstractArray{<:Number}, TH::SoftTH, t::Real)
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
function threshold!(x::AbstractArray{<:Number}, TH::SemiSoftTH, t::Real)
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
function threshold!(x::AbstractArray{<:Number}, TH::SteinTH, t::Real)
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
function threshold!(x::AbstractArray{<:Number}, TH::NegTH)
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
function threshold!(x::AbstractArray{<:Number}, TH::PosTH)
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
function threshold(x::AbstractArray{T}, TH::THType, t::Real) where T<:Number
    y = Vector{T}(undef, size(x))
    return threshold!(copyto!(y, x), TH, t)
end
function threshold(x::AbstractArray{T}, TH::THType) where T<:Number
    y = Vector{T}(undef, size(x))
    return threshold!(copyto!(y, x), TH)
end


# DENOISING

abstract type DNFT end

struct VisuShrink <: DNFT
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
# estnoise is (x::AbstractArray, wt::Union{DiscreteWavelet,Nothing})
function denoise(x::AbstractArray,
                wt::Union{DiscreteWavelet,Nothing}=DEFAULT_WAVELET;
                L::Int=min(maxtransformlevels(x),6),
                dnt::S=VisuShrink(size(x,1)),
                estnoise::Function=noisest,
                TI::Bool=false,
                nspin::Union{Int,Tuple}=tuple([8 for i=1:ndims(x)]...) ) where S<:DNFT
    iscube(x) || throw(ArgumentError("array must be square/cube"))
    sigma = estnoise(x, wt)

    if TI
        wt == nothing && error("TI not supported with wt=nothing")
        y = zeros(eltype(x), size(x))
        xt = similar(x)
        pns = prod(nspin)

        if ndims(x) == 1
            z = similar(x)
            for i = 1:pns
                shift = i - 1
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
        rmul!(y, 1/pns)
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
function noisest(x::AbstractArray, wt::Union{DiscreteWavelet,Nothing}=DEFAULT_WAVELET, L::Integer = 1)
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

# convert index i to a circshift array starting at 0 shift
nspin2circ(nspin::Int, i::Int) = nspin2circ((nspin,), i)
function nspin2circ(nspin::Tuple, i::Int)
    c1 = CartesianIndices(nspin)[i].I
    c = Vector{Int}(undef, length(c1))
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
        tmp = similar(x, N)
        ftr = similar(x, N)
        aphi = similar(x, length(x))
    end
    spat = zeros(eltype(x), length(y))  # sparse for atom computation
    nmax == -1 && (nmax = length(y))

    while norm(r) > tol && n <= nmax
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
