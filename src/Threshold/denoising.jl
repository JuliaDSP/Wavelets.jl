using ..Transforms



# DENOISING

abstract type DNFT end

struct VisuShrink <: DNFT
    th::THType      # threshold type
    t::Float64      # threshold for noise level sigma=1, use sigma*t in application
end
# define type for signal length n
function VisuShrink(n::Int)
    return VisuShrink(DEFAULT_TH, sqrt(2 * log(n)))
end

const DEFAULT_WAVELET = wavelet(WT.sym5, WT.Filter)    # default wavelet type

# denoise signal x by thresholding in wavelet space
# estnoise is (x::AbstractArray, wt::Union{DiscreteWavelet,Nothing})
function denoise(
        x::AbstractArray,
        wt::Union{DiscreteWavelet,Nothing}=DEFAULT_WAVELET;
        L::Int=min(maxtransformlevels(x), 6),
        dnt::S=VisuShrink(size(x, 1)),
        estnoise::Function=noisest,
        TI::Bool=false,
        nspin::Union{Int,Tuple}=ntuple(Returns(8), ndims(x))
    ) where {S<:DNFT}
    iscube(x) || throw(ArgumentError("array must be square/cube"))
    sigma = estnoise(x, wt)

    if TI
        isnothing(wt) && error("TI not supported with wt=nothing")
        y = zeros(eltype(x), size(x))
        xt = similar(x)
        pns = prod(nspin)

        if ndims(x) == 1
            z = similar(x)
            for i = 1:pns
                shift = i - 1
                circshift!(z, x, shift)

                Transforms.dwt_oop!(xt, z, wt, L)
                threshold!(xt, dnt.th, sigma * dnt.t)
                Transforms.idwt_oop!(z, xt, wt, L)

                circshift!(xt, z, -shift)
                arrayadd!(y, xt)
            end
        else # ndims > 1
            for i = 1:pns
                shift = nspin2circ(nspin, i)
                z = circshift(x, shift)

                Transforms.dwt_oop!(xt, z, wt, L)
                threshold!(xt, dnt.th, sigma * dnt.t)
                Transforms.idwt_oop!(z, xt, wt, L)

                z = circshift(z, -shift)
                arrayadd!(y, z)
            end
        end
        rmul!(y, 1 / pns)
    else # !TI
        if isnothing(wt)
            y = copy(x)
            threshold!(y, dnt.th, sigma * dnt.t)
        else
            y = dwt(x, wt, L)
            threshold!(y, dnt.th, sigma * dnt.t)
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
function noisest(x::AbstractArray, wt::Union{DiscreteWavelet,Nothing}=DEFAULT_WAVELET, L::Integer=1)
    if isnothing(wt)
        y = x
    else
        y = dwt(x, wt, L)
    end
    dr = y[detailrange(y, L)]
    return mad!(dr) / 0.6745
end
# Median absolute deviation
function mad!(y::AbstractArray)
    m = median!(y)
    for i in eachindex(y)
        y[i] = abs(y[i] - m)
    end
    return median!(y)
end

# convert index i to a circshift array starting at 0 shift
nspin2circ(nspin::Int, i::Int) = nspin2circ((nspin,), i)
function nspin2circ(nspin::Tuple, i::Int)
    c1 = CartesianIndices(nspin)[i].I
    c = Vector{Int}(undef, length(c1))
    for k in eachindex(c1)
        c[k] = c1[k] - 1
    end
    return c
end
