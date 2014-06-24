module LiftingTransforms
using ..Util
using ..LiftingSchemes
import ..FilterTransforms: dwt!
export dwt!

# FWT, Forward Wavelet Transform
# IWT, Inverse Wavelet Transform
# DWT, Discrete Wavelet Transform
# periodic, dyadic length (powers of 2)

# inplace transform of y, no vector allocation
function dwt!{T<:FloatingPoint}(y::AbstractVector{T}, L::Integer, scheme::GPLS, fw::Bool, tmp::Vector{T})

    n = length(y)
    J = nscales(n)
    n != 2^J && error("length not a power of 2")
    !(0 <= L <= J) && error("L out of bounds, use 0 <= L <= J")
    L == 0 && return y          # do nothing
    
    if fw
        jrange = (J-1):-1:(J-L)
        stepseq = scheme.step
        ns = n
        half = ns>>1
    else
        jrange = (J-L):(J-1)
        stepseq = reverse(scheme.step)
        ns = int(2^(jrange[1]+1))
        half = ns>>1
    end
    s = y

    for j in jrange
        if fw
            split!(s,tmp,ns)        # lazy transform
        else
            normalize!(s, half, ns, scheme.norm1, scheme.norm2, fw)
        end
        for step in stepseq
            if step.stept == 'p'
                predict!(s, half, step.coef, step.shift, fw)
            elseif step.stept == 'u'
                update!(s, half, step.coef, step.shift, fw)
            else
                error("step type ", step.stept," not supported")
            end
        end       
        if fw
            normalize!(s, half, ns, scheme.norm1, scheme.norm2, fw)
            ns = ns>>1 
            half = half>>1
        else
            merge!(s,tmp,ns)        # reverse split
            ns = ns<<1 
            half = half<<1
        end
    end
    return nothing
end

function normalize!{T<:FloatingPoint}(x::AbstractVector{T}, half::Integer, ns::Integer, n1::Real, n2::Real, fw::Bool)
    if !fw
        n1 = 1/n1
        n2 = 1/n2
    end
    for i = 1:half
        @inbounds x[i] *= n1
    end
    for i = half+1:ns
        @inbounds x[i] *= n2
    end

    return nothing
end

function predict!{T<:FloatingPoint}(x::AbstractVector{T}, half::Integer, c::Vector{T}, shift::Integer, fw::Bool)
    c1=c[1]
    for i = 1:half-1
        @inbounds x[i] -= c1*x[i+half-shift]
    end
    i=half
    x[half] -= c1*x[mod1(i-shift,half)+half]
    return nothing
end

function update!{T<:FloatingPoint}(x::AbstractVector{T}, half::Integer, c::Vector{T}, shift::Integer, fw::Bool)
    rhsis = -shift  # right hand side index shift
    nls = length(c)
    
    irange=2:half
    
    # periodic boundary
    for k = 1:nls
        x[1+half] -= c[k]*x[mod1(k+rhsis,half)]
    end
    # main loop
    if nls==1  # hard code the most common cases (1, 2, 3)
        c1 = c[1]
        for i in irange
            @inbounds x[i+half] -= c1*x[i+rhsis]
        end
    elseif nls==2
        c1,c2 = c[1],c[2]
        rhsisp1 = rhsis+1
        for i in irange
            @inbounds x[i+half] -= c1*x[i+rhsis] + c2*x[i+rhsisp1]
        end
    elseif nls==3
        c1,c2,c3 = c[1],c[2],c[3]
        rhsisp1 = rhsis+1
        rhsisp2 = rhsis+2
        for i = irange
            @inbounds x[i+half] -= c1*x[i+rhsis] + c2*x[i+rhsisp1] + c3*x[i+rhsisp2]
        end
    else
        for i in irange
            for k = 1:nls
                @inbounds x[i+half] -= c[k]*x[i+k-1+rhsis]
            end
        end
    end
    # periodic boundary
    ##
    
    return nothing
end





end
