module LiftingTransforms
using ..Util
using ..LiftingSchemes
export dwtlift!

# FWT, Forward Wavelet Transform
# IWT, Inverse Wavelet Transform
# DWT, Discrete Wavelet Transform
# dyadic length (powers of 2)

# inplace transform of y, no vector allocation
function dwtlift!{T<:FloatingPoint}(y::AbstractVector{T}, L::Integer, scheme::GLS, fw::Bool, tmp::Vector{T})

    n = length(y)
    J = nscales(n)
    n != 2^J && error("length not a power of 2")
    !(0 <= L <= J) && error("L out of bounds, use 0 <= L <= J")
    L == 0 && return y    # do nothing
    
    if fw
        jrange = (J-1):-1:(J-L)
    else
        jrange = (J-L):(J-1)
    end
    s = y
    ns = n
    half = ns>>1
    for j in jrange
        # lazy transform
        split!(s,tmp,ns)
        
        for step in scheme.step
            if step.stept == 'p'
                predict!(s, half, step.coef, step.shift, fw)
            elseif step.stept == 'u'
                update!(s, half, step.coef, step.shift, fw)
            else
                error("step type ", step.stept," not supported")
            end

        end
        # normalize
        for i = 1:half
            @inbounds s[i] *= n1
        end
        for i = half+1:ns
            @inbounds s[i] *= n2
        end

        ns = ns>>1
        half = half>>1
    end
    return y
end

function update!{T<:FloatingPoint}(x::AbstractVector{T}, half::Integer, c::Vector{T}, shift::Integer, fw::Bool)
    rhsis = -shift  # right hand side index shift
    nls=length(c)
    
    irange=2:half
    
    # periodic boundary
    for k = 1:nls
        x[1+half] -= c[k]*x[mod1(k+rhsis,half)]
    end
    # main loop
    if nls==1  # hard code the most common cases (1, 2, 3)
        c1 = c[1]
        for i = irange
            @inbounds x[i+half] -= c1*x[i+rhsis]
        end
    elseif nls==2
        ls1,ls2 = ls[1],ls[2]
        for i = irange
            @inbounds x[i+half] -= c1*x[i+rhsis] + c2*x[i+1+rhsis]
        end
    elseif nls==3
        c1,c2,c3 = c[1],c[2],c[3]
        for i = irange
            @inbounds x[i+half] -= c1*x[i+rhsis] + c2*x[i+1+rhsis] + c3*x[i+2+rhsis]
        end
    else
        for i = irange
            for k = 1:nls
                @inbounds x[i+half] -= c[k]*x[i+k-1+rhsis]
            end
        end
    end
    # periodic boundary
    ##
    
    return x
end





end
