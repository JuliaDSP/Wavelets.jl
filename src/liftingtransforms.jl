module LiftingTransforms
using ..Util
using ..LiftingSchemes
import ..FilterTransforms: fwt,iwt,dwt!
export dwt!

# FWT, Forward Wavelet Transform
# IWT, Inverse Wavelet Transform
# DWT, Discrete Wavelet Transform
# periodic, dyadic length (powers of 2)

# 1D
# Forward Wavelet Transform
# and
# Inverse Wavelet Transform
for (Xwt, fw) in ((:fwt,true),(:iwt,false))
@eval begin
    function $Xwt{T<:FloatingPoint}(x::AbstractVector{T}, L::Integer, scheme::GPLS)
        y = copy(x)
        tmp = Array(T,length(x)>>2)
        dwt!(y,L,scheme,$fw,tmp)
        return y
    end
    # assume full transform
    function $Xwt{T<:FloatingPoint}(x::AbstractVector{T}, scheme::GPLS)
        y = copy(x)
        tmp = Array(T,length(x)>>2)
        dwt!(y,nscales(length(x)),scheme,$fw,tmp)
        return y
    end
    $Xwt{T<:Integer}(x::AbstractVector{T}, L::Integer, scheme::GPLS) = $Xwt(float(x),L,scheme)
    $Xwt{T<:Integer}(x::AbstractVector{T}, scheme::GPLS) = $Xwt(float(x),nscales(length(x)),scheme)
end
end

# inplace transform of y, no vector allocation
# tmp: size at least n>>2
function dwt!{T<:FloatingPoint}(y::AbstractVector{T}, L::Integer, scheme::GPLS, fw::Bool, tmp::Vector{T}=Array(T,length(y)>>2))

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
        ns = 2^(jrange[1]+1)
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
                predictupdate!(s, half, step.coef, step.shift, fw, predict=true)
            elseif step.stept == 'u'
                predictupdate!(s, half, step.coef, step.shift, fw, predict=false)
            else
                error("step type ", step.stept," not supported")
            end
        end
        if fw
            normalize!(s, half, ns, scheme.norm1, scheme.norm2, fw)
            ns = ns>>1 
            half = half>>1
        else
            merge!(s,tmp,ns)        # inverse split
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

# predict or update lifting step inplace on x
# half: half of the length under consideration, shift: shift to left, c: coefs
function predictupdate!{T<:FloatingPoint}(x::AbstractVector{T}, half::Integer, c::Vector{T}, shift::Integer, fw::Bool; predict::Bool=true)

    nc = length(c)
    if predict
        rhsis = -shift+half  
        lhsis = 0
    else
        rhsis = -shift      # right hand side index shift
        lhsis = half        # left hand side index shift
    end
    # range limits for 1<=irange<=half without going over boundaries
    irmin = min(max(1, shift+1),  half)
    irmax = max(min(half, half+1+shift-nc),  1)
    irange = irmin:irmax
    # periodic boundary
    if length(irange)==0
        rhsr = irmin:half
    else
        rhsr = irmax+1:half
    end
    # the only difference between fw and !fw is the += and -= symbols
    # for predict=true, the modulus has to be shifted in the update case
    if fw
        # periodic boundary
        if predict
            for i in 1:irmin-1
                for k = 1:nc  
                    @inbounds x[i+lhsis] -= c[k]*x[mod1(i+k-1+rhsis-half,half)+half] 
                end
            end
        else
            for i in 1:irmin-1
                for k = 1:nc  
                    @inbounds x[i+lhsis] -= c[k]*x[mod1(i+k-1+rhsis,half)] 
                end
            end
        end
        # main loop
        if nc == 1  # hard code the most common cases (1, 2, 3) for speed
            c1 = c[1]
            for i in irange
                @inbounds x[i+lhsis] -= c1*x[i+rhsis] 
            end
        elseif nc == 2
            c1,c2 = c[1],c[2]
            rhsisp1 = rhsis+1
            for i in irange
                @inbounds x[i+lhsis] -= c1*x[i+rhsis] + c2*x[i+rhsisp1] 
            end
        elseif nc == 3
            c1,c2,c3 = c[1],c[2],c[3]
            rhsisp1 = rhsis+1
            rhsisp2 = rhsis+2
            for i = irange
                @inbounds x[i+lhsis] -= c1*x[i+rhsis] + c2*x[i+rhsisp1] + c3*x[i+rhsisp2] 
            end
        else
            for i in irange
                for k = 1:nc  
                    @inbounds x[i+lhsis] -= c[k]*x[i+k-1+rhsis] 
                end
            end
        end
        # periodic boundary
        if predict
            for i in rhsr
                for k = 1:nc
                    @inbounds x[i+lhsis] -= c[k]*x[mod1(i+k-1+rhsis-half,half)+half] 
                end
            end
        else
            for i in rhsr
                for k = 1:nc
                    @inbounds x[i+lhsis] -= c[k]*x[mod1(i+k-1+rhsis,half)] 
                end
            end
        end
        
    else  # !fw
        # periodic boundary
        if predict
            for i in 1:irmin-1
                for k = 1:nc  
                    @inbounds x[i+lhsis] += c[k]*x[mod1(i+k-1+rhsis-half,half)+half] 
                end
            end
        else
            for i in 1:irmin-1
                for k = 1:nc  
                    @inbounds x[i+lhsis] += c[k]*x[mod1(i+k-1+rhsis,half)] 
                end
            end
        end
        # main loop
        if nc == 1  # hard code the most common cases (1, 2, 3) for speed
            c1 = c[1]
            for i in irange
                @inbounds x[i+lhsis] += c1*x[i+rhsis] 
            end
        elseif nc == 2
            c1,c2 = c[1],c[2]
            rhsisp1 = rhsis+1
            for i in irange
                @inbounds x[i+lhsis] += c1*x[i+rhsis] + c2*x[i+rhsisp1] 
            end
        elseif nc == 3
            c1,c2,c3 = c[1],c[2],c[3]
            rhsisp1 = rhsis+1
            rhsisp2 = rhsis+2
            for i = irange
                @inbounds x[i+lhsis] += c1*x[i+rhsis] + c2*x[i+rhsisp1] + c3*x[i+rhsisp2] 
            end
        else
            for i in irange
                for k = 1:nc  
                    @inbounds x[i+lhsis] += c[k]*x[i+k-1+rhsis] 
                end
            end
        end
        # periodic boundary
        if predict
            for i in rhsr
                for k = 1:nc
                    @inbounds x[i+lhsis] += c[k]*x[mod1(i+k-1+rhsis-half,half)+half] 
                end
            end
        else
            for i in rhsr
                for k = 1:nc
                    @inbounds x[i+lhsis] += c[k]*x[mod1(i+k-1+rhsis,half)] 
                end
            end
        end
        
    end


    return nothing
end





end
