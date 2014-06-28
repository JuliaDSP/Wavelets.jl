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
    $Xwt{T<:Integer}(x::AbstractVector{T}, scheme::GPLS) = $Xwt(float(x),scheme)
end
end

# 2D MRA (multiresolution analysis)
# Forward Wavelet Transform
# and
# Inverse Wavelet Transform
for (Xwt, fw) in ((:fwt,true),(:iwt,false))
@eval begin
    function $Xwt{T<:FloatingPoint}(x::AbstractMatrix{T}, L::Integer, scheme::GPLS)
        n = size(x,1)
        n != size(x,2) && error("2D dwt: not a square matrix")   
        y = copy(x)
        dwt!(y,L,scheme,$fw)
        return y
    end
    # assume full transform
    function $Xwt{T<:FloatingPoint}(x::AbstractMatrix{T}, scheme::GPLS)
        n = size(x,1)
        n != size(x,2) && error("2D dwt: not a square matrix")   
        y = copy(x)
        dwt!(y,nscales(n),scheme,$fw)
        return y
    end
    $Xwt{T<:Integer}(x::AbstractMatrix{T}, L::Integer, scheme::GPLS) = $Xwt(float(x),L,scheme)
    $Xwt{T<:Integer}(x::AbstractMatrix{T}, scheme::GPLS) = $Xwt(float(x),scheme)
end
end

# 1D
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
            split!(s,tmp,ns)
            for step in stepseq
                if step.stept == 'p'
                    predictfw!(s, half, step.coef, step.shift)
                elseif step.stept == 'u'
                    updatefw!(s, half, step.coef, step.shift)
                else
                    error("step type ", step.stept," not supported")
                end
            end
            normalize!(s, half, ns, scheme.norm1, scheme.norm2)
            ns = ns>>1 
            half = half>>1
        else
            normalize!(s, half, ns, 1/scheme.norm1, 1/scheme.norm2)
            for step in stepseq
                if step.stept == 'p'
                    predictbw!(s, half, step.coef, step.shift)
                elseif step.stept == 'u'
                    updatebw!(s, half, step.coef, step.shift)
                else
                    error("step type ", step.stept," not supported")
                end
            end
            merge!(s,tmp,ns)        # inverse split
            ns = ns<<1 
            half = half<<1
        end
    end
    return y
end

# 2D
# inplace transform of y, no vector allocation
# tmp: size at least n>>2
# tmpvec: size at least n
function dwt!{T<:FloatingPoint}(y::AbstractMatrix{T}, L::Integer, scheme::GPLS, fw::Bool, tmp::Vector{T}=Array(T,size(y,1)>>2), tmpvec::Vector{T}=Array(T,size(y,1)))

    n = size(y,1)
    J = nscales(n)
    n != size(y,2) && error("2D dwt: not a square matrix")   
    n != 2^J && error("length not a power of 2")
    !(0 <= L <= J) && error("L out of bounds, use 0 <= L <= J")
    L == 0 && return y          # do nothing
    
    
    if fw
        jrange = (J-1):-1:(J-L)
        nsub = n
    else
        jrange = (J-L):(J-1)
        nsub = int(2^(J-L+1))
    end
    tmpsub = sub(tmpvec,1:nsub)
    for j in jrange
        
        if fw
            # rows
            xs = n
            for i=1:nsub
                xi = i
                xm = n*(nsub-1)+i
                ya = sub(y, xi:xs:xm)  # final dest and src
                # move to a dense array for speed
                copy!(tmpsub,1,ya,1,nsub)
                dwt!(tmpsub, 1, scheme, fw, tmp)
                copy!(ya,1,tmpsub,1,nsub)
                #dwt!(ya, 1, scheme, fw, tmp)
            end
            # columns
            for i=1:nsub
                xi = 1+(i-1)*n
                xm = xi+nsub-1
                ya = sub(y, xi:xm)
                dwt!(ya, 1, scheme, fw, tmp)
            end       
        else
            # columns
            for i=1:nsub
                xi = 1+(i-1)*n
                xm = xi+nsub-1
                ya = sub(y, xi:xm)
                dwt!(ya, 1, scheme, fw, tmp)
            end   
            # rows
            xs = n
            for i=1:nsub
                xi = i
                xm = n*(nsub-1)+i
                ya = sub(y, xi:xs:xm)  # final dest and src
                # move to a dense array for speed
                copy!(tmpsub,1,ya,1,nsub)
                dwt!(tmpsub, 1, scheme, fw, tmp)
                copy!(ya,1,tmpsub,1,nsub)
                #dwt!(ya, 1, scheme, fw, tmp)
            end

        end 

        fw  && (nsub = nsub>>1)
        !fw && (nsub = nsub<<1)
        fw && (tmpsub = sub(tmpvec,1:nsub))
        !fw && j != jrange[end] && (tmpsub = sub(tmpvec,1:nsub))
        #s = y
    end
    
    return y
end

function normalize!{T<:FloatingPoint}(x::AbstractVector{T}, half::Integer, ns::Integer, n1::Real, n2::Real)
    for i = 1:half
        @inbounds x[i] *= n1
    end
    for i = half+1:ns
        @inbounds x[i] *= n2
    end
    return x
end

# predict and update lifting step inplace on x, forward and backward
# half: half of the length under consideration, shift: shift to left, c: coefs
for (fname,op,puxind,iss,pred) in (  (:predictfw!,:-,:(mod1(i+k-1+rhsis-half,half)+half),:(rhsis = -shift+half; lhsis = 0),true),
                            (:predictbw!,:+,:(mod1(i+k-1+rhsis-half,half)+half),:(rhsis = -shift+half; lhsis = 0),true),
                            (:updatefw!, :-,:(mod1(i+k-1+rhsis,half)),:(rhsis = -shift-half;),false),
                            (:updatebw!, :+,:(mod1(i+k-1+rhsis,half)),:(rhsis = -shift-half;),false)
                            )
@eval begin
function ($fname){T<:FloatingPoint}(x::AbstractVector{T}, half::Integer, c::Vector{T}, shift::Integer)

    nc = length(c)
    # define index shift rhsis
    $iss
    # range limits for 1<=irange<=half without going over boundaries
    irmin = min(max(1, shift+1),  half)
    irmax = max(min(half, half+1+shift-nc),  1)
    irange = irmin:irmax
    # periodic boundary
    lhsr = 1:irmin-1
    if length(irange)==0
        rhsr = irmin:half
    else
        rhsr = irmax+1:half
    end
    if !($pred)  # shift ranges for update
    	irange += half
    	lhsr += half
    	rhsr += half
    end
    # periodic boundary
    for i in lhsr
        for k = 1:nc  
            x[i] = ($op)(x[i], c[k]*x[$puxind] )
        end
    end
    # main loop
    if nc == 1  # hard code the most common cases (1, 2, 3) for speed
        c1 = c[1]
        for i in irange
            x[i] = ($op)(x[i], c1*x[i+rhsis] )
        end
    elseif nc == 2
        c1,c2 = c[1],c[2]
        rhsisp1 = rhsis+1
        for i in irange
            x[i] = ($op)(x[i], c1*x[i+rhsis] + c2*x[i+rhsisp1] )
        end
    elseif nc == 3
        c1,c2,c3 = c[1],c[2],c[3]
        rhsisp1 = rhsis+1
        rhsisp2 = rhsis+2
        for i = irange
            @inbounds x[i] = ($op)(x[i], c1*x[i+rhsis] + c2*x[i+rhsisp1] + c3*x[i+rhsisp2] )
        end
    else
        for i in irange
            for k = 1:nc  
                @inbounds x[i] = ($op)(x[i], c[k]*x[i+k-1+rhsis] )
            end
        end
    end
    # periodic boundary
    for i in rhsr
        for k = 1:nc
            @inbounds x[i] = ($op)(x[i], c[k]*x[$puxind] )
        end
    end

    return x
end
end # eval begin
end # for



end

