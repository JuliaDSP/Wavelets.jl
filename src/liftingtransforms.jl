module LiftingTransforms
using ..Util
using ..LiftingSchemes
import ..FilterTransforms: fwt,iwt,dwt!
export fwt,iwt,dwt!

# FWT, Forward Wavelet Transform
# IWT, Inverse Wavelet Transform
# DWT, Discrete Wavelet Transform
# periodic boundaries, dyadic length (powers of 2)

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
end
end

# 1D
# inplace transform of y, no vector allocation
# tmp: size at least n>>2
# oopc: out of place computation (e.g. for a non-unit strided vector)
# oopv: the out of place location
function dwt!{T<:FloatingPoint}(y::AbstractVector{T}, L::Integer, scheme::GPLS, fw::Bool, tmp::Vector{T}=Array(T,length(y)>>2); oopc::Bool=false, oopv::Union(AbstractVector{T},Nothing)=nothing)

    n = length(y)
    J = nscales(n)
    n != 2^J && error("length not a power of 2")
    !(0 <= L <= J) && error("L out of bounds, use 0 <= L <= J")
    oopc && oopv == nothing && error("out of place vector not set")
    oopc && n != length(oopv) && error("out of place vector not same length as input")
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
            if oopc && j==jrange[1]
                split!(oopv, y, ns)
                s = oopv
            else
                split!(s, ns, tmp)
            end
            for step in stepseq
                if step.stept == 'p'
                    predictfw!(s, half, convert(Array{T}, step.coef), step.shift)
                elseif step.stept == 'u'
                    updatefw!(s, half, convert(Array{T}, step.coef), step.shift)
                else
                    error("step type ", step.stept," not supported")
                end
            end
            if oopc && L==1  # directly use out of place normalize
                normalize!(y, oopv, half, ns, scheme.norm1, scheme.norm2)
            elseif oopc && j==jrange[end]
                normalize!(s, half, ns, scheme.norm1, scheme.norm2)
                copy!(y, oopv)
            else
                normalize!(s, half, ns, scheme.norm1, scheme.norm2)
            end
            ns = ns>>1 
            half = half>>1
        else
            if oopc && L==1  # directly use out of place normalize
                normalize!(oopv, y, half, ns, 1/scheme.norm1, 1/scheme.norm2)
                s = oopv
            elseif oopc && j==jrange[1]
                copy!(oopv, y)
                s = oopv
                normalize!(s, half, ns, 1/scheme.norm1, 1/scheme.norm2)
            else
                normalize!(s, half, ns, 1/scheme.norm1, 1/scheme.norm2)
            end
            for step in stepseq
                if step.stept == 'p'
                    predictbw!(s, half, convert(Array{T}, step.coef), step.shift)
                elseif step.stept == 'u'
                    updatebw!(s, half, convert(Array{T}, step.coef), step.shift)
                else
                    error("step type ", step.stept," not supported")
                end
            end
            if oopc && j==jrange[end]
                merge!(y, oopv, ns)
            else
                merge!(s, ns, tmp)        # inverse split
            end
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
                # out of place in a dense array for speed
                dwt!(ya, 1, scheme, fw, tmp, oopc=true, oopv=tmpsub)
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
                # out of place in a dense array for speed
                dwt!(ya, 1, scheme, fw, tmp, oopc=true, oopv=tmpsub)
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
	n1 = convert(T, n1)
	n2 = convert(T, n2)
    for i = 1:half
        @inbounds x[i] *= n1
    end
    for i = half+1:ns
        @inbounds x[i] *= n2
    end
    return x
end
# out of place normalize from x to y
function normalize!{T<:FloatingPoint}(y::AbstractVector{T}, x::AbstractVector{T}, half::Integer, ns::Integer, n1::Real, n2::Real)
	n1 = convert(T, n1)
	n2 = convert(T, n2)
    for i = 1:half
        @inbounds y[i] = n1*x[i]
    end
    for i = half+1:ns
        @inbounds y[i] = n2*x[i]
    end
    return y
end

# predict and update lifting step inplace on x, forward and backward
# half: half of the length under consideration, shift: shift to left, c: coefs
for (fname,op,puxind,pred) in (  (:predictfw!,:-,:(mod1(i+k-1+rhsis-half,half)+half),true),
                            (:predictbw!,:+,:(mod1(i+k-1+rhsis-half,half)+half),true),
                            (:updatefw!, :-,:(mod1(i+k-1+rhsis,half)),false),
                            (:updatebw!, :+,:(mod1(i+k-1+rhsis,half)),false)
                            )
@eval begin
function ($fname){T<:FloatingPoint}(x::AbstractVector{T}, half::Integer, c::Vector{T}, shift::Integer)

    nc = length(c)
    # define index shift rhsis
    if $pred
    	rhsis = -shift+half
    else
        rhsis = -shift-half
    end
    # conditions for every element i in irange to be in bounds
    # 1 <= i <= half
    # 1 <= i+1-1-shift <= half
    # 1 <= i+nc-1-shift <= half
    irmin = max(shift+1, 1-nc+shift)
    irmax = min(half+1+shift-nc, half+shift)
    if irmin > half || irmax < 1
        irange = 1:0  # empty
    else
        irmin = max(irmin,1)
        irmax = min(irmax,half)
        irange = irmin:irmax
    end
    # periodic boundary
    if length(irange)==0
        lhsr = 1:half
        rhsr = 1:0
    else
        lhsr = 1:irmin-1
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
            @inbounds x[i] = ($op)(x[i], c[k]*x[$puxind] )
        end
    end
    # main loop
    if nc == 1  # hard code the most common cases (1, 2, 3) for speed
        c1 = c[1]
        for i in irange
            @inbounds x[i] = ($op)(x[i], c1*x[i+rhsis] )
        end
    elseif nc == 2
        c1,c2 = c[1],c[2]
        rhsisp1 = rhsis+1
        for i in irange
            @inbounds x[i] = ($op)(x[i], c1*x[i+rhsis] + c2*x[i+rhsisp1] )
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

