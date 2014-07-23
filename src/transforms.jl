module Transforms
using ..Util, ..WaveletTypes
export fwt, iwt, dwt!, fwtc, iwtc

# includes:
#  GENERAL TRANSFORM FUNCTIONS
#  PO FILTER TRANSFORMS 
#  LIFTING TRANSFORMS 

# FWT, Forward Wavelet Transform
# IWT, Inverse Wavelet Transform
# DWT, Discrete Wavelet Transform
# FWTC, Forward Wavelet Transform, Column-wise
# IWTC, Inverse Wavelet Transform, Column-wise

##################################################################################
#
#  GENERAL TRANSFORM FUNCTIONS
#
##################################################################################

# general functions for all types and dimensions
for Xwt in (:fwt, :iwt, :fwtc, :iwtc)
@eval begin
	# assume full transform
	$Xwt(x::AbstractArray, wt::WaveletType) = $Xwt(x,nscales(size(x,1)),wt)
	# int -> float
	$Xwt{T<:Integer}(x::AbstractArray{T}, L::Integer,  wt::WaveletType) = $Xwt(float(x),L,wt)
	# swap arguments
	$Xwt(x::AbstractArray, wt::WaveletType, L::Integer) = $Xwt(x, L, wt)
end
end

# 1-d
# 2-d MRA (multiresolution analysis)
# Forward Wavelet Transform
# and
# Inverse Wavelet Transform
for (Xwt, fw) in ((:fwt,true),(:iwt,false))
@eval begin
    function $Xwt{T<:FloatingPoint}(x::AbstractArray{T}, L::Integer, wt::WaveletType)
        if typeof(wt) == POfilter
            y = Array(T, size(x))
            dwt!(y, x, L, wt, $fw)
        elseif typeof(wt) == GPLS
            y = copy(x)
            dwt!(y, L, wt, $fw)
        else
            error("unknown wavelet type")
        end
        return y
    end
end
end

# column-wise transforms or color images, transform each x[:,...,:,i] separately
for Xwt in (:fwtc, :iwtc)
@eval begin
    function $Xwt{T<:FloatingPoint}(x::AbstractArray{T}, L::Integer, wt::WaveletType)
        dim = ndims(x)
        cn = size(x, dim)
        y = Array(eltype(x), size(x))
        
        ind = Array(Any, dim)
        for i = 1:dim
            ind[i] = 1:size(x, i)
        end
        
        for d = 1:cn
            ind[dim] = d
	        xc = reshape(x[ind...], size(x)[1:end-1]...)
	        y[ind...] = fwt(xc, L, wt)
        end
        return y
    end
end
end



##################################################################################
#
#  PO FILTER TRANSFORMS 
#  Periodic boundaries, Orthogonal, dyadic length (powers of 2)
#
##################################################################################

function makeqmf(h::AbstractVector, fw::Bool, T::Type=eltype(h))
	scfilter, dcfilter = makereverseqmf(h, fw, T)
    return reverse(scfilter), reverse(dcfilter)
end
function makereverseqmf(h::AbstractVector, fw::Bool, T::Type=eltype(h))
	h = convert(Vector{T}, h)
    if fw
    	scfilter = reverse(h)
        dcfilter = mirror(h)
    else
        scfilter = h
        dcfilter = reverse(mirror(h))
    end
    return scfilter, dcfilter
end

# 1-d
# writes to y
function dwt!{T<:FloatingPoint}(y::AbstractVector{T}, x::AbstractVector{T}, L::Integer, filter::OrthoFilter, fw::Bool)
    si = Array(T, filter.n-1)       # tmp filter vector
    scfilter, dcfilter = makereverseqmf(filter.qmf, fw, T)
    
    dwt!(y, x, L, filter, fw, dcfilter, scfilter, si)
    return y
end
# pseudo "inplace" by copying
function dwt!{T<:FloatingPoint}(x::AbstractArray{T}, L::Integer, filter::OrthoFilter, fw::Bool)
    y = Array(T, size(x))
    dwt!(y, x, L, filter, fw)
    copy!(x,y)
    return x
end
function dwt!{T<:FloatingPoint}(y::AbstractVector{T}, x::AbstractVector{T}, L::Integer, filter::POfilter, fw::Bool, dcfilter::Vector{T}, scfilter::Vector{T}, si::Vector{T}, snew::Vector{T} = Array(T, ifelse(L>1, length(x)>>1, 0)))
    n = length(x)
    J = nscales(n)
    @assert size(x) == size(y)
    @assert isdyadic(y)
    @assert 0 <= L <= J
    is(y,x) && error("input vector is output vector")
    
    L == 0 && return copy!(y,x)     # do nothing
    L > 1 && length(snew) != n>>1 && error("length of snew incorrect")
    s = x                           # s is currect scaling coefs location
    
    if fw
        jrange = (J-1):-1:(J-L)
    else
        jrange = (J-L):(J-1)
    end
    for j in jrange
        if fw
            # detail coefficients
            filtdown!(dcfilter, si, y, detailindex(j,1), detailn(j), s, 1,-filter.n+1,1)
            # scaling coefficients
            filtdown!(scfilter, si, y, 1,                detailn(j), s, 1, 0, 0)
        else
            # scaling coefficients
            filtup!(false,scfilter,  si, y, 1, detailn(j+1), s, 1, -filter.n+1, 0)
            # detail coefficients
            filtup!(true,dcfilter, si, y, 1, detailn(j+1), x, detailindex(j,1), 0, 1)
        end
        # if not final iteration: copy to tmp location
        fw  && j != jrange[end] && copy!(snew,1,y,1,detailn(j))
        !fw && j != jrange[end] && copy!(snew,1,y,1,detailn(j+1))
        L > 1 && (s = snew)
    end
    return y
end

# 2-d
# writes to y
function dwt!{T<:FloatingPoint}(y::AbstractMatrix{T}, x::AbstractMatrix{T}, L::Integer, filter::OrthoFilter, fw::Bool)
    n = size(x,1)
    si = Array(T, filter.n-1)       # tmp filter vector
    tmpvec = Array(T,n)             # tmp storage vector
    scfilter, dcfilter = makereverseqmf(filter.qmf, fw, T)
    
    dwt!(y, x, L, filter, fw, dcfilter, scfilter, si, tmpvec)
    return y
end
function dwt!{T<:FloatingPoint}(y::AbstractMatrix{T}, x::AbstractMatrix{T}, L::Integer, filter::OrthoFilter, fw::Bool, dcfilter::Vector{T}, scfilter::Vector{T}, si::Vector{T}, tmpvec::Vector{T})

    n = size(x,1)
    J = nscales(n)
    @assert size(x) == size(y)
    @assert iscube(y)
    @assert isdyadic(y)
    @assert 0 <= L <= J
    is(y,x) && error("input matrix is output matrix")
    
    L == 0 && return copy!(y,x)        # do nothing
    #s = x

    if fw
        jrange = (J-1):-1:(J-L)
        nsub = n
    else
        jrange = (J-L):(J-1)
        nsub = int(2^(J-L+1))
        copy!(y,x) # !!! not needed if we have an iwt with seperate arrays for s and d 
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
                if j != jrange[1]
                    copy!(tmpsub,1,ya,1,nsub)
                else  # use x in first iteration
                    xa = sub(x, xi:xs:xm)  # src
                    copy!(tmpsub,1,xa,1,nsub)
                end
                dwt!(ya, tmpsub, 1, filter, fw, dcfilter, scfilter, si)
            end
            # columns
            for i=1:nsub
                xi = 1+(i-1)*n
                xm = xi+nsub-1
                ya = sub(y, xi:xm)  # final dest and src
                copy!(tmpsub,1,ya,1,nsub)
                dwt!(ya, tmpsub, 1, filter, fw, dcfilter, scfilter, si)
            end       
        else
            # columns
            for i=1:nsub
                xi = 1+(i-1)*n
                xm = xi+nsub-1
                ya = sub(y, xi:xm)  # final dest and src
                copy!(tmpsub,1,ya,1,nsub)  # x has been copied to y
                dwt!(ya, tmpsub, 1, filter, fw, dcfilter, scfilter, si)
            end   
            # rows
            xs = n
            for i=1:nsub
                xi = i
                xm = n*(nsub-1)+i
                ya = sub(y, xi:xs:xm)  # final dest and src
                copy!(tmpsub,1,ya,1,nsub)
                dwt!(ya, tmpsub, 1, filter, fw, dcfilter, scfilter, si)
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


# periodic filter and downsampling (by 2)
# b : filter
# si: tmp array of length 1 less than b
# out : result gets written to out[iout:nout-1], where nout==nx/2
# x : filter convolved with x[ix:nx-1] (shifted by shift)
# sw : shift downsampling (0 or 1)
# based on Base.filt
function filtdown!{T<:FloatingPoint}(b::Vector{T}, si::Vector{T},
                              out::Union(AbstractVector{T}, T),  iout::Integer, nout::Integer, 
                              x::AbstractVector{T}, ix::Integer, shift::Integer=0, sw::Integer=0)
    nx = nout<<1
    silen = length(si)
    bs = length(b)
    @assert length(x) >= ix+nx-1
    @assert length(out) >= iout+nout-1
    @assert (nx-1)>>1 + iout <= length(out)
    @assert sw==0 || sw==1
    @assert silen == bs-1

    fill!(si,0.0)
    istart = bs + sw
    dsshift = (bs%2 + sw)%2  # is bs odd, and shift downsampling
    @inbounds begin
        for i = 1:(nx-1+istart)
            # periodic in the range [ix:ix+nx-1]
            xindex = mod(i-1+shift,nx) + ix
            xatind = x[xindex]
            # take every other value after istart
            if (i+dsshift)%2 == 0 && i >= istart
                out[(i-istart)>>1 + iout] = si[1] + b[1]*xatind
            end
            for j=1:(silen-1)
                si[j] = si[j+1] + b[j+1]*xatind
            end
            si[silen] = b[silen+1]*xatind
        end
    end

    return nothing
end

# periodic filter and upsampling (by 2)
# b : filter
# si: tmp array of length 1 less than b
# out : result gets written to out[iout:nout-1], where nout==nx*2
# x : filter convolved with x[ix:nx-1] upsampled (then shifted by shift)
# sw : shift upsampling (0 or 1)
# based on Base.filt
function filtup!{T<:FloatingPoint}(add2out::Bool,b::Vector{T}, si::Vector{T},
                              out::Union(AbstractVector{T}, T),  iout::Integer, nout::Integer, 
                              x::AbstractVector{T}, ix::Integer, shift::Integer=0, sw::Integer=0)
    nx = nout>>1
    silen = length(si)
    bs = length(b)
    @assert length(x) >= ix+nx-1                # check array size
    @assert length(out) >= iout+nout-1          # check array size
    @assert (nx<<1-1)>>1 + iout <= length(out)  # max for out index
    @assert sw==0 || sw==1
    @assert silen == bs-1
  
    fill!(si,0.0)
    istart = bs - shift%2
    dsshift = (sw)%2  # shift upsampling
    @inbounds begin
        for i = 1:(nout-1+istart)
            # periodic in the range [ix:ix+nx-1]
            if (i+dsshift)%2==0
                xatind = 0.0
            else
                xindex = mod((i-1)>>1+shift>>1,nx) + ix    #(i-1)>>1 increm. every other
                xatind = x[xindex]
            end
            if i >= istart
                if add2out
                    out[(i-istart) + iout] += si[1] + b[1]*xatind
                else
                    out[(i-istart) + iout] = si[1] + b[1]*xatind
                end
            end
            if (i+dsshift)%2==0
                for j=1:(silen-1)
                    si[j] = si[j+1]  # + 0
                end
            else
                for j=1:(silen-1)
                    si[j] = si[j+1] + b[j+1]*xatind
                end
            end
            si[silen] = b[silen+1]*xatind
        end
    end

    return nothing
end


##################################################################################
#
#  LIFTING TRANSFORMS 
#  Periodic boundaries, dyadic length (powers of 2)
#
##################################################################################

# 1-d
# inplace transform of y, no vector allocation
# tmp: size at least n>>2
# oopc: out of place computation (e.g. for a non-unit strided vector)
# oopv: the out of place location
function dwt!{T<:FloatingPoint}(y::AbstractVector{T}, L::Integer, scheme::GPLS, fw::Bool, tmp::Vector{T}=Array(T,length(y)>>2); oopc::Bool=false, oopv::Union(AbstractVector{T},Nothing)=nothing)

    n = length(y)
    J = nscales(n)
    @assert isdyadic(y)
    @assert 0 <= L <= J
    @assert !(oopc && oopv == nothing)
    @assert !(oopc && n != length(oopv))
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
# pseudo "out of place" by copying
function dwt!{T<:FloatingPoint}(y::AbstractArray{T}, x::AbstractArray{T}, L::Integer, scheme::NRLS, fw::Bool)
    copy!(y, x)
    dwt!(y, L, scheme, fw)
    return y
end

# 2-d
# inplace transform of y, no vector allocation
# tmp: size at least n>>2
# tmpvec: size at least n
function dwt!{T<:FloatingPoint}(y::AbstractMatrix{T}, L::Integer, scheme::NRLS, fw::Bool, tmp::Vector{T}=Array(T,size(y,1)>>2), tmpvec::Vector{T}=Array(T,size(y,1)))

    n = size(y,1)
    J = nscales(n)
    @assert iscube(y)
    @assert isdyadic(y)
    @assert 0 <= L <= J
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

