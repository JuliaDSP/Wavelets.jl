module POtransforms
using ..Util
using ..POfilters
export fwt, iwt, dwt!


# FWT, Forward Wavelet Transform
# IWT, Inverse Wavelet Transform
# DWT, Discrete Wavelet Transform
# periodic boundaries, orthogonal, dyadic length (powers of 2)
# L=0 does nothing, L=1 transforms 1 level, L=J=length(x) is a full transform

# 1D
# Forward Wavelet Transform
# and
# Inverse Wavelet Transform
for (Xwt, fw) in ((:fwt,true),(:iwt,false))
@eval begin
    function $Xwt{T<:FloatingPoint}(x::AbstractVector{T}, L::Integer, filter::POfilter)
        y = Array(T,length(x))
        dwt!(y,x,L,filter,$fw)
        return y
    end
    # assume full transform
    function $Xwt{T<:FloatingPoint}(x::AbstractVector{T}, filter::POfilter)
        y = Array(T,length(x))
        dwt!(y,x,nscales(length(x)),filter,$fw)
        return y
    end
    $Xwt{T<:Integer}(x::AbstractVector{T}, L::Integer, filter::POfilter) = $Xwt(float(x),L,filter)
    $Xwt{T<:Integer}(x::AbstractVector{T}, filter::POfilter) = $Xwt(float(x),nscales(length(x)),filter)
end
end

# 2D MRA (multiresolution analysis)
# Forward Wavelet Transform
# and
# Inverse Wavelet Transform
for (Xwt, fw) in ((:fwt,true),(:iwt,false))
@eval begin
    function $Xwt{T<:FloatingPoint}(x::AbstractMatrix{T}, L::Integer, filter::POfilter)
        n = size(x,1)
        n != size(x,2) && error("2D dwt: not a square matrix")   
        y = Array(T,n,n)
        dwt!(y,x,L,filter,$fw)
        return y
    end
    # assume full transform
    function $Xwt{T<:FloatingPoint}(x::AbstractMatrix{T}, filter::POfilter)
        n = size(x,1)
        n != size(x,2) && error("2D dwt: not a square matrix")   
        y = Array(T,n,n)
        dwt!(y,x,nscales(n),filter,$fw)
        return y
    end
    $Xwt{T<:Integer}(x::AbstractMatrix{T}, L::Integer, filter::POfilter) = $Xwt(float(x),L,filter)
    $Xwt{T<:Integer}(x::AbstractMatrix{T}, filter::POfilter) = $Xwt(float(x),nscales(size(x,1)),filter)
end
end



# ==== low level methods =====

# 1D
# writes to y
function dwt!{T<:FloatingPoint}(y::AbstractVector{T}, x::AbstractVector{T}, L::Integer, filter::POfilter, fw::Bool)
    si = Array(T, filter.n-1)       # tmp filter vector
    if fw
        dcfilter = mirror(convert(Vector{T},filter.qmf))  #mqmf
        scfilter = reverse(convert(Vector{T},filter.qmf))  #rqmf
    else
        scfilter = convert(Vector{T},filter.qmf)  #qmf
        dcfilter = reverse(mirror(convert(Vector{T},filter.qmf)))  #mrqmf
    end
    
    dwt!(y, x, L, filter, fw, dcfilter, scfilter, si)
    return nothing
end
function dwt!{T<:FloatingPoint}(y::AbstractVector{T}, x::AbstractVector{T}, L::Integer, filter::POfilter, fw::Bool, dcfilter::Vector{T}, scfilter::Vector{T}, si::Vector{T}, snew::Union(Vector{T},Nothing) = nothing)
    n = length(x)
    J = nscales(n)
    n != 2^J && error("length not a power of 2")
    n != length(y) && error("length of output does not match input")
    !(0 <= L <= J) && error("L out of bounds, use 0 <= L <= J")
    
    L == 0 && return copy!(y,x)     # do nothing
    L > 1 && snew == nothing && (snew = Array(T,n>>1)) # tmp storage of scaling coefs (not used for L<=1)
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
    return nothing
end

# 2D
# writes to y
function dwt!{T<:FloatingPoint}(y::AbstractMatrix{T}, x::AbstractMatrix{T}, L::Integer, filter::POfilter, fw::Bool)
    n = size(x,1)
    si = Array(T, filter.n-1)       # tmp filter vector
    tmpvec = Array(T,n)             # tmp storage vector
    if fw
        dcfilter = mirror(convert(Vector{T},filter.qmf))  #mqmf
        scfilter = reverse(convert(Vector{T},filter.qmf))  #rqmf
    else
        scfilter = convert(Vector{T},filter.qmf)  #qmf
        dcfilter = reverse(mirror(convert(Vector{T},filter.qmf)))  #mrqmf
    end
    
    dwt!(y, x, L, filter, fw, dcfilter, scfilter, si, tmpvec)
    return nothing
end
function dwt!{T<:FloatingPoint}(y::AbstractMatrix{T}, x::AbstractMatrix{T}, L::Integer, filter::POfilter, fw::Bool, dcfilter::Vector{T}, scfilter::Vector{T}, si::Vector{T}, tmpvec::Vector{T})

    n = size(x,1)
    J = nscales(n)
    n != size(x,2) && error("2D dwt: not a square matrix")
    size(x) != size(y) && error("input and output size mismatch")
    n != 2^J && error("length not a power of 2")
    !(0 <= L <= J) && error("L out of bounds, use 0 <= L <= J")
    
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
    return nothing
end



# periodic inplace filter and downsampling (by 2)
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
    (length(x) < ix+nx-1) && error("out of bounds")
    (length(out) < iout+nout-1) && error("out of bounds")
    ((nx-1)>>1 + iout > length(out)) && error("out of bounds")
    !(0 <= sw <= 1) && error("use sw = 0 or 1")
    silen != bs-1 && error("si: incorrect length")

    fill!(si,0.0)
    istart = bs + sw
    dsshift = (bs%2 + sw)%2  # is bs odd, and shift downsampling
    istart=istart
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

# periodic inplace filter and upsampling (by 2)
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
    (length(x) < ix+nx-1) && error("out of bounds")             # check array size
    (length(out) < iout+nout-1) && error("out of bounds")       # check array size
    ((nx<<1-1)>>1 + iout > length(out)) && error("out of bounds")  # max for out index
    !(0 <= sw <= 1) && error("use sw = 0 or 1")
    silen != bs-1 && error("si: incorrect length")
  
    fill!(si,0.0)
    dsshift = (sw)%2  # shift upsampling
    istart = bs - shift%2

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



end

