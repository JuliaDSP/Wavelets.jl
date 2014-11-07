
##################################################################################
#
#  ORTHOFILTER TRANSFORMS 
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
function dwt!{T<:FloatingPoint}(y::AbstractVector{T}, x::AbstractVector{T}, filter::OrthoFilter, L::Integer, fw::Bool)
    si = Array(T, length(filter)-1)       # tmp filter vector
    scfilter, dcfilter = makereverseqmf(filter.qmf, fw, T)
    
    dwt!(y, x, filter, L, fw, dcfilter, scfilter, si)
    return y
end
# pseudo "inplace" by copying
function dwt!{T<:FloatingPoint}(x::AbstractArray{T}, filter::OrthoFilter, L::Integer, fw::Bool)
    y = Array(T, size(x))
    dwt!(y, x, filter, L, fw)
    copy!(x,y)
    return x
end
function dwt!{T<:FloatingPoint}(y::AbstractVector{T}, x::AbstractVector{T}, filter::OrthoFilter, L::Integer, fw::Bool, dcfilter::Vector{T}, scfilter::Vector{T}, si::Vector{T}, snew::Vector{T} = Array(T, ifelse(L>1, length(x)>>1, 0)))
    n = length(x)
    J = nscales(n)
    @assert size(x) == size(y)
    @assert isdyadic(y)
    @assert 0 <= L <= J
    is(y,x) && error("input vector is output vector")
    
    L == 0 && return copy!(y,x)     # do nothing
    L > 1 && length(snew) != n>>1 && error("length of snew incorrect")
    s = x                           # s is currect scaling coefs location
    filtlen = length(filter)
    
    if fw
        jrange = (J-1):-1:(J-L)
    else
        jrange = (J-L):(J-1)
    end
    for j in jrange
        if fw
            # detail coefficients
            filtdown!(dcfilter, si, y, detailindex(j,1), detailn(j), s, 1,-filtlen+1,1)
            # scaling coefficients
            filtdown!(scfilter, si, y,                1, detailn(j), s, 1, 0, 0)
        else
            # scaling coefficients
            filtup!(false, scfilter, si, y, 1, detailn(j+1), s, 1, -filtlen+1, 0)
            # detail coefficients
            filtup!(true,  dcfilter, si, y, 1, detailn(j+1), x, detailindex(j,1), 0, 1)
        end
        # if not final iteration: copy to tmp location
        fw  && j != jrange[end] && copy!(snew,1,y,1,detailn(j))
        !fw && j != jrange[end] && copy!(snew,1,y,1,detailn(j+1))
        L > 1 && (s = snew)
    end
    return y
end
function unsafe_dwt1level!{T<:FloatingPoint}(y::AbstractVector{T}, x::AbstractVector{T}, filter::OrthoFilter, fw::Bool, dcfilter::Vector{T}, scfilter::Vector{T}, si::Vector{T})
    n = length(x)
    j = nscales(n) - 1
    filtlen = length(filter)

    if fw
        # detail coefficients
        filtdown!(dcfilter, si, y, detailindex(j,1), detailn(j), x, 1,-filtlen+1,1)
        # scaling coefficients
        filtdown!(scfilter, si, y,                1, detailn(j), x, 1, 0, 0)
    else
        # scaling coefficients
        filtup!(false, scfilter, si, y, 1, detailn(j+1), x, 1, -filtlen+1, 0)
        # detail coefficients
        filtup!(true,  dcfilter, si, y, 1, detailn(j+1), x, detailindex(j,1), 0, 1)
    end
    return y
end

# 2-d
# writes to y
function dwt!{T<:FloatingPoint}(y::Matrix{T}, x::AbstractMatrix{T}, filter::OrthoFilter, L::Integer, fw::Bool)
    n = size(x,1)
    si = Array(T, length(filter)-1)       # tmp filter vector
    tmpvec = Array(T,n<<1)             # tmp storage vector
    scfilter, dcfilter = makereverseqmf(filter.qmf, fw, T)
    
    dwt!(y, x, filter, L, fw, dcfilter, scfilter, si, tmpvec)
    return y
end
function dwt!{T<:FloatingPoint}(y::Matrix{T}, x::AbstractMatrix{T}, filter::OrthoFilter, L::Integer, fw::Bool, dcfilter::Vector{T}, scfilter::Vector{T}, si::Vector{T}, tmpvec::Vector{T})

    n = size(x,1)
    J = nscales(n)
    @assert size(x) == size(y)
    @assert iscube(y)
    @assert isdyadic(y)
    @assert 0 <= L <= J
    @assert length(tmpvec) >= n<<1
    is(y,x) && error("input matrix is output matrix")
    
    L == 0 && return copy!(y,x)        # do nothing
    xs = n
    #s = x

    if fw
        jrange = (J-1):-1:(J-L)
        nsub = n
    else
        jrange = (J-L):(J-1)
        nsub = int(2^(J-L+1))
        copy!(y,x)
    end

    for j in jrange
        tmpsub = unsafe_vectorslice(tmpvec, 1, nsub)
        tmpsub2 = unsafe_vectorslice(tmpvec, nsub+1, nsub)
        if fw
            # rows
            for i=1:nsub
                xi = i
                if j != jrange[1]
                    stridedcopy!(tmpsub, y, xi, xs, nsub)
                else  # use x in first iteration
                    stridedcopy!(tmpsub, x, xi, xs, nsub)
                end
                unsafe_dwt1level!(tmpsub2, tmpsub, filter, fw, dcfilter, scfilter, si)
                stridedcopy!(y, xi, xs, tmpsub2, nsub)
            end
            # columns
            for i=1:nsub
                xi = 1+(i-1)*n
                ya = unsafe_vectorslice(y, xi, nsub)
                copy!(tmpsub,1,ya,1,nsub)
                unsafe_dwt1level!(ya, tmpsub, filter, fw, dcfilter, scfilter, si)
            end       
        else
            # columns
            for i=1:nsub
                xi = 1+(i-1)*n
                ya = unsafe_vectorslice(y, xi, nsub)
                if j != jrange[1]
                    copy!(tmpsub,1,ya,1,nsub)
                    unsafe_dwt1level!(ya, tmpsub, filter, fw, dcfilter, scfilter, si)
                else  # use x in first iteration
                    #xsub = unsafe_vectorslice(x, xi, nsub)  # not safe for AbstractArray
                    copy!(tmpsub,1,x,xi,nsub)
                    unsafe_dwt1level!(ya, tmpsub, filter, fw, dcfilter, scfilter, si)
                end
            end   
            # rows
            for i=1:nsub
                xi = i
                stridedcopy!(tmpsub, y, xi, xs, nsub)
                unsafe_dwt1level!(tmpsub2, tmpsub, filter, fw, dcfilter, scfilter, si)
                stridedcopy!(y, xi, xs, tmpsub2, nsub)
            end

        end 

        fw  && (nsub = nsub>>1)
        !fw && (nsub = nsub<<1)
    end
    return y
end

macro filtermainloop(si, silen, b, val)
    quote
        @inbounds for j=2:$(esc(silen))
            $(esc(si))[j-1] = $(esc(si))[j] + $(esc(b))[j]*$(esc(val))
        end
        @inbounds si[$(esc(silen))] = $(esc(b))[$(esc(silen))+1]*$(esc(val))
    end
end
macro filtermainloopzero(si, silen)
    quote
        @inbounds for j=2:$(esc(silen))
            $(esc(si))[j-1] = $(esc(si))[j]
        end
        @inbounds si[$(esc(silen))] = 0.0
    end
end

# periodic filter and downsampling (by 2)
# b : filter
# si: tmp array of length 1 less than b
# out : result gets written to out[iout:nout-1], where nout==nx/2
# x : filter convolved with x[ix:nx-1] (shifted by shift)
# sw : shift downsampling (0 or 1)
# based on Base.filt
function filtdown!{T<:FloatingPoint}(b::Vector{T}, si::Vector{T},
                              out::AbstractVector{T},  iout::Integer, nout::Integer, 
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
        for i = 1:istart-1
            # periodic in the range [ix:ix+nx-1]
            xatind = x[mod(i-1+shift, nx) + ix]
            @filtermainloop(si, silen, b, xatind)
        end
        for i = istart:(nx-1+istart)
            # periodic in the range [ix:ix+nx-1]
            xatind = x[mod(i-1+shift, nx) + ix]
            # take every other value after istart
            if (i+dsshift)%2 == 0
                out[(i-istart)>>1 + iout] = si[1] + b[1]*xatind
            end
            @filtermainloop(si, silen, b, xatind)
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
                              out::AbstractVector{T},  iout::Integer, nout::Integer, 
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
        xatind = 0.0
        for i = 1:istart-1
            if (i+dsshift)%2==0
                @filtermainloopzero(si, silen)
            else
                # periodic in the range [ix:ix+nx-1]
                xindex = mod((i-1)>>1+shift>>1,nx) + ix    #(i-1)>>1 increm. every other
                xatind = x[xindex]
                @filtermainloop(si, silen, b, xatind)
            end
        end
        for i = istart:(nout-1+istart)
            if (i+dsshift)%2==0
                xatind = 0.0
            else
                xindex = mod((i-1)>>1+shift>>1,nx) + ix
                xatind = x[xindex]
            end
            if add2out
                out[(i-istart) + iout] += si[1] + b[1]*xatind
            else
                out[(i-istart) + iout] = si[1] + b[1]*xatind
            end
            if (i+dsshift)%2==0
                @filtermainloopzero(si, silen)
            else
                @filtermainloop(si, silen, b, xatind)
            end
        end
    end

    return nothing
end



