
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
            filtdown!(scfilter, si, y, 1,                detailn(j), s, 1, 0, 0)
        else
            # scaling coefficients
            filtup!(false,scfilter,  si, y, 1, detailn(j+1), s, 1, -filtlen+1, 0)
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
function dwt!{T<:FloatingPoint}(y::AbstractMatrix{T}, x::AbstractMatrix{T}, filter::OrthoFilter, L::Integer, fw::Bool)
    n = size(x,1)
    si = Array(T, length(filter)-1)       # tmp filter vector
    tmpvec = Array(T,n)             # tmp storage vector
    scfilter, dcfilter = makereverseqmf(filter.qmf, fw, T)
    
    dwt!(y, x, filter, L, fw, dcfilter, scfilter, si, tmpvec)
    return y
end
function dwt!{T<:FloatingPoint}(y::AbstractMatrix{T}, x::AbstractMatrix{T}, filter::OrthoFilter, L::Integer, fw::Bool, dcfilter::Vector{T}, scfilter::Vector{T}, si::Vector{T}, tmpvec::Vector{T})

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
                dwt!(ya, tmpsub, filter, 1, fw, dcfilter, scfilter, si)
            end
            # columns
            for i=1:nsub
                xi = 1+(i-1)*n
                xm = xi+nsub-1
                ya = sub(y, xi:xm)  # final dest and src
                copy!(tmpsub,1,ya,1,nsub)
                dwt!(ya, tmpsub, filter, 1, fw, dcfilter, scfilter, si)
            end       
        else
            # columns
            for i=1:nsub
                xi = 1+(i-1)*n
                xm = xi+nsub-1
                ya = sub(y, xi:xm)  # final dest and src
                copy!(tmpsub,1,ya,1,nsub)  # x has been copied to y
                dwt!(ya, tmpsub, filter, 1, fw, dcfilter, scfilter, si)
            end   
            # rows
            xs = n
            for i=1:nsub
                xi = i
                xm = n*(nsub-1)+i
                ya = sub(y, xi:xs:xm)  # final dest and src
                copy!(tmpsub,1,ya,1,nsub)
                dwt!(ya, tmpsub, filter, 1, fw, dcfilter, scfilter, si)
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

