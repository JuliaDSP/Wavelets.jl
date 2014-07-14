module Util
export detailindex, detailrange, detailn, nscales, maxlevel, mirror, split!, merge!, testfunction

# detail coef at level j location i (i starting at 1) -> vector index
detailindex(j::Integer,i::Integer) = 2^j+i
# the range of detail coefs for level j
detailrange(j::Integer) = (2^j+1):(2^(j+1))
# number of detail coefs at level j
detailn(j::Integer) = 2^j
# number of scales of dyadic length signal (n=2^J)
nscales(n::Integer) = int(log2(n))
# the largest detial scale
maxlevel(n::Integer) = nscales(n)-1
# mirror of filter
mirror{T<:Number}(f::Vector{T}) = f.*(-1).^(0:length(f)-1)

# inplace circular shift of vector by shift, such that anew[i]=aold[i-shift] (mod length(a))
function circshift!(a::AbstractVector, shift::Integer)
    atype = typeof(a)
    s = length(a)
    shift = mod(shift,s)
    shift == 0 && return a
    shift = ifelse(s>>1 < shift, shift-s, shift)  # shift a smaller distance if possible
    if shift < 0
        tmp = a[1:-shift]
        for i = 1:s+shift
            @inbounds a[i] = a[i-shift]
        end
        a[s+1+shift:s] = tmp
    else
        tmp = a[s+1-shift:s]
        for i = s:-1:1+shift
            @inbounds a[i] = a[i-shift]
        end
        a[1:shift] = tmp
    end
    return a::atype
end

# put odd elements into first half, even into second half
function split!{T<:Number}(a::AbstractVector{T})
    n = length(a)
    nt = n>>2 + (n>>1)%2
    tmp = Array(T,nt)
    split!(a,tmp)
    return a
end
# split only the range 1:n
function split!{T<:Number}(a::AbstractVector{T}, n::Integer, tmp::Vector{T})
    n > length(a) && error("n too big")
    n == 2 && return nothing
    n%2 == 1 && error("must be even length")
    nt = n>>2 + (n>>1)%2
    nt > length(tmp) && error("tmp vector to small")

    for i=1:nt # store evens
        @inbounds tmp[i] = a[i<<1]
    end
    for i=1:n>>1  # odds to first part
        @inbounds a[i] = a[(i-1)<<1 + 1]
    end
    for i=0:nt - 1  # evens to end
        @inbounds a[n-i] = a[n - 2*i]
    end
    copy!(a,n>>1+1,tmp,1,nt)
    return a
end
# out of place split from a to b, only the range 1:n
function split!{T<:Number}(b::AbstractVector{T}, a::AbstractVector{T}, n::Integer)
    (n > length(a) || n > length(b)) && error("n too big")
    if n == 2
    	b[1] = a[1]
    	b[2] = a[2]
    	return b
    end
    n%2 == 1 && error("must be even length")
    h=n>>1
    for i=1:h  # odds to b
        @inbounds b[i] = a[(i-1)<<1 + 1]
    end
    for i=h+1:n  # evens to b
        @inbounds b[i] = a[2*(i - h)]
    end

    return b
end

# inverse the operation of split!
function merge!{T<:Number}(a::AbstractVector{T})
    n = length(a)
    nt = n>>2 + (n>>1)%2
    tmp = Array(T,nt)
    merge!(a,tmp)
    return a
end
# merge only the range 1:n
function merge!{T<:Number}(a::AbstractVector{T}, n::Integer, tmp::Vector{T})
    n > length(a) && error("n too big")
    n == 2 && return a
    n%2 == 1 && error("must be even length")
    nt = n>>2 + (n>>1)%2
    nt > length(tmp) && error("tmp vector to small")

    copy!(tmp,1,a,n>>1+1,nt)
    for i=nt-1:-1:0  # evens from end
        @inbounds a[n - 2*i] = a[n-i]
    end
    for i=n>>1:-1:1  # odds from first part
        @inbounds a[(i-1)<<1 + 1] = a[i]
    end
    for i=nt:-1:1 # retrieve evens
        @inbounds a[i<<1] = tmp[i]
    end
    return a
end
# out of place merge from a to b, only the range 1:n
function merge!{T<:Number}(b::AbstractVector{T}, a::AbstractVector{T}, n::Integer)
    (n > length(a) || n > length(b)) && error("n too big")
    if n == 2
    	b[1] = a[1]
    	b[2] = a[2]
    	return b
    end
    n%2 == 1 && error("must be even length")
    h=n>>1
    for i=1:h  # odds to b
        @inbounds b[(i-1)<<1 + 1] = a[i]
    end
    for i=h+1:n  # evens to b
        @inbounds b[2*(i - h)] = a[i]
    end

    return b
end


function testfunction(n::Int, t::String)
	if t=="HeaviSine"
		return [4*sin(4*pi*t)-sign(t-0.3)-sign(0.72-t) for t=0:1/n:1-eps()]
	else
		error("test function not found")
	end
end

end

