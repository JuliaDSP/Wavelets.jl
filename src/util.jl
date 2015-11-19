module Util
export  dyadicdetailindex,
        dyadicdetailrange,
        dyadicscalingrange,
        dyadicdetailn,
        ndyadicscales,
        maxdyadiclevel,
        tl2dyadiclevel,
        dyadiclevel2tl,
        # non-dyadic
        detailindex,
        detailrange,
        detailn,
        maxtransformlevels,
        #
        mirror,
        upsample,
        downsample,
        iscube,
        isdyadic,
        sufficientpoweroftwo,
        wcount,
        stridedcopy!,
        isvalidtree,
        maketree,
        makewavelet,
        testfunction
using Compat

# WAVELET INDEXING AND SIZES

# Non-dyadic
# detail coef at level l location i (i starting at 1) -> vector index
detailindex(arraysize::Integer, l::Integer, i::Integer) = @compat round(Int, arraysize/2^l+i)
detailindex(x::AbstractArray, l::Integer, i::Integer) = detailindex(size(x,1), l ,i)
# the range of detail coefs at level l
detailrange(arraysize::Integer, l::Integer) = @compat round(Int, (arraysize/2^l+1)) : round(Int, arraysize/2^(l-1))
detailrange(x::AbstractArray, l::Integer) = detailrange(size(x,1), l)
# number of detail coefs at level l
detailn(arraysize::Integer, l::Integer) = @compat round(Int, arraysize/2^l)
detailn(x::AbstractArray, l::Integer) = detailn(size(x,1), l)
# max levels to transform
maxtransformlevels(x::AbstractArray) = maxtransformlevels(minimum(size(x)))
function maxtransformlevels(arraysize::Integer)
    arraysize > 1 || return 0
    tl = 0
    while sufficientpoweroftwo(arraysize, tl)
        tl += 1
    end
    return tl - 1
end

# Dyadic
# detail coef at level j location i (i starting at 1) -> vector index
dyadicdetailindex(j::Integer,i::Integer) = 2^j+i
# the range of detail coefs at level j
dyadicdetailrange(j::Integer) = (2^j+1):(2^(j+1))
# the range of scaling coefs at level j
dyadicscalingrange(j::Integer) = 1:2^j
# number of detail coefs at level j
dyadicdetailn(j::Integer) = 2^j
# number of scales of dyadic length signal (n=2^J)
ndyadicscales(n::Integer) = @compat round(Int, log2(n))
ndyadicscales(x::AbstractArray) = ndyadicscales(size(x,1))
# the largest detail level
maxdyadiclevel(n::Integer) = ndyadicscales(n)-1
maxdyadiclevel(x::AbstractArray) = maxdyadiclevel(size(x,1))
# convert number of transformed levels L to minimum level j
tl2dyadiclevel(n::Integer, L::Integer) = ndyadicscales(n)-L
tl2dyadiclevel(x::AbstractVector, L::Integer) = tl2dyadiclevel(length(x), L)
# convert maximum level j to number of transformed levels L
dyadiclevel2tl(arg...) = tl2dyadiclevel(arg...)



# UTILITY FUNCTIONS

# are all dimensions equal length?
function iscube(x::AbstractArray)
    for i = 1:ndims(x)
        size(x,1) != size(x,i) && return false
    end
    return true
end
# are all dimensions dyadic length?
function isdyadic(x::AbstractArray)
    for i = 1:ndims(x)
        isdyadic(size(x,i)) || return false
    end
    return true
end
isdyadic(n::Integer) = (n == 2^(ndyadicscales(n)))

# To perform a level L transform, the size of the signal in each dimension
# must have 2^L as a factor.
function sufficientpoweroftwo(x::AbstractArray, L::Integer)
    for i = 1:ndims(x)
        sufficientpoweroftwo(size(x,i), L) || return false
    end
    return true
end
sufficientpoweroftwo(n::Integer, L::Integer) = (n%(2^L) == 0)

# mirror of filter
mirror{T<:Number}(f::AbstractVector{T}) = f.*(-1).^(0:length(f)-1)
# upsample
function upsample(x::AbstractVector, sw::Int=0)
    @assert sw==0 || sw==1
    n = length(x)
    y = zeros(eltype(x), n<<1)
    sw -= 1

    for i = 1:n
        @inbounds y[i<<1 + sw] = x[i]
    end
    return y
end
# downsample
function downsample(x::AbstractVector, sw::Int=0)
    @assert sw==0 || sw==1
    n = length(x)
    @assert n%2 == 0
    y = zeros(eltype(x), n>>1)
    sw -= 1

    for i in eachindex(y)
        @inbounds y[i] = x[i<<1 + sw]
    end
    return y
end

# count coefficients above threshold t (>=), excluding coefficients in levels < level
# where level -1 is the x[1] coefficient
function wcount(x::AbstractVector, t::Real=0; level::Int=-1)
    @assert level >= -1
    c = 0
    si = 1
    level >= 0 && (si = 1 + 2^level)
    @inbounds for i = si:length(x)
        if abs(x[i]) >= t
            c += 1
        end
    end
    return c
end
# count coefficients above threshold t (>=)
function wcount(x::AbstractArray, t::Real=0)
    c = 0
    @inbounds for i in eachindex(x)
        if abs(x[i]) >= t
            c += 1
        end
    end
    return c
end

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
# out of place circular shift of vector by shift, such that b[i]=a[i-shift] (mod length(a))
function circshift!(b::AbstractVector, a::AbstractVector, shift::Integer)
    @assert length(a) == length(b)
    atype = typeof(a)
    s = length(a)
    shift = mod(shift,s)
    shift == 0 && return copy!(b,a)
    shift = ifelse(s>>1 < shift, shift-s, shift)  # shift a smaller distance if possible
    if shift < 0
        for i = 1:s+shift
            @inbounds b[i] = a[i-shift]
        end
        sh = s+shift
        for i = s+shift+1:s
            @inbounds b[i] = a[i-sh]
        end
    else
        for i = s:-1:1+shift
            @inbounds b[i] = a[i-shift]
        end
        sh = -s+shift
        for i = shift:-1:1
            @inbounds b[i] = a[i-sh]
        end
    end
    return b::atype
end

# put odd elements into first half, even into second half
function split!{T<:Number}(a::AbstractVector{T})
    n = length(a)
    nt = n>>2 + (n>>1)%2
    tmp = Array(T, nt)
    split!(a, n, tmp)
    return a
end

# split only the range 1:n
function split!{T<:Number}(a::AbstractVector{T}, n::Integer, tmp::Vector{T})
    @assert n <= length(a)
    @assert n%2 == 0
    n == 2 && return a
    nt = n>>2 + (n>>1)%2
    @assert nt <= length(tmp)

    for i=1:nt # store evens
        @inbounds tmp[i] = a[i<<1]
    end
    for i=1:n>>1  # odds to first part
        @inbounds a[i] = a[(i-1)<<1 + 1]
    end
    for i=0:nt - 1  # evens to end
        @inbounds a[n-i] = a[n - i<<1]
    end
    copy!(a,n>>1+1,tmp,1,nt)
    return a
end

# out of place split from a to b, only the range 1:n
function split!{T<:Number}(b::AbstractVector{T}, a::AbstractVector{T}, n::Integer)
    @assert n <= length(a) && n <= length(b)
    @assert n%2 == 0
    if n == 2
        b[1] = a[1]
        b[2] = a[2]
        return b
    end

    h = n>>1
    for i=1:h  # odds to b
        @inbounds b[i] = a[(i-1)<<1 + 1]
    end
    for i=h+1:n  # evens to b
        @inbounds b[i] = a[(i - h)<<1]
    end

    return b
end
# out of place split from a to b, only the range a[ia:inca:ia+(n-1)*inca] to b[1:n]
function split!{T<:Number}(b::AbstractVector{T}, a::AbstractArray{T}, ia::Integer, inca::Integer, n::Integer)
    @assert ia+(n-1)*inca <= length(a) && n <= length(b)
    @assert n%2 == 0
    if n == 2
        b[1] = a[ia]
        b[2] = a[ia + inca]
        return b
    end

    h = n>>1
    inca2 = inca<<1
    for i=1:h  # odds to b
        @inbounds b[i] = a[ia + (i-1)*inca2]
    end
    iainca = ia + inca
    hp1 = h + 1
    for i=h+1:n  # evens to b
        @inbounds b[i] = a[iainca + (i - hp1)*inca2]
    end

    return b
end

# inverse the operation of split!
function merge!{T<:Number}(a::AbstractVector{T})
    n = length(a)
    nt = n>>2 + (n>>1)%2
    tmp = Array(T,nt)
    merge!(a, n, tmp)
    return a
end

# merge only the range 1:n
function merge!{T<:Number}(a::AbstractVector{T}, n::Integer, tmp::Vector{T})
    @assert n <= length(a)
    @assert n%2 == 0
    n == 2 && return a
    nt = n>>2 + (n>>1)%2
    @assert nt <= length(tmp)

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
    @assert n <= length(a) && n <= length(b)
    @assert n%2 == 0
    if n == 2
        b[1] = a[1]
        b[2] = a[2]
        return b
    end

    h = n>>1
    for i=1:h  # odds to b
        @inbounds b[(i-1)<<1 + 1] = a[i]
    end
    for i=h+1:n  # evens to b
        @inbounds b[(i - h)<<1] = a[i]
    end

    return b
end
# out of place merge from a to b, only the range a[1:n] to b[ib:incb:ib+(n-1)*incb]
function merge!{T<:Number}(b::AbstractArray{T}, ib::Integer, incb::Integer, a::AbstractVector{T}, n::Integer)
    @assert n <= length(a) && ib+(n-1)*incb <= length(b)
    @assert n%2 == 0
    if n == 2
        b[ib] = a[1]
        b[ib + incb] = a[2]
        return b
    end

    h = n>>1
    incb2 = incb<<1
    for i=1:h  # odds to b
        @inbounds b[ib + (i-1)*incb2] = a[i]
    end
    ibincb = ib + incb
    hp1 = h + 1
    for i=h+1:n  # evens to b
        @inbounds b[ibincb + (i - hp1)*incb2] = a[i]
    end

    return b
end


function stridedcopy!{T<:Number}(b::AbstractVector{T}, a::AbstractArray{T}, ia::Integer, inca::Integer, n::Integer)
    @assert ia+(n-1)*inca <= length(a) && n <= length(b)

    @inbounds for i = 1:n
        b[i] = a[ia + (i-1)*inca]
    end
    return b
end
function stridedcopy!{T<:Number}(b::AbstractArray{T}, ib::Integer, incb::Integer, a::AbstractVector{T}, n::Integer)
    @assert ib+(n-1)*incb <= length(b) && n <= length(a)

    @inbounds for i = 1:n
        b[ib + (i-1)*incb] = a[i]
    end
    return b
end

# wavelet packet transforms WPT
# valid if 0 nodes have 0 children and length+1 is dyadic
# and has a nodes corresponding every transform level
function isvalidtree(x::AbstractVector, b::BitVector)
    ns = maxtransformlevels(x)
    nb = length(b)
    nb == 2^(ns)-1 || return false
    @assert (2^(ns-1)-1)<<1+1 <= nb

    for i in 1:2^(ns-1)-1
        @inbounds if !b[i] && (b[i<<1] || b[i<<1+1])
            return false
        end
    end
    return true
end
# return a tree (BitVector)
# s=:full, all nodes for first L levels equal 1, others 0
# s=:dwt, nodes corresponding to a dwt for first L levels equal 1, others 0
maketree(x::Vector, s::Symbol=:full) = maketree(length(x), maxtransformlevels(x), s)
function maketree(n::Int, L::Int, s::Symbol=:full)
    ns = maxtransformlevels(n)
    nb = 2^(ns)-1
    @assert 0 <= L <= ns
    @assert (2^(ns-1)-1)<<1+1 <= nb

    b = BitArray(nb)
    fill!(b, false)

    t = true
    if s == :full
        for i in 1:2^(L)-1
            @inbounds b[i] = t
        end
    elseif s == :dwt
        for i in 1:L
            @inbounds b[2^(i-1)] = t
        end
    else
        throw(ArgumentError("uknown symbol"))
    end
    return b
end

# return scaling and wavelet functions and location vector, made from filter h
# iterated with a cascade algorithm with N steps
makewavelet(h, arg...) = makewavelet(h.qmf, arg...)
function makewavelet(h::AbstractVector, N::Integer=8)
    @assert N>=0
    sc = norm(h)
    h = h*sqrt(2)/sc
    phi = copy(h)
    psi = mirror(reverse(h))

    for i=1:N
        phi = conv(upsample(phi), h)
        psi = conv(upsample(psi), h)
    end
    phi = phi[1:end-2^(N)+1]
    psi = psi[1:end-2^(N)+1]
    return scale!(phi,sc/sqrt(2)), scale!(psi,sc/sqrt(2)), linspace(0,length(h)-1,length(psi))
end

# return a vector of test function values on [0,1), see
#Donoho, D.L.; I.M. Johnstone (1994), "Ideal spatial adaptation by wavelet shrinkage," Biometrika, vol. 81, pp. 425–455.
function testfunction(n::Int, ft::AbstractString)
    @assert n >= 1
    f = Array(Float64,n)
    range = 0:1/n:1-eps()
    i = 1
    if ft=="Blocks"
        tj = [0.1,0.13,0.15,0.23,0.25,0.4,0.44,0.65,0.76,0.78,0.81]
        hj = [4,-5,3,-4,5,-4.2,2.1,4.3,-3.1,2.1,-4.2]
        for t in range
            f[i] = 0
            for j = 1:11
                f[i] += hj[j]*(1+sign(t-tj[j]))/2
            end
            i += 1
        end
    elseif ft=="Bumps"
        tj = [0.1,0.13,0.15,0.23,0.25,0.4,0.44,0.65,0.76,0.78,0.81]
        hj = [4,5,3,4,5,4.2,2.1,4.3,3.1,5.1,4.2]
        wj = [0.005,0.005,0.006,0.01,0.01,0.03,0.01,0.01,0.005,0.008,0.005]
        for t in range
            f[i] = 0
            for j = 1:11
                f[i] += hj[j]/(1+abs((t-tj[j])/wj[j]))^4
            end
            i += 1
        end
    elseif ft=="HeaviSine"
        for t in range
            f[i] = 4*sin(4*pi*t) - sign(t-0.3) - sign(0.72-t)
            i += 1
        end
    elseif ft=="Doppler"
        for t in range
            f[i] = sqrt(t*(1-t))*sin(2*pi*1.05/(t+0.05))
            i += 1
        end
    else
        throw(ArgumentError("unknown test function"))
    end
    return f
end


function copygeneral2!(y::AbstractArray, x::AbstractMatrix, t::Range, d::Integer)
    k = 1
    for i in t
        @inbounds y[k] = x[i,d]
        k +=1
    end
    return y
end
function copygeneral2!(y::AbstractArray, x::AbstractMatrix, d::Integer, t::Range)
    k = 1
    for i in t
        @inbounds y[k] = x[d,i]
        k +=1
    end
    return y
end
function copygeneral1!(y::AbstractMatrix, t::Range, d::Integer, x::AbstractArray)
    k = 1
    for i in t
        @inbounds y[i,d] = x[k]
        k +=1
    end
    return y
end
function copygeneral1!(y::AbstractMatrix, d::Integer, t::Range, x::AbstractArray)
    k = 1
    for i in t
        @inbounds y[d,i] = x[k]
        k +=1
    end
    return y
end


end # module
