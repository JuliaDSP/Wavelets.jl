# UTILITY FUNCTIONS

# are all dimensions equal length?
iscube(x::AbstractArray) = allequal(size(x))
# are all dimensions dyadic length?
isdyadic(x::AbstractArray) = all(isdyadic, size(x))
isdyadic(n::Integer) = ispow2(n)

# To perform a level L transform, the size of the signal in each dimension
# must have 2^L as a factor.
sufficientpoweroftwo(x::AbstractArray, L::Integer) =
    all(Base.Fix2(sufficientpoweroftwo, L), size(x))
sufficientpoweroftwo(n::Integer, L::Integer) = (n % (1<<L) == 0)

# mirror of filter
function mirror(h::AbstractVector{<:Number})
    g = Vector{eltype(h)}(undef, length(h))
    for i in eachindex(g, h)
        g[i] = iseven(i) ? -h[i] : h[i]
    end
    return g
end

# upsample
function upsample(x::AbstractVector, to_evens::Bool=false)
    n = length(x)
    y = zeros(eltype(x), 2n)
    sw = - 1 + to_evens

    for i in eachindex(x)
        y[2i+sw] = x[i]
    end
    return y
end
upsample(x::AbstractVector, i::Int) = upsample(x, Bool(i))
# downsample
function downsample(x::AbstractVector, from_evens::Bool=false)
    n = length(x)
    @assert iseven(n)
    y = zeros(eltype(x), n >> 1)
    sw = - 1 + from_evens

    for i in eachindex(y)
        y[i] = x[2i+sw]
    end
    return y
end
downsample(x::AbstractVector, i::Int) = downsample(x, Bool(i))

# count coefficients above threshold t (>=), excluding coefficients in levels < level
# where level -1 is the x[1] coefficient
function wcount(x::AbstractVector, t::Real=0; level::Int=-1)
    @assert level >= -1
    c = 0
    si = 1
    level >= 0 && (si = 1 + 2^level)
    for i = si:length(x)
        if abs(x[i]) >= t
            c += 1
        end
    end
    return c
end
# count coefficients above threshold t (>=)
wcount(x::AbstractArray, t::Real=0) = count(y -> abs(y) >= t, x)

# put odd elements into first half, even into second half
function split!(a::AbstractVector{T}) where T<:Number
    n = length(a)
    nt = n >> 2 + (n >> 1) % 2
    tmp = Vector{T}(undef, nt)
    split!(a, n, tmp)
    return a
end

# split only the range 1:n
function split!(a::AbstractVector{T}, n::Integer, tmp::Vector{T}) where T<:Number
    @assert n <= length(a)
    @assert iseven(n)
    n == 2 && return a
    nt = n >> 2 + (n >> 1) % 2
    @assert nt <= length(tmp)

    for i = 1:nt    # store evens
        @inbounds tmp[i] = a[2i]
    end
    for i = 1:n>>1  # odds to first part
        @inbounds a[i] = a[2i-1]
    end
    for i = 0:nt-1  # evens to end
        @inbounds a[n-i] = a[n-2i]
    end
    copyto!(a, n >> 1 + 1, tmp, 1, nt)
    return a
end

# out of place split from a to b, only the range 1:n
function split!(b::AbstractVector{T}, a::AbstractVector{T}, n::Integer) where T<:Number
    @assert n <= length(a) && n <= length(b)
    @assert iseven(n)
    if n == 2
        b[1] = a[1]
        b[2] = a[2]
        return b
    end

    h = n >> 1
    for i = 1:h     # odds to b
        b[i] = a[2i-1]
    end
    for i = h+1:n   # evens to b
        b[i] = a[2*(i-h)]
    end

    return b
end
# out of place split from a to b, only the range a[ia:inca:ia+(n-1)*inca] to b[1:n]
function split!(b::AbstractVector{T}, a::AbstractArray{T}, ia::Integer, inca::Integer, n::Integer) where T<:Number
    @assert ia + (n - 1) * inca <= length(a) && n <= length(b)
    @assert iseven(n)
    if n == 2
        b[1] = a[ia]
        b[2] = a[ia+inca]
        return b
    end

    h = n >> 1
    inca2 = 2inca
    for i = 1:h     # odds to b
        b[i] = a[ia+(i-1)*inca2]
    end
    iainca = ia + inca
    hp1 = h + 1
    for i = h+1:n   # evens to b
        b[i] = a[iainca+(i-hp1)*inca2]
    end

    return b
end

# inverse the operation of split!
function merge!(a::AbstractVector{T}) where T<:Number
    n = length(a)
    nt = n >> 2 + (n >> 1) % 2
    tmp = Vector{T}(undef, nt)
    merge!(a, n, tmp)
    return a
end

# merge only the range 1:n
function merge!(a::AbstractVector{T}, n::Integer, tmp::Vector{T}) where T<:Number
    @assert n <= length(a)
    @assert iseven(n)
    n == 2 && return a
    nt = n >> 2 + (n >> 1) % 2
    @assert nt <= length(tmp)

    copyto!(tmp, 1, a, n >> 1 + 1, nt)
    for i = nt-1:-1:0   # evens from end
        @inbounds a[n-2i] = a[n-i]
    end
    for i = n>>1:-1:1   # odds from first part
        @inbounds a[2i-1] = a[i]
    end
    for i = nt:-1:1     # retrieve evens
        @inbounds a[2i] = tmp[i]
    end
    return a
end

# out of place merge from a to b, only the range 1:n
function merge!(b::AbstractVector{T}, a::AbstractVector{T}, n::Integer) where T<:Number
    @assert n <= length(a) && n <= length(b)
    @assert iseven(n)
    if n == 2
        b[1] = a[1]
        b[2] = a[2]
        return b
    end

    h = n >> 1
    for i = 1:h     # odds to b
        b[2i-1] = a[i]
    end
    for i = h+1:n   # evens to b
        b[2*(i-h)] = a[i]
    end

    return b
end
# out of place merge from a to b, only the range a[1:n] to b[ib:incb:ib+(n-1)*incb]
function merge!(b::AbstractArray{T}, ib::Integer, incb::Integer, a::AbstractVector{T}, n::Integer) where T<:Number
    @assert n <= length(a) && ib + (n - 1) * incb <= length(b)
    @assert iseven(n)
    if n == 2
        b[ib] = a[1]
        b[ib+incb] = a[2]
        return b
    end

    h = n >> 1
    incb2 = incb << 1
    for i = 1:h     # odds to b
        b[ib+(i-1)*incb2] = a[i]
    end
    ibincb = ib + incb
    hp1 = h + 1
    for i = h+1:n   # evens to b
        b[ibincb+(i-hp1)*incb2] = a[i]
    end

    return b
end


function stridedcopy!(b::AbstractVector{<:Number}, a::AbstractArray{<:Number}, ia::Integer, inca::Integer, n::Integer)
    @assert ia + (n - 1) * inca <= length(a) && n <= length(b)

    for i = 1:n
        b[i] = a[ia+(i-1)*inca]
    end
    return b
end
function stridedcopy!(b::AbstractArray{<:Number}, ib::Integer, incb::Integer, a::AbstractVector{<:Number}, n::Integer)
    @assert ib + (n - 1) * incb <= length(b) && n <= length(a)

    for i = 1:n
        b[ib+(i-1)*incb] = a[i]
    end
    return b
end

# wavelet packet transforms WPT
# valid if 0 nodes have 0 children and length+1 is dyadic
# and has a nodes corresponding every transform level
function isvalidtree(x::AbstractVector, b::BitVector)
    ns = maxtransformlevels(x)
    nb = length(b)
    nb == 2^(ns) - 1 || return false
    @assert (2^(ns - 1) - 1) << 1 + 1 <= nb

    for i in 1:2^(ns-1)-1
        if !b[i] && (b[2i] || b[2i+1])
            return false
        end
    end
    return true
end
@doc """
    maketree(x::Vector, s::Symbol=:full)
    maketree(n::Int, L::Int, s::Symbol=:full)
return a tree (BitVector)
s=:full, all nodes for first L levels equal 1, others 0
s=:dwt, nodes corresponding to a dwt for first L levels equal 1, others 0
"""
maketree(x::AbstractVector, s::Symbol=:full) = maketree(length(x), maxtransformlevels(x), s)
function maketree(n::Int, L::Int, s::Symbol=:full)
    ns = maxtransformlevels(n)
    nb = 2^(ns) - 1
    @assert 0 <= L <= ns
    @assert (2^(ns - 1) - 1) << 1 + 1 <= nb

    b = BitArray(undef, nb)
    fill!(b, false)

    t = true
    if s == :full
        for i in 1:2^(L)-1
            b[i] = t
        end
    elseif s == :dwt
        for i in 1:L
            b[2^(i-1)] = t
        end
    else
        throw(ArgumentError("unknown symbol"))
    end
    return b
end

@doc """
    makewavelet(h::AbstractVector, N::Integer=8)
return scaling and wavelet functions and location vector, made from filter h
iterated with a cascade algorithm with N steps
"""
makewavelet(h, arg...) = makewavelet(h.qmf, arg...)
function makewavelet(h::AbstractVector, N::Integer=8)
    @assert N >= 0
    sc = norm(h)
    h = h * sqrt(2) / sc
    phi = copy(h)
    psi = mirror(reverse(h))

    for _ = 1:N
        phi = conv(upsample(phi), h)
        psi = conv(upsample(psi), h)
    end
    phi = phi[1:end-2^(N)+1]
    psi = psi[1:end-2^(N)+1]
    return rmul!(phi, sc / sqrt(2)), rmul!(psi, sc / sqrt(2)), range(0; stop=length(h) - 1, length=length(psi))
end

"""
    testfunction(n::Int, ft::AbstractString)
return a vector of test function values on [0,1), see Donoho, D.L.; I.M. Johnstone (1994), "Ideal spatial adaptation by wavelet shrinkage," Biometrika, vol. 81, pp. 425–455.

Options for ft are
* Blocks
* Bumps
* HeaviSine
* Doppler
"""
function testfunction(n::Int, ft::AbstractString)
    @assert n >= 1
    f = Vector{Float64}(undef, n)
    range = 0:1/n:1-eps()
    i = 1
    if ft == "Blocks"
        tj = [0.1, 0.13, 0.15, 0.23, 0.25, 0.4, 0.44, 0.65, 0.76, 0.78, 0.81]
        hj = [4, -5, 3, -4, 5, -4.2, 2.1, 4.3, -3.1, 2.1, -4.2]
        for t in range
            f[i] = 0
            for j = 1:11
                f[i] += hj[j] * (1 + sign(t - tj[j])) / 2
            end
            i += 1
        end
    elseif ft == "Bumps"
        tj = [0.1, 0.13, 0.15, 0.23, 0.25, 0.4, 0.44, 0.65, 0.76, 0.78, 0.81]
        hj = [4, 5, 3, 4, 5, 4.2, 2.1, 4.3, 3.1, 5.1, 4.2]
        wj = [0.005, 0.005, 0.006, 0.01, 0.01, 0.03, 0.01, 0.01, 0.005, 0.008, 0.005]
        for t in range
            f[i] = 0
            for j = 1:11
                f[i] += hj[j] / (1 + abs((t - tj[j]) / wj[j]))^4
            end
            i += 1
        end
    elseif ft == "HeaviSine"
        for t in range
            f[i] = 4 * sin(4 * pi * t) - sign(t - 0.3) - sign(0.72 - t)
            i += 1
        end
    elseif ft == "Doppler"
        for t in range
            f[i] = sqrt(t * (1 - t)) * sin(2 * pi * 1.05 / (t + 0.05))
            i += 1
        end
    else
        throw(ArgumentError("unknown test function"))
    end
    return f
end
