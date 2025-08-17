using ..Transforms


# ENTROPY MEASURES

abstract type Entropy end
struct ShannonEntropy <: Entropy end  #Coifman-Wickerhauser
struct LogEnergyEntropy <: Entropy end

# Entropy measures: Additive with coefentropy(0) = 0
# all coefs assumed to be on [-1,1] after normalization with nrm
# given x and y, where x has "more concentrated energy" than y
# then coefentropy(x, et, norm) <= coefentropy(y, et, norm) should be satisfied.

function coefentropy(x::T, et::ShannonEntropy, nrm::T) where {T<:AbstractFloat}
    s = (x / nrm)^2
    if s == 0.0
        return -zero(T)
    else
        return -s * log(s)
    end
end
function coefentropy(x::T, et::LogEnergyEntropy, nrm::T) where {T<:AbstractFloat}
    s = (x / nrm)^2
    if s == 0.0
        return -zero(T)
    else
        return -log(s)
    end
end
function coefentropy(x::AbstractArray{T}, et::Entropy, nrm::T=norm(x)) where {T<:AbstractFloat}
    @assert nrm >= 0
    sum = zero(T)
    nrm == sum && return sum
    for i in eachindex(x)
        @inbounds sum += coefentropy(x[i], et, nrm)
    end
    return sum
end


# find the best tree that is a subset of the input tree (use :full to find the best tree)
# for wpt
function bestbasistree(
    y::AbstractVector{T}, wt::DiscreteWavelet,
    L::Integer=maxtransformlevels(y),
    et::Entropy=ShannonEntropy()
) where {T<:AbstractFloat}
    bestbasistree(y, wt, maketree(length(y), L, :full), et)
end
function bestbasistree(
    y::AbstractVector{T}, wt::DiscreteWavelet,
    tree::BitVector,
    et::Entropy=ShannonEntropy()
) where {T<:AbstractFloat}

    isvalidtree(y, tree) || throw(ArgumentError("invalid tree"))

    x = copy(y)
    n = length(y)
    tmp = Vector{T}(undef, n)
    ntree = length(tree)
    entr_bf = Vector{T}(undef, ntree)
    nrm = norm(y)

    Lmax = maxtransformlevels(n)
    L = Lmax
    k = 1
    while L > 0
        ix = 1
        Lfw = Lmax - L
        nj = detailn(n, Lfw)

        @assert nj <= n
        dtmp = Transforms.unsafe_vectorslice(tmp, 1, nj)
        while ix <= n
            @assert nj + ix - 1 <= n
            dx = Transforms.unsafe_vectorslice(x, ix, nj)

            entr_bf[k] = coefentropy(dx, et, nrm)

            dwt!(dtmp, dx, wt, 1)
            copyto!(dx, dtmp)

            ix += nj
            k += 1
        end
        L -= 1
    end

    # entropy of fully transformed signal (end nodes)
    n_af = 2^(Lmax - 1)
    entr_af = Vector{T}(undef, n_af)
    n_coef_af = div(n, n_af)
    for i in 1:n_af
        range = (i-1)*n_coef_af+1:i*n_coef_af
        entr_af[i] = coefentropy(x[range], et, nrm)
    end

    # make the best tree
    besttree = copy(tree)
    for i in 1:ntree
        if (i > 1 && !besttree[i>>1]) || !tree[i]  # parent is 0 or input tree-node is 0
            besttree[i] = false
        else
            besttree[i] = (entr_bf[i] > bestsubtree_entropy(entr_bf, entr_af, i))
        end
    end

    @assert isvalidtree(y, besttree)
    return besttree::BitVector
end

# the entropy of best subtree
function bestsubtree_entropy(entr_bf::Array, entr_af::Array, i::Int)
    n = length(entr_bf)
    n_af = length(entr_af)
    @assert isdyadic(n + 1)
    @assert isdyadic(n_af)
    @assert n + 1 == n_af << 1

    if n < (i << 1)  # bottom of tree
        sum = entr_af[i-n_af+1]
    else
        sum = bestsubtree_entropy(entr_bf, entr_af, i << 1)
        sum += bestsubtree_entropy(entr_bf, entr_af, i << 1 + 1)
    end

    lowestentropy = min(entr_bf[i], sum)
    return lowestentropy
end
