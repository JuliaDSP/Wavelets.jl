
# Dyadic
# detail coef at level j location i (i starting at 1) -> vector index
dyadicdetailindex(j::Integer, i::Integer) = 2^j + i
# the range of detail coefs at level j
dyadicdetailrange(j::Integer) = (2^j+1):(2^(j+1))
# the range of scaling coefs at level j
dyadicscalingrange(j::Integer) = 1:2^j
# number of detail coefs at level j
dyadicdetailn(j::Integer) = 2^j
# number of scales of dyadic length signal (n=2^J)
ndyadicscales(n::Integer) = round(Int, log2(n))
ndyadicscales(x::AbstractArray) = ndyadicscales(size(x, 1))
# the largest detail level
maxdyadiclevel(n::Integer) = ndyadicscales(n) - 1
maxdyadiclevel(x::AbstractArray) = maxdyadiclevel(size(x, 1))
# convert number of transformed levels L to minimum level j
tl2dyadiclevel(n::Integer, L::Integer) = ndyadicscales(n) - L
tl2dyadiclevel(x::AbstractVector, L::Integer) = tl2dyadiclevel(length(x), L)
# convert maximum level j to number of transformed levels L
dyadiclevel2tl(arg...) = tl2dyadiclevel(arg...)
