module Util
export detailindex, detailrange, detailn, nscales, mirror

# detail coef at level j location i (i starting at 1) -> vector index
detailindex(j::Integer,i::Integer) = 2^j+i
# the range of detail coefs for level j
detailrange(j::Integer) = (2^j+1):(2^(j+1))
# number of detail coefs at level j
detailn(j::Integer) = 2^j
# number of scales of dyadic length signal (n=2^J)
nscales(n::Integer) = int(log2(n))
# mirror of filter
mirror{T<:Number}(f::Vector{T}) = f.*(-1).^(0:length(f)-1)


end

