using ..Util

# WAVELET INDEXING AND SIZES

# Non-dyadic
# detail coef at level l location i (i starting at 1) -> vector index
detailindex(arraysize::Integer, l::Integer, i::Integer) = round(Int, arraysize / 2^l + i)
detailindex(x::AbstractArray, l::Integer, i::Integer) = detailindex(size(x, 1), l, i)
# the range of detail coefs at level l
detailrange(arraysize::Integer, l::Integer) = round(Int, (arraysize / 2^l + 1)):round(Int, arraysize / 2^(l - 1))
detailrange(x::AbstractArray, l::Integer) = detailrange(size(x, 1), l)
# number of detail coefs at level l
detailn(arraysize::Integer, l::Integer) = round(Int, arraysize / 2^l)
detailn(x::AbstractArray, l::Integer) = detailn(size(x, 1), l)
# max levels to transform
maxtransformlevels(x::AbstractArray) = minimum(maxtransformlevels.(size(x)))
function maxtransformlevels(arraysize::Integer)
    arraysize > 1 || return 0
    tl = 0
    while sufficientpoweroftwo(arraysize, tl)
        tl += 1
    end
    return tl - 1
end

maxmodwttransformlevels(x::AbstractArray) = floor(Int, log2(length(x)))
maxmodwttransformlevels(arraysize::Integer) = floor(Int, log2(arraysize))
