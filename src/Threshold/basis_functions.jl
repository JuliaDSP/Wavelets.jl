using ..Threshold
using LinearAlgebra: norm


# BASIS FUNCTIONS


"""
### Matching Pursuit
see: _Mallat (2009) p.642 "A wavelet tour of signal processing"_

Find sparse vector y such that ||x - f(y)|| < tol approximately

- `f` is the operation of a M by N (M<N) dictionary/matrix
- `ft` is a function defining the transpose of `f`
"""
function matchingpursuit(
    x::AbstractVector,
    f::Function,
    ft::Function,
    tol::Real,
    nmax::Int=-1,
    oop::Bool=false,
    N::Int=0
)

    @assert nmax >= -1
    @assert tol > 0
    r = x
    n = 1

    if !oop
        y = zeros(eltype(x), length(ft(x)))
    else # out of place functions f and ft
        y = zeros(eltype(x), N)
        tmp = similar(x, N)
        ftr = similar(x, N)
        aphi = similar(x, length(x))
    end
    spat = zeros(eltype(x), length(y))  # sparse for atom computation
    nmax == -1 && (nmax = length(y))

    while norm(r) > tol && n <= nmax
        # find largest inner product
        !oop && (ftr = ft(r))
        oop && ft(ftr, r, tmp)
        i = findmaxabs(ftr)

        # project on i-th atom
        spat[i] = ftr[i]
        !oop && (aphi = f(spat))
        oop && f(aphi, spat, tmp)
        spat[i] = 0

        # update residual, r = r - aphi
        r = r .- aphi

        y[i] += ftr[i]
        n += 1
    end
    return y
end


function findmaxabs(x::AbstractVector)
    m = abs(x[1])
    k = 1
    @inbounds for i in eachindex(x)
        if abs(x[i]) > m
            k = i
            m = abs(x[i])
        end
    end
    return k
end
