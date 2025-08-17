
# BASIS FUNCTIONS

# Matching Pursuit
# see: Mallat (2009) p.642 "A wavelet tour of signal processing"
# find sparse vector y such that ||x - f(y)|| < tol approximately
# f is the operation of a M by N (M<N) dictionary/matrix
# ft is a function defining the transpose of f
function matchingpursuit(x::AbstractVector, f::Function, ft::Function, tol::Real, nmax::Int=-1, oop::Bool=false, N::Int=0)
    @assert nmax >= -1
    @assert tol > 0
    r = x
    n = 1

    if oop  # out of place functions f and ft
        y = zeros(eltype(x), N)
        tmp = similar(x, N)
        ftr = similar(x, N)
        aphi = similar(x, length(x))
    else
        y = zeros(eltype(x), length(ft(x)))
    end
    spat = zeros(eltype(x), length(y))  # sparse for atom computation
    nmax == -1 && (nmax = length(y))

    while norm(r) > tol && n <= nmax
        # find largest inner product
        if oop
            ft(ftr, r, tmp)  # compute f^T(r) in place
        else
            ftr = ft(r)  # compute f^T(r)
        end
        ftri, i = findmax(abs, ftr)

        # project on i-th atom
        spat[i] = ftri
        if oop
            f(aphi, spat, tmp)  # compute aphi = f(spat) in place
        else
            aphi = f(spat)  # compute aphi = f(spat)
        end
        spat[i] = 0

        # update residual, r = r - aphi
        broadcast!(-, r, r, aphi)

        y[i] += ftri
        n += 1
    end
    return y
end
