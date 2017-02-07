module Plot
using ..Util, ..WT, ..Transforms
export wplotdots, wplotim

# PLOTTING UTILITIES

# return levels and detail coefficient centers on the interval [0,r) above (>=) threshold t
# as tuple (d,l)
function wplotdots(x::AbstractVector, t::Real=0, r::Real=1)
    isdyadic(x) || throw(ArgumentError("array must be of dyadic size"))
	n = length(x)
	c = wcount(x, t, level=0)
    d = Array(Float64, c)
    l = Array(Int, c)
    range = 0:1/n:1-eps()
    range *= r

    J = ndyadicscales(n)
    k = 1
    @inbounds begin
        for j = 0:J-1
            rind = 2^(J-1-j):2^(J-j):n
            for i in 1:dyadicdetailn(j)
                if abs(x[dyadicdetailindex(j,i)]) >= t
                    d[k] = range[rind[i]]
                    l[k] = j
                    k += 1
                end
            end
        end
    end
    return (d,l)
end

# return a matrix of detail coefficient values where row j+1 is level j
function wplotim(x::AbstractVector)
    isdyadic(x) || throw(ArgumentError("array must be of dyadic size"))
    n = length(x)
    J = ndyadicscales(n)
    A = zeros(Float64, J, n)

    @inbounds begin
        for j = 0:J-1
	        dr = dyadicdetailrange(j)
	        m = 2^(J-j)
	        for i in eachindex(dr)
		        A[j+1,1+(i-1)*m:i*m] = x[dr[i]]
	        end
        end
    end
    return A
end

# return an array of scaled detail coefficients and unscaled scaling coefficients
# ready to be plotted as an image
function wplotim(x::AbstractArray, L::Integer, wt::Union{DiscreteWavelet,Void}=nothing;
                wabs::Bool=true, power::Real=0.7, pnorm::Real=1)
    isdyadic(x) || throw(ArgumentError("array must be of dyadic size"))
    dim = ndims(x)
    (dim == 2 || dim == 3) || throw(ArgumentError("dimension $(dim) not supported"))
    n = size(x,1)
    cn = size(x,3)  # color dimension
    n == size(x,2) || throw(ArgumentError("array must be square"))
    (cn == 1 || cn == 3) || throw(ArgumentError("third dimension $(cn)  not supported"))
    J = ndyadicscales(n)
    nsc = 2^(J-L)

    # do wavelet transform
    if wt != nothing
        if size(x,3)>1
            x = dwtc(x, wt, L)
        else
            x = dwt(x, wt, L)
        end
    end

    # scaling coefs
    scs = x[1:nsc,1:nsc,:]
    scale01!(scs)

    # detail coefs
    xts = x
    wabs && (broadcast!(abs,xts,xts))
    xts[1:nsc,1:nsc,:] = 0
    scale01!(xts)
    for j=1:n, i=1:n
        @inbounds xts[i,j,:] = vecnorm(xts[i,j,:],pnorm).^(power)
    end

    # merge and reshape final image
    scale01!(xts)
    xts[1:nsc,1:nsc,:] = scs
    return xts
end
# scale elements of z to the interval [0,1]
function scale01!(z)
    mi = minimum(z)
    ma = maximum(z)
    for i in eachindex(z)
        @inbounds z[i] = (z[i] - mi)/(ma-mi)
    end
    return z
end



end
