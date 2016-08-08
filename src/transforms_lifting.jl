
##################################################################################
#
#  LIFTING TRANSFORMS
#  Periodic boundaries
#
##################################################################################

# for split and merge
reqtmplength(x::AbstractArray) = (size(x,1)>>2) + (size(x,1)>>1)%2

# return scheme parameters adjusted for direction and type
function makescheme{T<:Number}(::Type{T}, scheme::GLS, fw::Bool)
    n = length(scheme.step)
    stepseq = Array(WT.LSStep{T}, n)
    for i = 1:n
        j = fw ? i : n+1-i
        stepseq[i] = WT.LSStep(scheme.step[j].steptype,
                            convert(Vector{T}, scheme.step[j].param.coef),
                            scheme.step[j].param.shift)
    end
    norm1, norm2 =  convert(T, fw ? scheme.norm1 : 1/scheme.norm1),
                    convert(T, fw ? scheme.norm2 : 1/scheme.norm2)
    return stepseq, norm1, norm2
end


# 1-D
# inplace transform of y, no vector allocation
function _dwt!{T<:Number}(y::AbstractVector{T}, scheme::GLS, L::Integer, fw::Bool, tmp::Vector{T}=Array(T,reqtmplength(y)))

    n = length(y)
    0 <= L ||
        throw(ArgumentError("L must be positive"))
    sufficientpoweroftwo(y, L) ||
        throw(ArgumentError("size must have a sufficient power of 2 factor"))
    length(tmp) >= n>>2 ||
        throw(ArgumentError("length of tmp incorrect"))

    if L == 0
        return y
    end

    if fw
        lrange = 1:L
        ns = n
    else
        lrange = L:-1:1
        ns = detailn(n, L-1)
    end
    half = ns>>1
    s = y
    stepseq, norm1, norm2 = makescheme(T, scheme, fw)

    for l in lrange
        if fw
            Util.split!(s, ns, tmp)
            for step in stepseq
                liftfw!(s, half, step.param, step.steptype)
            end
            normalize!(s, half, ns, norm1, norm2)
            ns = ns>>1
            half = half>>1
        else
            normalize!(s, half, ns, norm1, norm2)
            for step in stepseq
                liftbw!(s, half, step.param, step.steptype)
            end
            Util.merge!(s, ns, tmp)        # inverse split
            ns = ns<<1
            half = half<<1
        end
    end
    return y
end
# 1-D, unsafe simple 1 level transform
# inplace transform of y, no vector allocation
# tmp: size at least n>>2
# oopc: use oop computation, if false iy and incy are assumed to be 1
# oopv: the out of place location
function unsafe_dwt1level!{T<:Number}(y::AbstractArray{T}, iy::Integer, incy::Integer, oopc::Bool, oopv::Vector{T}, scheme::GLS, fw::Bool, stepseq::Vector, norm1::T, norm2::T, tmp::Vector{T})
    if !oopc
        oopv = y
    end
    ns = length(oopv)
    half = ns>>1

    if fw
        if oopc
            Util.split!(oopv, y, iy, incy, ns)
        else
            Util.split!(oopv, ns, tmp)
        end
        for step in stepseq
            liftfw!(oopv, half, step.param, step.steptype)
        end
        if oopc
            normalize!(y, iy, incy, oopv, half, ns, norm1, norm2)
        else
            normalize!(oopv, half, ns, norm1, norm2)
        end
    else
        if oopc
            normalize!(oopv, y, iy, incy, half, ns, norm1, norm2)
        else
            normalize!(oopv, half, ns, norm1, norm2)
        end
        for step in stepseq
            liftbw!(oopv, half, step.param, step.steptype)
        end
        if oopc
            Util.merge!(y, iy, incy, oopv, ns)
        else
            Util.merge!(oopv, ns, tmp)
        end
    end

    return y
end

# 2-D
# inplace transform of y, no vector allocation
# tmp: size at least n>>2
# tmpvec: size at least n
function _dwt!{T<:Number}(y::Matrix{T}, scheme::GLS, L::Integer, fw::Bool, tmp::Vector{T}=Array(T,reqtmplength(y)), tmpvec::Vector{T}=Array(T,size(y,1)))

    n = size(y,1)
    iscube(y) ||
        throw(ArgumentError("array must be square/cube"))
    0 <= L ||
        throw(ArgumentError("L must be positive"))
    sufficientpoweroftwo(y, L) ||
        throw(ArgumentError("size must have a sufficient power of 2 factor"))
    (length(tmp) >= n>>2) ||
        throw(ArgumentError("length of tmp incorrect"))
    (length(tmpvec) >= n) ||
        throw(ArgumentError("length of tmpvec incorrect"))

    if L == 0
        return y
    end
    row_stride = n
    height_stride = n*n

    if fw
        lrange = 1:L
        nsub = n
    else
        lrange = L:-1:1
        nsub = div(n,2^(L-1))
    end
    stepseq, norm1, norm2 = makescheme(T, scheme, fw)

    # transforms with stride are out of place in a dense array for speed
    for l in lrange
        tmpsub = unsafe_vectorslice(tmpvec, 1, nsub)
        if fw
            # rows
            for i in 1:nsub
                xi = row_idx(i, n)
                unsafe_dwt1level!(y, xi, row_stride, true, tmpsub, scheme, fw,
                                    stepseq, norm1, norm2, tmp)
            end
            # columns
            for i in 1:nsub
                xi = col_idx(i, n)
                ya = unsafe_vectorslice(y, xi, nsub)
                unsafe_dwt1level!(ya, 1, 1, false, tmpsub, scheme, fw,
                                    stepseq, norm1, norm2, tmp)
            end
        else
            # columns
            for i in 1:nsub
                xi = col_idx(i, n)
                ya = unsafe_vectorslice(y, xi, nsub)
                unsafe_dwt1level!(ya, 1, 1, false, tmpsub, scheme, fw,
                                    stepseq, norm1, norm2, tmp)
            end
            # rows
            for i in 1:nsub
                xi = row_idx(i, n)
                unsafe_dwt1level!(y, xi, row_stride, true, tmpsub, scheme, fw,
                                    stepseq, norm1, norm2, tmp)
            end
        end

        nsub = (fw ? nsub>>1 : nsub<<1)
    end

    return y
end

# 3-D
# inplace transform of y, no vector allocation
# tmp: size at least n>>2
# tmpvec: size at least n
function _dwt!{T<:Number}(y::Array{T,3}, scheme::GLS, L::Integer, fw::Bool, tmp::Vector{T}=Array(T,reqtmplength(y)), tmpvec::Vector{T}=Array(T,size(y,1)))

    n = size(y,1)
    iscube(y) ||
        throw(ArgumentError("array must be square/cube"))
    0 <= L ||
        throw(ArgumentError("L must be positive"))
    sufficientpoweroftwo(y, L) ||
        throw(ArgumentError("size must have a sufficient power of 2 factor"))
    (length(tmp) >= n>>2) ||
        throw(ArgumentError("length of tmp incorrect"))
    (length(tmpvec) >= n) ||
        throw(ArgumentError("length of tmpvec incorrect"))

    if L == 0
        return y
    end
    row_stride = n
    height_stride = n*n

    if fw
        lrange = 1:L
        nsub = n
    else
        lrange = L:-1:1
        nsub = div(n,2^(L-1))
    end
    stepseq, norm1, norm2 = makescheme(T, scheme, fw)

    # transforms with stride are out of place in a dense array for speed
    for l in lrange
        tmpsub = unsafe_vectorslice(tmpvec, 1, nsub)
        if fw
            # heights
            for i in 1:nsub, j in 1:nsub
                xi = hei_idx(i, j, n)
                unsafe_dwt1level!(y, xi, height_stride, true, tmpsub, scheme, fw,
                                    stepseq, norm1, norm2, tmp)
            end
            # rows
            for i in 1:nsub, j in 1:nsub
                xi = row_idx(i, j, n)
                unsafe_dwt1level!(y, xi, row_stride, true, tmpsub, scheme, fw,
                                    stepseq, norm1, norm2, tmp)
            end
            # columns
            for i in 1:nsub, j in 1:nsub
                xi = col_idx(i, j, n)
                ya = unsafe_vectorslice(y, xi, nsub)
                unsafe_dwt1level!(ya, 1, 1, false, tmpsub, scheme, fw,
                                    stepseq, norm1, norm2, tmp)
            end
        else
            # columns
            for i in 1:nsub, j in 1:nsub
                xi = col_idx(i, j, n)
                ya = unsafe_vectorslice(y, xi, nsub)
                unsafe_dwt1level!(ya, 1, 1, false, tmpsub, scheme, fw,
                                    stepseq, norm1, norm2, tmp)
            end
            # rows
            for i in 1:nsub, j in 1:nsub
                xi = row_idx(i, j, n)
                unsafe_dwt1level!(y, xi, row_stride, true, tmpsub, scheme, fw,
                                    stepseq, norm1, norm2, tmp)
            end
            # heights
            for i in 1:nsub, j in 1:nsub
                xi = hei_idx(i, j, n)
                unsafe_dwt1level!(y, xi, height_stride, true, tmpsub, scheme, fw,
                                    stepseq, norm1, norm2, tmp)
            end
        end

        nsub = (fw ? nsub>>1 : nsub<<1)
    end

    return y
end


# WPT
# 1-D
function _wpt!{T<:Number}(y::AbstractVector{T}, scheme::GLS, tree::BitVector, fw::Bool, tmp::Vector{T}=Array(T,reqtmplength(y)))

    n = length(y)
    isvalidtree(y, tree) ||
        throw(ArgumentError("invalid tree"))
    length(tmp) >= n>>2 ||
        throw(ArgumentError("length of tmp incorrect"))

    if !tree[1]
        return y
    end

    stepseq, norm1, norm2 = makescheme(T, scheme, fw)

    Lmax = maxtransformlevels(n)
    L = Lmax
    while L > 0
        ix = 1
        k = 1
        fw  && (Lfw = Lmax-L)
        !fw && (Lfw = L-1)
        nj = detailn(n, Lfw)
        treeind = 2^(Lfw)-1

        while ix <= n
            if tree[treeind+k]
                dy = unsafe_vectorslice(y, ix, nj)
                unsafe_dwt1level!(dy, 1, 1, false, tmp, scheme, fw, stepseq, norm1, norm2, tmp)
            end
            ix += nj
            k += 1
        end
        L -= 1
    end

    return y
end



function normalize!{T<:Number}(x::AbstractVector{T}, half::Int, ns::Int, n1::T, n2::T)
    for i = 1:half
        @inbounds x[i] *= n1
    end
    for i = half+1:ns
        @inbounds x[i] *= n2
    end
    return x
end
# out of place normalize from x to y
function normalize!{T<:Number}(y::AbstractVector{T}, x::AbstractVector{T}, half::Int, ns::Int, n1::T, n2::T)
    for i = 1:half
        @inbounds y[i] = n1*x[i]
    end
    for i = half+1:ns
        @inbounds y[i] = n2*x[i]
    end
    return y
end
function normalize!{T<:Number}(y::AbstractArray{T}, iy::Int, incy::Int, x::AbstractVector{T}, half::Int, ns::Int, n1::T, n2::T)
    for i = 1:half
        @inbounds y[iy + (i-1)*incy] = n1*x[i]
    end
    for i = half+1:ns
        @inbounds y[iy + (i-1)*incy] = n2*x[i]
    end
    return y
end
function normalize!{T<:Number}(y::AbstractVector{T}, x::AbstractArray{T}, ix::Int, incx::Int, half::Int, ns::Int, n1::T, n2::T)
    for i = 1:half
        @inbounds y[i] = n1*x[ix + (i-1)*incx]
    end
    for i = half+1:ns
        @inbounds y[i] = n2*x[ix + (i-1)*incx]
    end
    return y
end


# predict and update lifting steps inplace on x, forward and backward
# half: half of the length under consideration
# For predict: writes to range 1:half, reads from 1:2*half
# For update : writes to range half+1:2*half, reads from 1:2*half
for (fname, lift_inb, lift_bound) in (	(:liftfw!, :lift_inboundsfw!, :lift_perboundaryfw!),
										(:liftbw!, :lift_inboundsbw!, :lift_perboundarybw!) ),
    step_type in (WT.PredictStep, WT.UpdateStep)
@eval begin
function ($fname){T<:Number}(x::AbstractVector{T}, half::Int,
                                    param::WT.LSStepParam{T}, steptype::$step_type)
    lhsr, irange, rhsr, rhsis = getliftranges(half, length(param), param.shift, steptype)
    coefs = param.coef
    # left boundary
    ($lift_bound)(x, half, coefs, lhsr, rhsis, steptype)
    # main loop
    ($lift_inb)(x, coefs, irange, rhsis)
    # right boundary
    ($lift_bound)(x, half, coefs, rhsr, rhsis, steptype)
    return x
end
end # eval begin
end # for

function getliftranges(half::Int, nc::Int, shift::Int, steptype::WT.StepType)
    # define index shift rhsis
    pred = isa(steptype, typeof(WT.Predict))
    if pred
        rhsis = -shift+half
    else
        rhsis = -shift-half
    end
    # conditions for every element i in irange to be in bounds
    # 1 <= i <= half
    # 1 <= i+1-1-shift <= half
    # 1 <= i+nc-1-shift <= half
    irmin = max(shift+1, 1-nc+shift)
    irmax = min(half+1+shift-nc, half+shift)
    if irmin > half || irmax < 1
        irange = 1:0  # empty
    else
        irmin = max(irmin,1)
        irmax = min(irmax,half)
        irange = irmin:irmax
    end
    # periodic boundary
    if length(irange)==0
        lhsr = 1:half
        rhsr = 1:0
    else
        lhsr = 1:irmin-1
        rhsr = irmax+1:half
    end
    if !pred  # shift ranges for update
        irange += half
        lhsr += half
        rhsr += half
    end
    return (lhsr, irange, rhsr, rhsis)
end

# periodic boundary
for (fname, op) in ( (:lift_perboundaryfw!, :-), (:lift_perboundarybw!, :+) ),
	(step_type, puxind) in ((WT.PredictStep, :(mod1(i+k-1+rhsis-half,half)+half)),
                            (WT.UpdateStep,  :(mod1(i+k-1+rhsis,half))) )
@eval begin
function ($fname){T<:Number}(x::AbstractVector{T}, half::Int,
									c::Vector{T}, irange::Range, rhsis::Int, ::$step_type)
    nc = length(c)
    for i in irange
        for k = 1:nc
            @inbounds x[i] = ($op)(x[i], c[k]*x[$puxind] )
        end
    end
    return x
end
end # eval begin
end # for


# main lift loop
for (fname,op) in ( (:lift_inboundsfw!,:-), (:lift_inboundsbw!,:+,) )
@eval begin
function ($fname){T<:Number}(x::AbstractVector{T}, c::Vector{T}, irange::Range, rhsis::Int)
    nc = length(c)
    if nc == 1  # hard code the most common cases (1, 2, 3) for speed
        c1 = c[1]
        for i in irange
            @inbounds x[i] = ($op)(x[i], c1*x[i+rhsis] )
        end
    elseif nc == 2
        c1,c2 = c[1],c[2]
        rhsisp1 = rhsis+1
        for i in irange
            @inbounds x[i] = ($op)(x[i], c1*x[i+rhsis] + c2*x[i+rhsisp1] )
        end
    elseif nc == 3
        c1,c2,c3 = c[1],c[2],c[3]
        rhsisp1 = rhsis+1
        rhsisp2 = rhsis+2
        for i = irange
            @inbounds x[i] = ($op)(x[i], c1*x[i+rhsis] + c2*x[i+rhsisp1] + c3*x[i+rhsisp2] )
        end
    else
        for i in irange
            for k = 0:nc-1
                @inbounds x[i] = ($op)(x[i], c[k]*x[i+k+rhsis] )
            end
        end
    end
    return x
end
end # eval begin
end # for
