
##################################################################################
#
#  LIFTING TRANSFORMS 
#  Periodic boundaries, dyadic length (powers of 2)
#
##################################################################################

for (Xwt) in (:dwt!, :wpt!)
@eval begin
    # pseudo "out of place" by copying
    function $Xwt{T<:FloatingPoint}(y::AbstractArray{T}, x::AbstractArray{T}, scheme::GLS, L::Union(Integer, BitVector), fw::Bool)
        copy!(y, x)
        $Xwt(y, scheme, L, fw)
        return y
    end
end
end

# return scheme parameters adjusted for direction and type
function makescheme{T<:FloatingPoint}(::Type{T}, scheme::GLS, fw::Bool)
    if fw
        return scheme.step, convert(T, scheme.norm1), convert(T, scheme.norm2)
    else
        return reverse(scheme.step), convert(T, 1/scheme.norm1), convert(T, 1/scheme.norm2)
    end
end

# 1-D
# inplace transform of y, no vector allocation
# tmp: size at least n>>2
function dwt!{T<:FloatingPoint}(y::AbstractVector{T}, scheme::GLS, L::Integer, fw::Bool, tmp::Vector{T}=Array(T,length(y)>>2))

    n = length(y)
    #J = nscales(n)
    #@assert isdyadic(y)
    @assert sufficientpowersoftwo(y,L)
    @assert 0 <= L #<= J
    @assert length(tmp) >= n>>2
    L == 0 && return y          # do nothing
    
    if fw
        #jrange = (J-1):-1:(J-L)
        lrange = 1:L
        ns = n
        half = ns>>1
    else
        #jrange = (J-L):(J-1)
        lrange = L:-1:1
        #ns = 2^(jrange[1]+1)
        ns = div(n,(2^(lrange[1]-1)))
        half = ns>>1
    end
    s = y
    stepseq, norm1, norm2 = makescheme(T, scheme, fw)

    #for j in jrange
    for l in lrange
        if fw
            split!(s, ns, tmp)
            for step in stepseq
                stepcoef = convert(Array{T}, step.coef)
                if step.stept == 'p'
                    predictfw!(s, half, stepcoef, step.shift)
                elseif step.stept == 'u'
                    updatefw!(s, half, stepcoef, step.shift)
                end
            end
            normalize!(s, half, ns, norm1, norm2)
            ns = ns>>1 
            half = half>>1
        else
            normalize!(s, half, ns, norm1, norm2)
            for step in stepseq
                stepcoef = convert(Array{T}, step.coef)
                if step.stept == 'p'
                    predictbw!(s, half, stepcoef, step.shift)
                elseif step.stept == 'u'
                    updatebw!(s, half, stepcoef, step.shift)
                end
            end
            merge!(s, ns, tmp)        # inverse split
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
function unsafe_dwt1level!{T<:FloatingPoint}(y::AbstractArray{T}, iy::Integer, incy::Integer, oopc::Bool, oopv::Vector{T}, scheme::GLS, fw::Bool, stepseq::Vector, norm1::T, norm2::T, tmp::Vector{T})
    if !oopc
        oopv = y
    end
    ns = length(oopv)
    half = ns>>1

    if fw
        if oopc
            split!(oopv, y, iy, incy, ns)
        else
            split!(oopv, ns, tmp)
        end
        for step in stepseq
            stepcoef = convert(Vector{T}, step.coef)
            if step.stept == 'p'
                predictfw!(oopv, half, stepcoef, step.shift)
            else # 'u'
                updatefw!(oopv, half, stepcoef, step.shift)
            end
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
            stepcoef = convert(Vector{T}, step.coef)
            if step.stept == 'p'
                predictbw!(oopv, half, stepcoef, step.shift)
            else # 'u'
                updatebw!(oopv, half, stepcoef, step.shift)
            end
        end
        if oopc
            merge!(y, iy, incy, oopv, ns)
        else
            merge!(oopv, ns, tmp)
        end
    end

    return y
end

# 2-D
# inplace transform of y, no vector allocation
# tmp: size at least n>>2
# tmpvec: size at least n
function dwt!{T<:FloatingPoint}(y::Matrix{T}, scheme::GLS, L::Integer, fw::Bool, tmp::Vector{T}=Array(T,size(y,1)>>2), tmpvec::Vector{T}=Array(T,size(y,1)))

    n = size(y,1)
    #J = nscales(n)
    @assert iscube(y)
    #@assert isdyadic(y)
    @assert sufficientpowersoftwo(y,L)
    @assert 0 <= L #<= J
    @assert length(tmp) >= n>>2
    @assert length(tmpvec) >= n
    L == 0 && return y          # do nothing
    
    if fw
        #jrange = (J-1):-1:(J-L)
        lrange = 1:L
        nsub = n
    else
        #jrange = (J-L):(J-1)
        lrange = L:-1:1
        #nsub = int(2^(J-L+1))
        nsub = div(n,2^(L-1))
    end
    stepseq, norm1, norm2 = makescheme(T, scheme, fw)
    
    xm = 0
    xs = n
    #for j in jrange
    for l in lrange
        tmpsub = unsafe_vectorslice(tmpvec, 1, nsub)
        if fw
            # rows
            for i=1:nsub
                xi = i
                # out of place in a dense array for speed
                unsafe_dwt1level!(y, xi, xs, true, tmpsub, scheme, fw, 
                                    stepseq, norm1, norm2, tmp)
            end
            # columns
            for i=1:nsub
                xi = 1+(i-1)*n
                ya = unsafe_vectorslice(y, xi, nsub)
                unsafe_dwt1level!(ya, 1, 1, false, tmpsub, scheme, fw, 
                                    stepseq, norm1, norm2, tmp)
            end       
        else
            # columns
            for i=1:nsub
                xi = 1+(i-1)*n
                ya = unsafe_vectorslice(y, xi, nsub)
                unsafe_dwt1level!(ya, 1, 1, false, tmpsub, scheme, fw, 
                                    stepseq, norm1, norm2, tmp)
            end   
            # rows
            for i=1:nsub
                xi = i
                # out of place in a dense array for speed
                unsafe_dwt1level!(y, xi, xs, true, tmpsub, scheme, fw, 
                                    stepseq, norm1, norm2, tmp)
            end
        end 

        fw  && (nsub = nsub>>1)
        !fw && (nsub = nsub<<1)
    end
    
    return y
end



# WPT
# 1-D
function wpt!{T<:FloatingPoint}(y::AbstractVector{T}, scheme::GLS, L::Integer, fw::Bool, tmp::Vector{T}=Array(T,length(y)>>2))
    wpt!(y, scheme, maketree(length(y), L, :full), fw, tmp)
end
function wpt!{T<:FloatingPoint}(y::AbstractVector{T}, scheme::GLS, tree::BitVector, fw::Bool, tmp::Vector{T}=Array(T,length(y)>>2))

    n = length(y)
    J = nscales(n)
    @assert isdyadic(y)
    @assert isvalidtree(y, tree)
    @assert length(tmp) >= n>>2
    tree[1] || return y          # do nothing
    
    stepseq, norm1, norm2 = makescheme(T, scheme, fw)
    
    L = J
    while L > 0
        ix = 1
        k = 1
        fw  && (Lfw = J-L)
        !fw && (Lfw = L-1)
        nj = detailn(tl2level(n, Lfw))
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



function normalize!{T<:FloatingPoint}(x::AbstractVector{T}, half::Int, ns::Int, n1::T, n2::T)
    for i = 1:half
        @inbounds x[i] *= n1
    end
    for i = half+1:ns
        @inbounds x[i] *= n2
    end
    return x
end
# out of place normalize from x to y
function normalize!{T<:FloatingPoint}(y::AbstractVector{T}, x::AbstractVector{T}, half::Int, ns::Int, n1::T, n2::T)
    for i = 1:half
        @inbounds y[i] = n1*x[i]
    end
    for i = half+1:ns
        @inbounds y[i] = n2*x[i]
    end
    return y
end
function normalize!{T<:FloatingPoint}(y::AbstractArray{T}, iy::Int, incy::Int, x::AbstractVector{T}, half::Int, ns::Int, n1::T, n2::T)
    for i = 1:half
        @inbounds y[iy + (i-1)*incy] = n1*x[i]
    end
    for i = half+1:ns
        @inbounds y[iy + (i-1)*incy] = n2*x[i]
    end
    return y
end
function normalize!{T<:FloatingPoint}(y::AbstractVector{T}, x::AbstractArray{T}, ix::Int, incx::Int, half::Int, ns::Int, n1::T, n2::T)
    for i = 1:half
        @inbounds y[i] = n1*x[ix + (i-1)*incx]
    end
    for i = half+1:ns
        @inbounds y[i] = n2*x[ix + (i-1)*incx]
    end
    return y
end


# predict and update lifting steps inplace on x, forward and backward
# half: half of the length under consideration, shift: shift to left, c: lift coefs
# For predict: writes to range 1:half, reads from 1:2*half
# For update : writes to range half+1:2*half, reads from 1:2*half
for (fname,lift_inb,lift_bound,pred) in (
                        (:predictfw!,:lift_inboundsfw!, :liftp_perboundaryfw!, true),
                        (:predictbw!,:lift_inboundsbw!, :liftp_perboundarybw!, true),
                        (:updatefw!, :lift_inboundsfw!, :liftu_perboundaryfw!, false),
                        (:updatebw!, :lift_inboundsbw!, :liftu_perboundarybw!, false) )
@eval begin
function ($fname){T<:FloatingPoint}(x::AbstractVector{T}, half::Int, c::Vector{T}, shift::Int)
    lhsr, irange, rhsr, rhsis = getliftranges(half, length(c), shift, $pred)
    
    # left boundary
    ($lift_bound)(x, half, c, lhsr, rhsis)
    # main loop
    ($lift_inb)(x, c, irange, rhsis)
    # right boundary
    ($lift_bound)(x, half, c, rhsr, rhsis)
    return x
end
end # eval begin
end # for

function getliftranges(half::Int, nc::Int, shift::Int, pred::Bool)
    # define index shift rhsis
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
    if !(pred)  # shift ranges for update
        irange += half
        lhsr += half
        rhsr += half
    end
    return (lhsr, irange, rhsr, rhsis)
end

# periodic boundary
for (fname,op,puxind) in (  (:liftp_perboundaryfw!,:-,:(mod1(i+k-1+rhsis-half,half)+half)),
                            (:liftp_perboundarybw!,:+,:(mod1(i+k-1+rhsis-half,half)+half)),
                            (:liftu_perboundaryfw!, :-,:(mod1(i+k-1+rhsis,half))),
                            (:liftu_perboundarybw!, :+,:(mod1(i+k-1+rhsis,half)))
                            )
@eval begin
function ($fname){T<:FloatingPoint}(x::AbstractVector{T}, half::Int, c::Vector{T}, irange::Range, rhsis::Int)
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
function ($fname){T<:FloatingPoint}(x::AbstractVector{T}, c::Vector{T}, irange::Range, rhsis::Int)
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



