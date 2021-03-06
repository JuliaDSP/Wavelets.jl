
##################################################################################
#
#  ORTHOFILTER TRANSFORMS
#  Periodic boundaries, Orthogonal
#
##################################################################################


# DWT
# 1-D
# writes to y
function _dwt!(y::AbstractVector{Ty}, x::AbstractVector{Tx},
                filter::OrthoFilter, L::Integer, fw::Bool) where {Tx<:Number, Ty<:Number}
    T = promote_type(Tx, Ty)
    si = Vector{T}(undef, length(filter)-1) # tmp filter vector
    scfilter, dcfilter = WT.makereverseqmfpair(filter, fw, T)
    return _dwt!(y, x, filter, L, fw, dcfilter, scfilter, si)
end
function _dwt!(y::AbstractVector{<:Number}, x::AbstractVector{<:Number},
                filter::OrthoFilter, L::Integer, fw::Bool,
                dcfilter::Vector{T}, scfilter::Vector{T},
                si::Vector{T}, snew::Vector{T} = Vector{T}(undef, ifelse(L>1, length(x)>>1, 0))) where T<:Number
    n = length(x)
    size(x) == size(y) ||
        throw(DimensionMismatch("in and out array size must match"))
    0 <= L ||
        throw(ArgumentError("L must be positive"))
    sufficientpoweroftwo(y, L) ||
        throw(ArgumentError("size must have a sufficient power of 2 factor"))
    y === x &&
        throw(ArgumentError("in array is out array"))
    L > 1 && length(snew) != n>>1 &&
        throw(ArgumentError("length of snew incorrect"))

    if L == 0
        return copyto!(y,x)
    end
    s = x                           # s is current scaling coefs location
    filtlen = length(filter)

    lrange = 1:L
    !fw && (lrange = reverse(lrange))

    for l in lrange
        if fw
            # detail coefficients
            filtdown!(dcfilter, si, y, detailindex(n,l,1), detailn(n,l), s, 1,-filtlen+1, true)
            # scaling coefficients
            filtdown!(scfilter, si, y, 1, detailn(n,l), s, 1, 0, false)
        else
            # scaling coefficients
            filtup!(false, scfilter, si, y, 1, detailn(n,l-1), s, 1, -filtlen+1, false)
            # detail coefficients
            filtup!(true,  dcfilter, si, y, 1, detailn(n,l-1), x, detailindex(n,l,1), 0, true)
        end
        # if not final iteration: copy to tmp location
        l != lrange[end] && copyto!(snew,1,y,1,detailn(n, fw ? l : l-1))
        L > 1 && (s = snew)
    end
    return y
end
function unsafe_dwt1level!(y::AbstractVector{<:Number}, x::AbstractVector{<:Number},
                            filter::OrthoFilter, fw::Bool,
                            dcfilter::FVector{T}, scfilter::FVector{T},
                            si::FVector{T}) where T<:Number
    n = length(x)
    l = 1
    filtlen = length(filter)

    if fw
        # detail coefficients
        filtdown!(dcfilter, si, y, detailindex(n,l,1), detailn(n,l), x, 1,-filtlen+1, true)
        # scaling coefficients
        filtdown!(scfilter, si, y, 1, detailn(n,l), x, 1, 0, false)
    else
        # scaling coefficients
        filtup!(false, scfilter, si, y, 1, detailn(n,l-1), x, 1, -filtlen+1, false)
        # detail coefficients
        filtup!(true,  dcfilter, si, y, 1, detailn(n,l-1), x, detailindex(n,l,1), 0, true)
    end
    return y
end

function dwt_transform_strided!(y::AbstractArray{<:Number}, x::AbstractArray{<:Number},
                            msub::Int, nsub::Int, stride::Int, idx_func::Function,
                            tmpvec::FVector{T}, tmpvec2::FVector{T},
                            filter::OrthoFilter, fw::Bool,
                            dcfilter::FVector{T}, scfilter::FVector{T}, si::FVector{T}) where T<:Number
    for i=1:msub
        xi = idx_func(i)
        stridedcopy!(tmpvec, x, xi, stride, nsub)
        unsafe_dwt1level!(tmpvec2, tmpvec, filter, fw, dcfilter, scfilter, si)
        stridedcopy!(y, xi, stride, tmpvec2, nsub)
    end
end

function dwt_transform_cols!(y::AbstractArray{<:Number}, x::AbstractArray{<:Number},
                            msub::Int, nsub::Int, idx_func::Function,
                            tmpvec::FVector{T},
                            filter::OrthoFilter, fw::Bool,
                            dcfilter::FVector{T}, scfilter::FVector{T}, si::FVector{T}) where T<:Number
    for i=1:nsub
        xi = idx_func(i)
        copyto!(tmpvec, 1, x, xi, msub)
        ya = unsafe_vectorslice(y, xi, msub)
        unsafe_dwt1level!(ya, tmpvec, filter, fw, dcfilter, scfilter, si)
    end
end

# 2-D
# writes to y
function _dwt!(y::AbstractMatrix{Ty}, x::AbstractMatrix{Tx},
                filter::OrthoFilter, L::Integer, fw::Bool) where {Tx<:Number, Ty<:Number}
    m, n = size(x)
    T = promote_type(Tx, Ty)
    si = Vector{T}(undef, length(filter)-1)       # tmp filter vector
    tmpbuffer = Vector{T}(undef, max(n<<1, m))    # tmp storage vector
    scfilter, dcfilter = WT.makereverseqmfpair(filter, fw, T)

    return _dwt!(y, x, filter, L, fw, dcfilter, scfilter, si, tmpbuffer)
end
function _dwt!(y::AbstractMatrix{<:Number}, x::AbstractMatrix{<:Number},
                filter::OrthoFilter, L::Integer, fw::Bool,
                dcfilter::Vector{T}, scfilter::Vector{T},
                si::Vector{T}, tmpbuffer::Vector{T}) where T<:Number

    m, n = size(x)
    size(x) == size(y) ||
        throw(DimensionMismatch("in and out array size must match"))
    0 <= L ||
        throw(ArgumentError("L must be positive"))
    sufficientpoweroftwo(y, L) ||
        throw(ArgumentError("size must have a sufficient power of 2 factor"))
    y === x &&
        throw(ArgumentError("in array is out array"))
    length(tmpbuffer) >= max(n<<1,m) ||
        throw(ArgumentError("length of tmpbuffer incorrect"))

    if L == 0
        return copyto!(y,x)
    end
    row_stride = m
    #s = x

    if fw
        lrange = 1:L
        nsub = n
        msub = m
    else
        lrange = L:-1:1
        nsub = div(n,2^(L-1))
        msub = div(m,2^(L-1))
        copyto!(y,x)
    end

    inputArray = x

    row_idx_func = i -> row_idx(i, m)
    col_idx_func = i -> col_idx(i, m)
    for l in lrange
        tmpvec  = unsafe_vectorslice(tmpbuffer, 1, msub)
        tmpvec2 = unsafe_vectorslice(tmpbuffer, 1, nsub)
        tmpvec3 = unsafe_vectorslice(tmpbuffer, nsub+1, nsub)
        if fw
            # rows
            dwt_transform_strided!(y, inputArray, msub, nsub, row_stride, row_idx_func,
                tmpvec2, tmpvec3, filter, fw, dcfilter, scfilter, si)
            l == lrange[1] && (inputArray = y)

            # columns
            dwt_transform_cols!(y, y, msub, nsub, col_idx_func,
                tmpvec, filter, fw, dcfilter, scfilter, si)
        else
            # columns
            dwt_transform_cols!(y, inputArray, msub, nsub, col_idx_func,
                tmpvec, filter, fw, dcfilter, scfilter, si)
            l == lrange[1] && (inputArray = y)

            # rows
            dwt_transform_strided!(y, y, msub, nsub, row_stride, row_idx_func,
                tmpvec2, tmpvec3, filter, fw, dcfilter, scfilter, si)
        end
        msub = (fw ? msub>>1 : msub<<1)
        nsub = (fw ? nsub>>1 : nsub<<1)
    end
    return y
end

# 3-D
# writes to y
function _dwt!(y::AbstractArray{Ty, 3}, x::AbstractArray{Tx, 3},
               filter::OrthoFilter, L::Integer, fw::Bool) where {Tx<:Number, Ty<:Number}
    m, n, d = size(x)
    T = promote_type(Tx, Ty)
    si = Vector{T}(undef, length(filter)-1)            # tmp filter vector
    tmpbuffer = Vector{T}(undef, max(m, n<<1, d<<1))   # tmp storage vector
    scfilter, dcfilter = WT.makereverseqmfpair(filter, fw, T)

    return _dwt!(y, x, filter, L, fw, dcfilter, scfilter, si, tmpbuffer)
end
function _dwt!(y::AbstractArray{<:Number, 3}, x::AbstractArray{<:Number, 3},
                filter::OrthoFilter, L::Integer, fw::Bool,
                dcfilter::Vector{T}, scfilter::Vector{T},
                si::Vector{T}, tmpbuffer::Vector{T}) where T<:Number

    m, n, d = size(x)
    size(x) == size(y) ||
        throw(DimensionMismatch("in and out array size must match"))
    0 <= L ||
        throw(ArgumentError("L must be positive"))
    sufficientpoweroftwo(y, L) ||
        throw(ArgumentError("size must have a sufficient power of 2 factor"))
    y === x &&
        throw(ArgumentError("in array is out array"))
    length(tmpbuffer) >= n<<1 ||
        throw(ArgumentError("length of tmpbuffer incorrect"))

    if L == 0
        return copyto!(y,x)
    end
    row_stride = m
    plane_stride = m*n

    if fw
        lrange = 1:L
        msub = m
        nsub = n
        dsub = d
    else
        lrange = L:-1:1
        msub = div(m,2^(L-1))
        nsub = div(n,2^(L-1))
        dsub = div(d,2^(L-1))
        copyto!(y,x)
    end

    inputArray = x

    for l in lrange
        tmpcol  = unsafe_vectorslice(tmpbuffer, 1, msub)
        tmprow  = unsafe_vectorslice(tmpbuffer, 1, nsub)
        tmprow2 = unsafe_vectorslice(tmpbuffer, nsub+1, nsub)
        tmphei  = unsafe_vectorslice(tmpbuffer, 1, dsub)
        tmphei2 = unsafe_vectorslice(tmpbuffer, dsub+1, dsub)
        if fw
            # planes
            for j in 1:nsub
                plane_idx_func = i -> plane_idx(i, j, m)
                dwt_transform_strided!(y, inputArray, msub, dsub, plane_stride, plane_idx_func,
                    tmphei, tmphei2, filter, fw, dcfilter, scfilter, si)
            end
            l == lrange[1] && (inputArray = y)

            # rows
            for j in 1:dsub
                row_idx_func = i -> row_idx(i, j, m, n)
                dwt_transform_strided!(y, y, msub, nsub, row_stride, row_idx_func,
                    tmprow, tmprow2, filter, fw, dcfilter, scfilter, si)
            end
            # columns
            for j in 1:dsub
                col_idx_func = i -> col_idx(i, j, m, n)
                dwt_transform_cols!(y, y, msub, nsub, col_idx_func,
                    tmpcol, filter, fw, dcfilter, scfilter, si)
            end
        else
            # columns
            for j in 1:dsub
                col_idx_func = i -> col_idx(i, j, m, n)
                dwt_transform_cols!(y, inputArray, msub, nsub, col_idx_func,
                    tmpcol, filter, fw, dcfilter, scfilter, si)
            end
            l == lrange[1] && (inputArray = y)

            # rows
            for j in 1:dsub
                row_idx_func = i -> row_idx(i, j, m, n)
                dwt_transform_strided!(y, y, msub, nsub, row_stride, row_idx_func,
                    tmprow, tmprow2, filter, fw, dcfilter, scfilter, si)
            end
            # planes
            for j in 1:nsub
                plane_idx_func = i -> plane_idx(i, j, m)
                dwt_transform_strided!(y, y, msub, dsub, plane_stride, plane_idx_func,
                    tmphei, tmphei2, filter, fw, dcfilter, scfilter, si)
            end
        end
        msub = (fw ? msub>>1 : msub<<1)
        nsub = (fw ? nsub>>1 : nsub<<1)
        dsub = (fw ? dsub>>1 : dsub<<1)
    end
    return y
end



# WPT
# 1-D
# writes to y
function _wpt!(y::AbstractVector{T}, x::AbstractVector{T}, filter::OrthoFilter, tree::BitVector, fw::Bool) where T<:Number
    si = Vector{T}(undef, length(filter)-1)
    ns = ifelse(fw, length(x)>>1, length(x))
    snew = Vector{T}(undef, ns)
    scfilter, dcfilter = WT.makereverseqmfpair(filter, fw, T)

    return _wpt!(y, x, filter, tree, fw, dcfilter, scfilter, si, snew)
end
function _wpt!(y::AbstractVector{T}, x::AbstractVector{T}, filter::OrthoFilter, tree::BitVector, fw::Bool, dcfilter::Vector{T}, scfilter::Vector{T}, si::Vector{T}, snew::Vector{T}) where T<:Number

    size(x) == size(y) ||
        throw(DimensionMismatch("in and out array size must match"))
    y === x &&
        throw(ArgumentError("in array is out array"))
    isvalidtree(y, tree) ||
        throw(ArgumentError("invalid tree"))
    if tree[1] && length(snew) < ifelse(fw, length(x)>>1, length(x))
        throw(ArgumentError("length of snew incorrect"))
    end

    if !tree[1]
        return copyto!(y,x)
    end

    first = true
    n = length(x)
    Lmax = maxtransformlevels(n)
    L = Lmax
    while L > 0
        ix = 1
        k = 1
        Lfw = (fw ? Lmax-L : L-1)
        nj = detailn(n, Lfw)
        treeind = 2^(Lfw)-1
        dx = first ? x : unsafe_vectorslice(snew, 1, nj) # dx will be overwritten if first

        while ix <= n
            if tree[treeind+k]
                dy = unsafe_vectorslice(y, ix, nj)
                if first
                    dx = unsafe_vectorslice(x, ix, nj)
                else
                    copyto!(dx, dy)
                end
                unsafe_dwt1level!(dy, dx, filter, fw, dcfilter, scfilter, si)
            elseif first
                dy = unsafe_vectorslice(y, ix, nj)
                dx = unsafe_vectorslice(x, ix, nj)
                copyto!(dy, dx)
            end
            ix += nj
            k += 1
        end
        L -= 1
        first = false
    end

    return y
end


macro filtermainloop(si, silen, b, val)
    quote
        @inbounds for j=2:$(esc(silen))
            $(esc(si))[j-1] = $(esc(si))[j] + $(esc(b))[j]*$(esc(val))
        end
        @inbounds $(esc(si))[$(esc(silen))] = $(esc(b))[$(esc(silen))+1]*$(esc(val))
    end
end
macro filtermainloopzero(si, silen)
    quote
        @inbounds for j=2:$(esc(silen))
            $(esc(si))[j-1] = $(esc(si))[j]
        end
        @inbounds $(esc(si))[$(esc(silen))] = 0.0
    end
end

# periodic filter and downsampling (by 2)
# apply f to x[ix:ix+nx-1] and write to out[iout:iout+nout-1]
# f : filter
# si: tmp array of length 1 less than f
# out : result gets written to out[iout:iout+nout-1]
# x : filter convolved with x[ix:ix+nx-1], where nx=nout*2 (shifted by shift)
# ss : shift downsampling
# based on Base.filt
function filtdown!(f::AbstractVector{T}, si::AbstractVector{T},
                  out::AbstractVector{<:Number}, iout::Integer, nout::Integer,
                  x::AbstractVector{<:Number}, ix::Integer,
                  shift::Integer=0, ss::Bool=false) where T<:Number
    nx = nout<<1
    silen = length(si)
    flen = length(f)
    @assert length(x) >= ix+nx-1
    @assert length(out) >= iout+nout-1
    @assert (nx-1)>>1 + iout <= length(out)
    @assert silen == flen-1
    @assert shift <= 0

    fill!(si,0.0)
    istart = flen + Int(ss)
    dsshift = (flen%2 + Int(ss))%2  # is flen odd, and shift downsampling

    rout1, rin, rout2 = splitdownrangeper(istart, ix, nx, shift)
    @inbounds begin
        for i in rout1  # rout1 assumed to be in 1:istart-1
            # periodic in the range [ix:ix+nx-1]
            xatind = x[mod(i-1+shift, nx) + ix]
            @filtermainloop(si, silen, f, xatind)
        end
        # rin is inbounds in [ix:ix+nx-1]
        ixsh = -1 + shift + ix
        for i in rin
            xatind = x[i + ixsh]
            # take every other value after istart
            if (i+dsshift)%2 == 0 && i >= istart
                out[(i-istart)>>1 + iout] = si[1] + f[1]*xatind
            end
            @filtermainloop(si, silen, f, xatind)
        end
        for i in rout2
            # periodic in the range [ix:ix+nx-1]
            xatind = x[mod(i-1+shift, nx) + ix]
            # take every other value after istart
            if (i+dsshift)%2 == 0 && i >= istart
                out[(i-istart)>>1 + iout] = si[1] + f[1]*xatind
            end
            @filtermainloop(si, silen, f, xatind)
        end
    end

    return nothing
end
# find part of range which is inbounds for [ix:ix+nx-1] (ixsh = -1 + shift + ix)
# where mod(i-1+shift, nx) + ix == i + ixsh
function splitdownrangeper(istart, ix, nx, shift)
    inxi = 0
    ixsh = -1 + shift + ix
    if mod(shift, nx) + ix == 1 + ixsh   # shift likely 0
        inxi = 1
        iend = nx-1
        while mod(iend - 1 + shift, nx) + ix != iend + ixsh
            iend -= 1
        end
        return (0:-1, inxi:iend, (iend+1):(nx-1+istart))
    elseif mod(istart - 1 + shift, nx) + ix == istart + ixsh
        inxi = istart
        iend = nx-1+istart
        while mod(iend - 1 + shift, nx) + ix != iend + ixsh
            iend -= 1
        end
        return (1:inxi-1, inxi:iend, (iend+1):(nx-1+istart))
    else
        return (0:-1, 0:-1, 1:(nx-1+istart))
    end
end


# periodic filter and upsampling (by 2)
# apply f to x[ix:ix+nx-1] upsampled and write to out[iout:iout+nout-1]
# f : filter
# si: tmp array of length 1 less than f
# out : result gets written to out[iout:iout+nout-1]
# x : filter convolved with x[ix:ix+nx-1] upsampled, where nout==nx*2 (then shifted by shift)
# ss : shift upsampling
# based on Base.filt
function filtup!(add2out::Bool, f::Vector{T}, si::Vector{T},
              out::AbstractVector{<:Number}, iout::Integer, nout::Integer,
              x::AbstractVector{<:Number}, ix::Integer,
              shift::Integer=0, ss::Bool=false) where T<:Number
    nx = nout>>1
    silen = length(si)
    flen = length(f)
    @assert length(x) >= ix+nx-1                # check array size
    @assert length(out) >= iout+nout-1          # check array size
    @assert (nx<<1-1)>>1 + iout <= length(out)  # max for out index
    @assert silen == flen-1
    @assert shift <= 0

    fill!(si,0.0)
    istart = flen - shift%2
    dsshift = Int(ss)%2  # shift upsampling

    rout1, rin, rout2 = splituprangeper(istart, ix, nx, nout, shift)
    @inbounds begin
        xatind = 0.0
        for i in rout1  # rout1 assumed to be in 1:istart-1
            if (i+dsshift)%2 == 0
                @filtermainloopzero(si, silen)
            else
                # periodic in the range [ix:ix+nx-1]
                xindex = mod((i-1)>>1+shift>>1,nx) + ix    #(i-1)>>1 increm. every other
                xatind = x[xindex]
                @filtermainloop(si, silen, f, xatind)
            end
        end
        ixsh = shift>>1 + ix
        for i in rin
            if (i+dsshift)%2 == 0
                xatind = 0.0
            else
                xindex = (i-1)>>1 + ixsh
                xatind = x[xindex]
            end
            if i >= istart
                if add2out
                    out[(i-istart) + iout] += si[1] + f[1]*xatind
                else
                    out[(i-istart) + iout] = si[1] + f[1]*xatind
                end
            end
            if (i+dsshift)%2 == 0
                @filtermainloopzero(si, silen)
            else
                @filtermainloop(si, silen, f, xatind)
            end
        end
        for i in rout2
            if (i+dsshift)%2 == 0
                xatind = 0.0
            else
                xindex = mod((i-1)>>1+shift>>1,nx) + ix
                xatind = x[xindex]
            end
            if i >= istart
                if add2out
                    out[(i-istart) + iout] += si[1] + f[1]*xatind
                else
                    out[(i-istart) + iout] = si[1] + f[1]*xatind
                end
            end
            if (i+dsshift)%2 == 0
                @filtermainloopzero(si, silen)
            else
                @filtermainloop(si, silen, f, xatind)
            end
        end
    end

    return nothing
end
# find part of range which is inbounds for [ix:ix+nx-1] (ixsh = shift>>1 + ix)
# where mod((i-1)>>1+shift>>1, nx) + ix == (i-1)>>1 + ixsh
function splituprangeper(istart, ix, nx, nout, shift)
    inxi = 0
    ixsh = shift>>1 + ix
    if mod(shift>>1, nx) + ix == ixsh   # shift likely 0
        inxi = 1
        iend = nout-1
        while mod((iend-1)>>1+shift>>1, nx) + ix != (iend-1)>>1 + ixsh
            iend -= 1
        end
        return (0:-1, inxi:iend, (iend+1):(nout-1+istart))
    elseif mod((istart-1)>>1+shift>>1, nx) + ix == (istart-1)>>1 + ixsh
        inxi = istart
        iend = nout-1+istart
        while mod((iend-1)>>1+shift>>1, nx) + ix != (iend-1)>>1 + ixsh
            iend -= 1
        end
        return (1:inxi-1, inxi:iend, (iend+1):(nout-1+istart))
    else
        return (0:-1, 0:-1, 1:(nout-1+istart))
    end
end
