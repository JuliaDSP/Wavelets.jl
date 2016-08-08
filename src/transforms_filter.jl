
##################################################################################
#
#  ORTHOFILTER TRANSFORMS
#  Periodic boundaries, Orthogonal
#
##################################################################################


# DWT
# 1-D
# writes to y
function _dwt!{T<:Number}(y::AbstractVector{T}, x::AbstractVector{T},
                                filter::OrthoFilter, L::Integer, fw::Bool)
    si = Array(T, length(filter)-1)       # tmp filter vector
    scfilter, dcfilter = WT.makereverseqmfpair(filter, fw, T)
    return _dwt!(y, x, filter, L, fw, dcfilter, scfilter, si)
end
function _dwt!{T<:Number}(y::AbstractVector{T}, x::AbstractVector{T},
                                filter::OrthoFilter, L::Integer, fw::Bool,
                                dcfilter::Vector{T}, scfilter::Vector{T},
                                si::Vector{T}, snew::Vector{T} = Array(T, ifelse(L>1, length(x)>>1, 0)))
    n = length(x)
    size(x) == size(y) ||
        throw(DimensionMismatch("in and out array size must match"))
    0 <= L ||
        throw(ArgumentError("L must be positive"))
    sufficientpoweroftwo(y, L) ||
        throw(ArgumentError("size must have a sufficient power of 2 factor"))
    is(y,x) &&
        throw(ArgumentError("in array is out array"))
    L > 1 && length(snew) != n>>1 &&
        throw(ArgumentError("length of snew incorrect"))

    if L == 0
        return copy!(y,x)
    end
    s = x                           # s is currect scaling coefs location
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
        l != lrange[end] && copy!(snew,1,y,1,detailn(n, fw ? l : l-1))
        L > 1 && (s = snew)
    end
    return y
end
function unsafe_dwt1level!{T<:Number}(y::AbstractVector{T}, x::AbstractVector{T},
                                            filter::OrthoFilter, fw::Bool,
                                            dcfilter::Vector{T}, scfilter::Vector{T},
                                            si::Vector{T})
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


function dwt_transform_strided!{T<:Number}(y::Array{T}, x::AbstractArray{T},
                            nsub::Int, stride::Int, idx_func::Function,
                            tmpvec::Vector{T}, tmpvec2::Vector{T},
                            filter::OrthoFilter, fw::Bool,
                            dcfilter::Vector{T}, scfilter::Vector{T}, si::Vector{T})
    for i=1:nsub
        xi = idx_func(i)
        stridedcopy!(tmpvec, x, xi, stride, nsub)
        unsafe_dwt1level!(tmpvec2, tmpvec, filter, fw, dcfilter, scfilter, si)
        stridedcopy!(y, xi, stride, tmpvec2, nsub)
    end
end

function dwt_transform_cols!{T<:Number}(y::Array{T}, x::AbstractArray{T},
                            nsub::Int, idx_func::Function,
                            tmpvec::Vector{T},
                            filter::OrthoFilter, fw::Bool,
                            dcfilter::Vector{T}, scfilter::Vector{T}, si::Vector{T})
    for i=1:nsub
        xi = idx_func(i)
        copy!(tmpvec, 1, x, xi, nsub)
        ya = unsafe_vectorslice(y, xi, nsub)
        unsafe_dwt1level!(ya, tmpvec, filter, fw, dcfilter, scfilter, si)
    end
end

# 2-D
# writes to y
function _dwt!{T<:Number}(y::Matrix{T}, x::AbstractMatrix{T},
                                filter::OrthoFilter, L::Integer, fw::Bool)
    n = size(x,1)
    si = Array(T, length(filter)-1)       # tmp filter vector
    tmpbuffer = Array(T,n<<1)             # tmp storage vector
    scfilter, dcfilter = WT.makereverseqmfpair(filter, fw, T)

    return _dwt!(y, x, filter, L, fw, dcfilter, scfilter, si, tmpbuffer)
end
function _dwt!{T<:Number}(y::Matrix{T}, x::AbstractMatrix{T},
                                filter::OrthoFilter, L::Integer, fw::Bool,
                                dcfilter::Vector{T}, scfilter::Vector{T},
                                si::Vector{T}, tmpbuffer::Vector{T})

    n = size(x,1)
    size(x) == size(y) ||
        throw(DimensionMismatch("in and out array size must match"))
    iscube(x) ||
        throw(ArgumentError("array must be square/cube"))
    0 <= L ||
        throw(ArgumentError("L must be positive"))
    sufficientpoweroftwo(y, L) ||
        throw(ArgumentError("size must have a sufficient power of 2 factor"))
    is(y,x) &&
        throw(ArgumentError("in array is out array"))
    length(tmpbuffer) >= n<<1 ||
        throw(ArgumentError("length of tmpbuffer incorrect"))

    if L == 0
        return copy!(y,x)
    end
    row_stride = n
    #s = x

    if fw
        lrange = 1:L
        nsub = n
    else
        lrange = L:-1:1
        nsub = div(n,2^(L-1))
        copy!(y,x)
    end

    inputArray = x

    for l in lrange
        tmpvec = unsafe_vectorslice(tmpbuffer, 1, nsub)
        tmpvec2 = unsafe_vectorslice(tmpbuffer, nsub+1, nsub)
        row_idx_func = i -> row_idx(i, n)
        col_idx_func = i -> col_idx(i, n)
        if fw
            # rows
            dwt_transform_strided!(y, inputArray, nsub, row_stride, row_idx_func,
                tmpvec, tmpvec2, filter, fw, dcfilter, scfilter, si)
            l == lrange[1] && (inputArray = y)

            # columns
            dwt_transform_cols!(y, y, nsub, col_idx_func,
                tmpvec, filter, fw, dcfilter, scfilter, si)
        else
            # columns
            dwt_transform_cols!(y, inputArray, nsub, col_idx_func,
                tmpvec, filter, fw, dcfilter, scfilter, si)
            l == lrange[1] && (inputArray = y)

            # rows
            dwt_transform_strided!(y, y, nsub, row_stride, row_idx_func,
                tmpvec, tmpvec2, filter, fw, dcfilter, scfilter, si)
        end
        nsub = (fw ? nsub>>1 : nsub<<1)
    end
    return y
end

# 3-D
# writes to y
function _dwt!{T<:Number}(y::Array{T, 3}, x::AbstractArray{T, 3},
                                filter::OrthoFilter, L::Integer, fw::Bool)
    n = size(x,1)
    si = Array(T, length(filter)-1)       # tmp filter vector
    tmpbuffer = Array(T,n<<1)             # tmp storage vector
    scfilter, dcfilter = WT.makereverseqmfpair(filter, fw, T)

    return _dwt!(y, x, filter, L, fw, dcfilter, scfilter, si, tmpbuffer)
end
function _dwt!{T<:Number}(y::Array{T, 3}, x::AbstractArray{T, 3},
                                filter::OrthoFilter, L::Integer, fw::Bool,
                                dcfilter::Vector{T}, scfilter::Vector{T},
                                si::Vector{T}, tmpbuffer::Vector{T})

    n = size(x,1)
    size(x) == size(y) ||
        throw(DimensionMismatch("in and out array size must match"))
    iscube(x) ||
        throw(ArgumentError("array must be square/cube"))
    0 <= L ||
        throw(ArgumentError("L must be positive"))
    sufficientpoweroftwo(y, L) ||
        throw(ArgumentError("size must have a sufficient power of 2 factor"))
    is(y,x) &&
        throw(ArgumentError("in array is out array"))
    length(tmpbuffer) >= n<<1 ||
        throw(ArgumentError("length of tmpbuffer incorrect"))

    if L == 0
        return copy!(y,x)
    end
    row_stride = n
    height_stride = n*n

    if fw
        lrange = 1:L
        nsub = n
    else
        lrange = L:-1:1
        nsub = div(n,2^(L-1))
        copy!(y,x)
    end

    inputArray = x

    for l in lrange
        tmpvec = unsafe_vectorslice(tmpbuffer, 1, nsub)
        tmpvec2 = unsafe_vectorslice(tmpbuffer, nsub+1, nsub)
        if fw
            # heights
            for j in 1:nsub
                hei_idx_func = i -> hei_idx(i, j, n)
                dwt_transform_strided!(y, inputArray, nsub, height_stride, hei_idx_func,
                    tmpvec, tmpvec2, filter, fw, dcfilter, scfilter, si)
            end
            l == lrange[1] && (inputArray = y)

            # rows
            for j in 1:nsub
                row_idx_func = i -> row_idx(i, j, n)
                dwt_transform_strided!(y, y, nsub, row_stride, row_idx_func,
                    tmpvec, tmpvec2, filter, fw, dcfilter, scfilter, si)
            end
            # columns
            for j in 1:nsub
                col_idx_func = i -> col_idx(i, j, n)
                dwt_transform_cols!(y, y, nsub, col_idx_func,
                    tmpvec, filter, fw, dcfilter, scfilter, si)
            end
        else
            # columns
            for j in 1:nsub
                col_idx_func = i -> col_idx(i, j, n)
                dwt_transform_cols!(y, inputArray, nsub, col_idx_func,
                    tmpvec, filter, fw, dcfilter, scfilter, si)
            end
            l == lrange[1] && (inputArray = y)

            # rows
            for j in 1:nsub
                row_idx_func = i -> row_idx(i, j, n)
                dwt_transform_strided!(y, y, nsub, row_stride, row_idx_func,
                    tmpvec, tmpvec2, filter, fw, dcfilter, scfilter, si)
            end
            # heights
            for j in 1:nsub
                hei_idx_func = i -> hei_idx(i, j, n)
                dwt_transform_strided!(y, y, nsub, height_stride, hei_idx_func,
                    tmpvec, tmpvec2, filter, fw, dcfilter, scfilter, si)
            end
        end
        nsub = (fw ? nsub>>1 : nsub<<1)
    end
    return y
end



# WPT
# 1-D
# writes to y
function _wpt!{T<:Number}(y::AbstractVector{T}, x::AbstractVector{T}, filter::OrthoFilter, tree::BitVector, fw::Bool)
    si = Array(T, length(filter)-1)
    ns = ifelse(fw, length(x)>>1, length(x))
    snew = Array(T, ns)
    scfilter, dcfilter = WT.makereverseqmfpair(filter, fw, T)

    return _wpt!(y, x, filter, tree, fw, dcfilter, scfilter, si, snew)
end
function _wpt!{T<:Number}(y::AbstractVector{T}, x::AbstractVector{T}, filter::OrthoFilter, tree::BitVector, fw::Bool, dcfilter::Vector{T}, scfilter::Vector{T}, si::Vector{T}, snew::Vector{T})

    size(x) == size(y) ||
        throw(DimensionMismatch("in and out array size must match"))
    is(y,x) &&
        throw(ArgumentError("in array is out array"))
    isvalidtree(y, tree) ||
        throw(ArgumentError("invalid tree"))
    if tree[1] && length(snew) < ifelse(fw, length(x)>>1, length(x))
        throw(ArgumentError("length of snew incorrect"))
    end

    if !tree[1]
        return copy!(y,x)
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
        dx = unsafe_vectorslice(snew, 1, nj)

        while ix <= n
            if tree[treeind+k]
                dy = unsafe_vectorslice(y, ix, nj)
                if first
                    dx = unsafe_vectorslice(x, ix, nj)
                else
                    copy!(dx, dy)
                end
                unsafe_dwt1level!(dy, dx, filter, fw, dcfilter, scfilter, si)
            elseif first
                dy = unsafe_vectorslice(y, ix, nj)
                dx = unsafe_vectorslice(x, ix, nj)
                copy!(dy, dx)
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
function filtdown!{T<:Number}(f::Vector{T}, si::Vector{T},
                              out::AbstractVector{T},  iout::Integer, nout::Integer,
                              x::AbstractVector{T}, ix::Integer, shift::Integer=0, ss::Bool=false)
    nx = nout<<1
    silen = length(si)
    flen = length(f)
    @assert length(x) >= ix+nx-1
    @assert length(out) >= iout+nout-1
    @assert (nx-1)>>1 + iout <= length(out)
    @assert silen == flen-1
    @assert shift <= 0

    fill!(si,0.0)
    istart = @compat flen + Int(ss)
    dsshift = @compat (flen%2 + Int(ss))%2  # is flen odd, and shift downsampling

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
function filtup!{T<:Number}(add2out::Bool, f::Vector{T}, si::Vector{T},
                              out::AbstractVector{T},  iout::Integer, nout::Integer,
                              x::AbstractVector{T}, ix::Integer, shift::Integer=0, ss::Bool=false)
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
    dsshift = @compat Int(ss)%2  # shift upsampling

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
