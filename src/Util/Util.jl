module Util

using LinearAlgebra: rmul!, norm
using DSP: conv



include("non_dyadic.jl")
include("dyadic.jl")
include("_util.jl")



export
    # dyadic
    dyadicdetailindex,
    dyadicdetailrange,
    dyadicscalingrange,
    dyadicdetailn,
    ndyadicscales,
    maxdyadiclevel,
    tl2dyadiclevel,
    dyadiclevel2tl,
    # non-dyadic
    detailindex,
    detailrange,
    detailn,
    maxtransformlevels,
    maxmodwttransformlevels,
    # other util functions
    mirror,
    upsample,
    downsample,
    iscube,
    isdyadic,
    sufficientpoweroftwo,
    wcount,
    stridedcopy!,
    isvalidtree,
    maketree,
    makewavelet,
    testfunction

end
