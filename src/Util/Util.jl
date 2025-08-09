module Util

include("non_dyadic.jl")
include("dyadic.jl")
include("util.jl")


export
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

    # core utils functions
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


end # module
