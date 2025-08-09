module Wavelets

# Include all submodules
include("Util/Util.jl")
include("WT/WT.jl")
include("Transforms/Transforms.jl")
include("Threshold/Threshold.jl")
include("Plot/Plot.jl")

# Re-export everything (no prefixes required)
using .Util
using .Threshold
using .Transforms
using .WT
using .Plot

using SpecialFunctions
# TODO : no "using" of external packages in module declaration


export
    # From Util
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
    testfunction,
    # From Threshold
    # threshold
    threshold!,
    threshold,
    HardTH,
    SoftTH,
    SemiSoftTH,
    SteinTH,
    BiggestTH,
    PosTH,
    NegTH,
    # denoising
    DNFT,
    VisuShrink,
    denoise,
    noisest,
    # basis functions
    matchingpursuit,
    bestbasistree,
    # entropy
    coefentropy,
    Entropy,
    ShannonEntropy,
    LogEnergyEntropy,
    # From Transforms
    dwt,
    idwt,
    dwt!,
    idwt!,
    wpt,
    iwpt,
    wpt!,
    iwpt!,
    modwt,
    imodwt,
    # From WT
    DiscreteWavelet,
    FilterWavelet,
    LSWavelet,
    OrthoFilter,
    GLS,
    wavelet,
    # From Plot
    wplotdots,
    wplotim



end # module Wavelets
