module Wavelets


include("Util/Util.jl")
include("WT/WT.jl")
include("Transforms/Transforms.jl")
include("Threshold/Threshold.jl")
include("Plot/Plot.jl")

using .Util
using .WT
using .Transforms
using .Threshold
using .Plot


export
    # Submodules (except Util => would lead to conflict with DSP.jl (?)
    WT, Transforms, Threshold, Plot,


    # Util submodule :
    # dyadic
    dyadicdetailindex, dyadicdetailrange, dyadicscalingrange,
    dyadicdetailn, ndyadicscales, maxdyadiclevel,
    tl2dyadiclevel, dyadiclevel2tl,

    # non-dyadic
    detailindex, detailrange, detailn,
    maxtransformlevels, maxmodwttransformlevels,

    # other util functions
    mirror, upsample, downsample, iscube, isdyadic,
    sufficientpoweroftwo, wcount, stridedcopy!,
    isvalidtree, maketree, makewavelet, testfunction,


    # WT submodule :
    DiscreteWavelet, FilterWavelet, LSWavelet, OrthoFilter, GLS, wavelet,


    # Transforms submodule :
    dwt, idwt, dwt!, idwt!, wpt, iwpt, wpt!, iwpt!, modwt, imodwt,


    # Threshold submodule :

    # threshold
    threshold!, threshold, HardTH, SoftTH, SemiSoftTH, SteinTH,
    BiggestTH, PosTH, NegTH,

    # denoising
    DNFT, VisuShrink, denoise, noisest,

    # basis functions
    matchingpursuit, bestbasistree,

    # entropy
    coefentropy, Entropy, ShannonEntropy, LogEnergyEntropy,


    # Plot submodule :
    wplotdots, wplotim


end
