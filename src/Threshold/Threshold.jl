module Threshold

include("threshold.jl")
include("entropy.jl")
include("denoising.jl")
include("basis_functions.jl")

export
    # threshold functions and types
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
    LogEnergyEntropy


end
