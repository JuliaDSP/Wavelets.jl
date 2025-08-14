module Threshold

using LinearAlgebra: rmul!, norm
using Statistics: median!


include("threshold.jl")
include("basis_functions.jl")
include("denoising.jl")
include("entropy.jl")


export
    # threshold
    threshold!, threshold, HardTH, SoftTH, SemiSoftTH, SteinTH,
    BiggestTH, PosTH, NegTH,

    # denoising
    DNFT, VisuShrink, denoise, noisest,

    # basis functions
    matchingpursuit, bestbasistree,

    # entropy
    coefentropy, Entropy, ShannonEntropy, LogEnergyEntropy

end
