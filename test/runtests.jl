using Wavelets
using Base.Test

# v0.2
if !isdefined(:vecnorm)
    function vecnorm(s)
        norm(s[:])
    end
end

include("util.jl")
include("transforms.jl")


