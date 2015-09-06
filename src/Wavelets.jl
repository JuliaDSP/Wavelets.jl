VERSION >= v"0.4.0-dev+6521" && __precompile__()

module Wavelets

include("util.jl")

include("wt.jl")
include("transforms.jl")

include("threshold.jl")
include("plot.jl")

using Reexport
@reexport using .Util, .WT, .Transforms, .Threshold, .Plot

end

