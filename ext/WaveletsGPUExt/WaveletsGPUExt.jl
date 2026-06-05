module WaveletsGPUExt

using Wavelets
using GPUArrays
using GPUArrays: AbstractGPUVector, AbstractGPUMatrix, AbstractGPUArray
using KernelAbstractions
using KernelAbstractions: Adapt

import Wavelets: WT, Util
using Wavelets.WT: OrthoFilter, GLS
import Wavelets.Transforms: _dwt!, _wpt!, modwt, imodwt, modwt_step, imodwt_step
import Wavelets.Util: detailindex, detailn, detailrange,
    maxtransformlevels, maxmodwttransformlevels,
    sufficientpoweroftwo,
    isvalidtree, maketree, iscube

include("utils_gpu.jl")
include("filter_transforms_gpu.jl")
include("lifting_transforms_gpu.jl")
include("modwt_gpu.jl")

end # module
