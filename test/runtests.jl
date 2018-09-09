using Wavelets
using Compat.Test
using Compat.LinearAlgebra
using Compat.DelimitedFiles
using Compat: ComplexF64, copyto!, range
import Compat

# modified from Base.Test
function vecnorm_eq(va, vb, Eps, astr="a", bstr="b")
    if length(va) != length(vb)
        #error("lengths of ", astr, " and ", bstr, " do not match: ",
        #      "\n  ", astr, " (length $(length(va))) = ", va,
        #      "\n  ", bstr, " (length $(length(vb))) = ", vb)
        return false
    end

    diff = Compat.norm(va - vb)
    if !isnan(Eps) && diff > Eps
        #sdiff = string("|", astr, " - ", bstr, "| <= ", Eps)
        #error("assertion failed: ", sdiff,
        #      "\n  difference = ", diff, " > ", Eps)
        return false
    end
    return true
end
macro vecnorm_eq_eps(a, b, c)
    :(vecnorm_eq($(esc(a)), $(esc(b)), $(esc(c)), $(string(a)), $(string(b))))
end

@testset "Util" begin include("util.jl") end
@testset "Transforms" begin include("transforms.jl") end
@testset "Threshold" begin include("threshold.jl") end
@testset "Plot" begin include("plot.jl") end
