using Wavelets
using Base.Test

# modified from Base.Test
function test_vecnorm_eq(va, vb, Eps, astr, bstr)
    if length(va) != length(vb)
        error("lengths of ", astr, " and ", bstr, " do not match: ",
              "\n  ", astr, " (length $(length(va))) = ", va,
              "\n  ", bstr, " (length $(length(vb))) = ", vb)
    end
    diff = vecnorm(va - vb)

    if !isnan(Eps) && !(diff <= Eps)
        sdiff = string("|", astr, " - ", bstr, "| <= ", Eps)
        error("assertion failed: ", sdiff,
              "\n  difference = ", diff, " > ", Eps)
    end
end
macro test_vecnorm_eq_eps(a, b, c)
    :(test_vecnorm_eq($(esc(a)), $(esc(b)), $(esc(c)), $(string(a)), $(string(b))))
end



include("util.jl")
include("transforms.jl")

print("\ntesting: success\n")

