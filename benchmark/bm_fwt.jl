
include("setup_1d.jl")

println("fwt by filtering (N=",N,"), ", L, " levels")
f(x0, L, wf) = for i = 1:tn; fwt(x0, L, wf); end
f(x0, L, wf);
@time f(x0, L, wf);

