
include("setup_2d.jl")

println("dwt by filtering (N=",N,"x",N,"), ", L, " levels")
f(x0, L, wf) = for i = 1:tn; dwt(x0, wf, L); end
f(x0, L, wf);
@time f(x0, L, wf);

