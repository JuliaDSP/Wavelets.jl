
include("setup_2d.jl")

println("fwt (N=",N,"x",N,"), ", L, " levels")
f(x0, L, wf) = for i = 1:tn; fwt(x0, L, wf); end
f(x0, L, wf);
@time f(x0, L, wf);

