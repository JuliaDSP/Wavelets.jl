
include("setup_1d.jl")

println("idwt by filtering (N=",N,"), ", L, " levels")
f(x0, L, wf) = for i = 1:tn; idwt(x0, wf, L); end
f(x0, L, wf);
@time f(x0, L, wf);

