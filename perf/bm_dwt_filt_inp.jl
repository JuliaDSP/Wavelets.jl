
include("setup_1d.jl")

y = copy(x0)

println("dwt! by filtering (N=",N,"), ", L, " levels")
f(y, x0, L, wf) = for i = 1:tn; dwt!(y, x0, wf, L); end
f(y, x0, L, wf);
@time f(y, x0, L, wf);

