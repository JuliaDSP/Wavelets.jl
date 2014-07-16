
include("setup_1d.jl")

println("fwt by lifting (N=",N,"), ", L, " levels")
f(x0, L, wl) = for i = 1:tn; fwt(x0, L, wl); end
f(x0, L, wl);
@time f(x0, L, wl);

#@profile f(x0,L,wl);
#Profile.print()

