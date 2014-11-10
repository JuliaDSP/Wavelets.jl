
include("setup_1d.jl")

println("idwt by lifting (N=",N,"), ", L, " levels")
f(x0, L, wl) = for i = 1:tn; idwt(x0, wl, L); end
f(x0, L, wl);
@time f(x0, L, wl);

#@profile f(x0,L,wl);
#Profile.print()

