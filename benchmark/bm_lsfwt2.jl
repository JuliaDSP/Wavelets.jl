
include("setup_2d.jl")

println("fwt by lifting (N=",N,"x",N,"), ", L, " levels")
f(x0, L, wl) = for i = 1:tn; fwt(x0, L, wl); end
f(x0, L, wl);
@time f(x0, L, wl);

#@profile for i=1:tn dwt!(x0, L, wl, true); end
#@profile for i=1:tn fwt(x0, L, wl); end
#Profile.print()

