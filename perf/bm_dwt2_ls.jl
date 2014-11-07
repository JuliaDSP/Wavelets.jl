
include("setup_2d.jl")

println("dwt by lifting (N=",N,"x",N,"), ", L, " levels")
f(x0, L, wl) = for i = 1:tn; dwt(x0, wl, L); end
f(x0, L, wl);
@time f(x0, L, wl);

#@profile for i=1:tn dwt!(x0, L, wl, true); end
#let
#    @profile dwt(x0, wl, L);
#    Profile.print()
#end

