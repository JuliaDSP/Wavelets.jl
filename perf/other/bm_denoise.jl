
#include("setup_1d.jl")
using Wavelets

wf = waveletfilter(WT.db2)
wl = waveletls(WT.db2)
N = 8*1024;
x0 = rand(N);
L = nscales(N)  # int(log2(N)-2)
tn = 30

wt = wl
TI = false

println("denoise (N=",N,"), ", L, " levels")
f(x0, L, wt,TI) = for i = 1:tn; denoise(x0, wt=wt, TI=TI); end
f(x0, L, wt,TI);
@time f(x0, L, wt,TI);

#@profile f(x0,L,wt,TI);
#Profile.print()

