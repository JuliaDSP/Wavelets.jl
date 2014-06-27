
using Wavelets

wf = GPLS("cdf9/7")
N = 1024;
x0 = rand(N,N);
L = int(log2(N))
tn = 10
#y=x0*0.0

println("fwt by lifting (N=",N,"x",N,"), ", L, " levels")
#f(x0, L, wf) = for i = 1:tn; fwt(x0, L, wf); end
f(x0, L, wf) = for i = 1:tn; fwt(x0, L, wf); end
f(x0, L, wf);
@time f(x0, L, wf);

#@profile for i=1:tn dwt!(x0, L, wf, true); end
#@profile for i=1:tn fwt(x0, L, wf); end
#Profile.print()

