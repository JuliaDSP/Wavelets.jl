
using Wavelets

wf = POfilter("db4")
N = 1024;
x0 = rand(N,N);
L = int(log2(N)-2)
tn = 10
#y=x0*0.0

println("iwt (N=",N,"x",N,"), ", L, " levels")
f(x0, L, wf) = for i = 1:tn; iwt(x0, L, wf); end
f(x0, L, wf);
@time f(x0, L, wf);



