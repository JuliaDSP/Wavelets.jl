
using Wavelets

wf = waveletfilter(WT.db6)
N = 1024*32;
x0 = rand(N);
L = int(log2(N)-2)
tn = 1000

println("iwt (N=",N,"), ", L, " levels")
f(x0, L, wf) = for i = 1:tn; idwt(x0, wf, L); end
f(x0, L, wf);
@time f(x0, L, wf);



