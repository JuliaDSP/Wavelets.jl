
using Wavelets

wf = POfilter("db12")
N = 1024*32;
x0 = rand(N);
L = int(log2(N)-2)
tn = 1000

y=x0*0.0

println("fwt (N=",N,"), ", L, " levels")
f(x0, L, wf) = for i = 1:tn; fwt(x0, L, wf); end
#f(x0, L, wf) = for i = 1:tn; dwt!(y,x0, L, wf,true); end
f(x0, L, wf);
@time f(x0, L, wf);



