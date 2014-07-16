
using Wavelets

wf = POfilter("db6")
N = 1024*32;
x0 = rand(N);
L = int(log2(N)-2)
tn = 1000
fw=true

ft=Float64
si = Array(ft, wf.n-1)       # tmp filter vector
if fw
    dcfilter = mirror(convert(Vector{ft},wf.qmf))  #mqmf
    scfilter = reverse(convert(Vector{ft},wf.qmf))  #rqmf
else
    scfilter = convert(Vector{ft},wf.qmf)  #qmf
    dcfilter = reverse(mirror(convert(Vector{ft},wf.qmf)))  #mrqmf
end
snew = Array(ft, length(x0)>>1)


y=x0*0.0

println("dwt! (N=",N,"), ", L, " levels")
#f(x0, L, wf) = for i = 1:tn; fwt(x0, L, wf); end
f(x0, L, wf) = for i = 1:tn; dwt!(y, x0, L, wf, fw, dcfilter, scfilter, si, snew); end
f(x0, L, wf);
@time f(x0, L, wf);



