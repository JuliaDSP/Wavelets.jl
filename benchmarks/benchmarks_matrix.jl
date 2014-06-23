
using Wavelets

wf = POfilter("db6")
N = 1024*32;
N2 = 50
x0 = rand(N);
L = int(log2(N)-2)
tn = 10
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

A=zeros(N,N2)
y=x0*0.0


println("fwt (N=",N,"), ", L, " levels")
#f(x0, L, wf) = for i = 1:tn; fwt(x0, L, wf); end
f(x0, L, wf) = for i = 1:tn; 
for j=1:N2
	A[:,j] = fwt(x0, L, wf);
end	
end
f(x0, L, wf);
@time f(x0, L, wf);

println("dwt! (N=",N,"), ", L, " levels")
#f(x0, L, wf) = for i = 1:tn; fwt(x0, L, wf); end
f(x0, L, wf) = for i = 1:tn; 
for j=1:N2
	y = sub(A,(j-1)*N+1:j*N)
	dwt!(y, x0, L, wf, fw, dcfilter, scfilter, si, snew); 
end	
end
f(x0, L, wf);
@time f(x0, L, wf);





