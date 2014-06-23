
using Wavelets

wf = POfilter("db2")
N = 1024*32;
#N=16
x0 = rand(N);
#x0=zeros(N)
#x0[1:2]=[1,-3.3]
L = int(log2(N)-2)
#L = 1
tn = 100
fw=true

y=x0*0.0

println("dwt! (N=",N,"), ", L, " levels")
#f(x0, L, wf) = for i = 1:tn; fwt(x0, L, wf); end
f(y, x0, L, wf, fw) = for i = 1:tn; dwt!(y, x0, L, wf, fw); end
f(y, x0, L, wf, fw);
@time f(y, x0, L, wf, fw);

println(y[1:16])
#println(y[N/2+1:N/2+3])
println(norm(y))

println("dwtd4! (N=",N,"), ", L, " levels")
#f(x0, L, wf) = for i = 1:tn; fwt(x0, L, wf); end
f(y, x0, L, fw) = for i = 1:tn; dwtd4!(y, x0, L, fw); end
f(y, x0, L, fw);
@time f(y, x0, L, fw);

println(y[1:16])
#println(y[N/2+1:N/2+3])
println(norm(y))
println(norm(x0))
copy!(y,x0)

tmp = Array(eltype(x0),N>>2)
println("dwtliftd4! (N=",N,"), ", L, " levels")
#f(x0, L, wf) = for i = 1:tn; fwt(x0, L, wf); end
f(y, x0, L, fw,tmp) = for i = 1:tn; dwtliftd4!(y, y, L, fw, tmp); end
f(y, x0, L, fw,tmp);
@time f(y, x0, L, fw,tmp);
copy!(y,x0)
dwtliftd4!(y, y, L, fw, tmp);

println(y[1:16])
println(norm(y))
#println(y[N/2+1:N/2+3])

println("dwtd4! (N=",N,"), ", L, " levels")
#f(x0, L, wf) = for i = 1:tn; fwt(x0, L, wf); end
f(y, x0, L, fw) = for i = 1:tn; dwtd4!(y, x0, L, fw); end
f(y, x0, L, fw);
@time f(y, x0, L, fw);

