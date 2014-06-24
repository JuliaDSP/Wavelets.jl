
using Wavelets

wf = POfilter("db2")
N = 1024*32;
#N=16
x0 = rand(N);
#x0=zeros(N)
#x0[1:2]=[1,-3.3]
L = int(log2(N)-2)
L = 1
tn = 1000
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

# LIFTING

copy!(y,x0)
tmp = Array(eltype(x0),N>>2)
ls = GPLS([ LSstep('p',[-1.7320508075688772],0),
            LSstep('u',[-0.0669872981077807,0.4330127018922193],1),
            LSstep('p',[1.0],-1)],
            0.5176380902050414,1.9318516525781364,"db2")

println("dwt! (N=",N,"), ", L, " levels (lifting)")
#f(x0, L, wf) = for i = 1:tn; fwt(x0, L, wf); end
f(y, L, ls, fw, tmp) = for i = 1:tn; dwt!(y, L, ls, fw, tmp); end
f(y, L, ls, fw, tmp);
@time f(y, L, ls, fw, tmp);

copy!(y,x0)
dwt!(y, L, ls, fw, tmp);

println(y[1:16])
println(norm(y))
#println(y[N/2+1:N/2+3])

@profile f(y, L, ls, fw, tmp);
Profile.print()

