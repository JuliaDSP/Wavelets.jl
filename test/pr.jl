using Wavelets


wf = POfilter("db2")
N = 16;
x0 = [1:N]*1.0;
x0[3]=-2.2
L = int(log2(N)-2)
L = 2

println(round(fwt(x0,L,wf),3))
wf = GPLS("cdf9/7")
println(round(fwt(x0,L,wf),3))



