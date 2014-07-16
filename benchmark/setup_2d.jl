
using Wavelets

wf = POfilter("db4")
wl = GPLS("cdf9/7")
N = 1024;
x0 = rand(N,N);
L = nscales(N)  # int(log2(N)-2)
tn = 10

