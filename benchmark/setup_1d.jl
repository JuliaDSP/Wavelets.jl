
using Wavelets

wf = POfilter("db2")
wl = GPLS("db2")
N = 1024*1024;
x0 = rand(N);
L = nscales(N)  # int(log2(N)-2)
tn = 10

