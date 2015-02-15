
using Wavelets

wf = waveletfilter(WT.db2)
wl = waveletls(WT.db2)
N = 1024*1024;
x0 = rand(N);
L = nscales(N)  # int(log2(N)-2)
tn = 10

