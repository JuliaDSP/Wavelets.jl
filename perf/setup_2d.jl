
using Wavelets

wf = waveletfilter(WT.db4)
wl = waveletls(WT.cdf97)
N = 1024;
x0 = rand(N,N);
L = nscales(N)  # int(log2(N)-2)
tn = 10

