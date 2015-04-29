
using Wavelets

wf = wavelet(WT.db4, WT.Filter)
wl = wavelet(WT.cdf97, WT.Lifting)
N = 1024;
x0 = rand(N,N);
L = nscales(N)  # int(log2(N)-2)
tn = 10

