
using Wavelets

wf = wavelet(WT.db2, WT.Filter)
wl = wavelet(WT.db2, WT.Lifting)
N = 1024*1024;
x0 = rand(N);
L = maxtransformlevels(N)  # int(log2(N)-2)
tn = 10

