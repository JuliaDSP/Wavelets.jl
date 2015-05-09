
print("plot: execute functions ...\n")

J = 8
n = 2^J
x = testfunction(n, "Bumps")
y = dwt(x, wavelet(WT.cdf97, WT.Lifting))
d,l = wplotdots(y, 0.1, n)
A = wplotim(y)


x = randn(64,64)
L = 2
xts = wplotim(x, L, wavelet(WT.db3, WT.Filter))

