
print("threshold: execute functions ...\n")

# denoising
n = 2^8
x0 = testfunction(n, "Doppler")
x = x0 + 0.05*randn(n)
y = denoise(x, TI=true)
y = denoise(x, TI=false)
y = denoise(x, nothing)

# best basis
wt = wavelet(WT.db4)
x = sin(4*linspace(0,2*pi-eps(),1024))
tree = bestbasistree(x, wt)
xtb = wpt(x, wt, tree)
xt = dwt(x, wt)

