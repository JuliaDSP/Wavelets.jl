
print("threshold: execute functions ...\n")

#threshold
x = randn(200)*2
for th in (BiggestTH(), HardTH(), SoftTH(), SemiSoftTH(), SteinTH())
    threshold(x, th, 2)
end
for th in (PosTH(), NegTH())
    threshold(x, th)
end

# denoising
vs = VisuShrink(10)
n = 2^8
x0 = testfunction(n, "Doppler")
x = x0 + 0.05*randn(n)
y = denoise(x, TI=true)
y = denoise(x, TI=true, nspin = 8)
y = denoise(x, TI=false)
y = denoise(x, nothing)
y = denoise(randn(32,32), TI=true)

# best basis
wt = wavelet(WT.db4)
x = sin.(4*linspace(0,2*pi-eps(),1024))
tree = bestbasistree(x, wt)
xtb = wpt(x, wt, tree)
@test iwpt(xtb, wt, tree) ≈ x

x = sin.(4*linspace(0,2*pi-eps(),5*64)) # non-dyadic
tree = bestbasistree(x, wt)
xtb = wpt(x, wt, tree)
@test iwpt(xtb, wt, tree) ≈ x

#matching pursuit
N = 128
M = 64
y=randn(N)
Arr = randn(M,N)
func(a::AbstractVector) = Arr*a
funct(a::AbstractVector) = Arr'*a
x = func(y)
matchingpursuit(x, func, funct, 0.1)
