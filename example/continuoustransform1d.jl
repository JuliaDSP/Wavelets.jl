using Plots
using Wavelets
WT.paul10
supertype(typeof(WT.paul10))
wavelet(WT.morl)
wavelet(WT.morl,4)
plotlyjs()
J=11; n = 2^J
x = testfunction(n,"Bumps")
supertype(typeof(WT.dog0))
supertype(typeof(WT.paul5))
supertype(typeof(WT.morl))
EmphasizeFrequencyInfo = WT.Morlet(20)
y = cwt(x, wavelet(WT.morl))
heatmap(abs.(y)); plot!(x+90,label="")
y = cwt(x, wavelet(WT.dog1))
heatmap(abs.(y)); plot!(x+90,label="")
y = cwt(x, wavelet(WT.paul4))
heatmap(abs.(y)); plot!(x+90,label="")
