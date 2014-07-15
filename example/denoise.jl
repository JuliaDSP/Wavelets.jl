
using Wavelets
using PyPlot

n = 2^11;
x0 = testfunction(n,"Doppler")
x = x0 + 0.05*randn(n)
y = denoise(x,TI=true)
f, ax = plt.subplots(3,1, sharex=true)

cc=0.6
ax[1][:plot](x0,"b")
ax[1][:set_title]("original")
ax[2][:plot](x,"g")
ax[2][:set_title]("noisy")
ax[3][:plot](y,"r")
ax[3][:set_title]("denoised (TI VisuShrink)")
for i=1:length(ax)
    ax[i][:set_xlim]([0,n])
    ax[i][:set_ylim]([-cc,cc])
end

savefig("denoise_doppler.png")


