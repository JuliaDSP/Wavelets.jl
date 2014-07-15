
using Wavelets
using PyPlot

n = 2^11
x0 = testfunction(n,"Bumps")
y = fwt(x0,GPLS("cdf9/7"))
f, ax = plt.subplots(2,1, sharex=true)

cc=0.6
ax[1][:plot](x0,"b")
ax[1][:set_title]("signal",fontsize=11)
ax[2][:plot](y,"k")
ax[2][:set_title]("CDF 9/7 transform",fontsize=11)
for i=1:length(ax)
    ax[i][:set_xlim]([0,n])
end
f[:set_size_inches](6,3)

savefig("transform1d_bumps.png")

