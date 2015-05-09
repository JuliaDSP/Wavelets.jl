
using Wavelets
using PyPlot

J = 11
n = 2^J
x = testfunction(n,"Bumps")
y = dwt(x, wavelet(WT.cdf97, WT.Filter))
d,l = wplotdots(y, 0.1, n)
A = wplotim(y)

f, ax = plt.subplots(3, 1, sharex=true)
ax[1][:plot](x, "k")
ax[1][:set_xlim]([0,n])
ax[1][:set_ylabel](L"signal $x$")
ax[2][:scatter](d,l)
ax[2][:margins](0.05)
ax[2][:set_ylim](ax[2][:get_ylim]()[2:-1:1])
ax[2][:set_ylabel](L"level $j$")
ax[3][:imshow](A,aspect="auto",interpolation="none")
ax[3][:set_ylim]([J-0.5,-0.5])
ax[3][:set_ylabel](L"level $j$")

f[:set_size_inches](6,5)

savefig("transform1d_bumps.png")


