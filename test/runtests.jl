using Wavelets
using Base.Test


# ============= accuracy tests ================

include("data_fwt.jl")

# 1D
wf = POfilter("db6")
L = nscales(n1)
yw = fwt(x1, L, wf)
xtt = iwt(yw,L,wf)

@assert vecnorm(yw-y1)<4e-4
@assert abs(vecnorm(x1)-vecnorm(yw))<1e-12
@test_approx_eq xtt x1

# 2D
wf = POfilter("db6")
L = nscales(n2)
yw = fwt(x2, L, wf)
xtt = iwt(yw,L,wf)

@assert vecnorm(yw-y2)<4e-4
@assert abs(vecnorm(x2)-vecnorm(yw))<1e-12
@test_approx_eq xtt x2


# ============= type tests ================

n = 8
wf = POfilter("db4")
L = 2
ft = Float64; x = ft[1:n]; y=fwt(x,wf);
@test Array{ft,1} == typeof(y) && length(y) == n
ft = Float32; x = ft[1:n]; y=fwt(x,wf);
@test Array{ft,1} == typeof(y) && length(y) == n
ft = Int64; x = ft[1:n]; y=fwt(x,wf);
@test Array{typeof(1.0),1} == typeof(y) && length(y) == n

n = 8
wf = POfilter("db4")
L = 2
ft = Float64; x = ft[1:n]; y=fwt(x*x',wf);
@test Array{ft,2} == typeof(y) && length(y) == n*n
ft = Float32; x = ft[1:n]; y=fwt(x*x',wf);
@test Array{ft,2} == typeof(y) && length(y) == n*n
ft = Int64; x = ft[1:n]; y=fwt(int(x*x'),wf);
@test Array{typeof(1.0),2} == typeof(y) && length(y) == n*n




