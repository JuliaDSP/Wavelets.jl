using Wavelets
using Base.Test

include("util.jl")

# ============= accuracy tests ================

include("data_fwt.jl")

# 1-d
wf = POfilter("db3")
L = nscales(n1)
yw = fwt(x1, L, wf)
xtt = iwt(yw,L,wf)

@test vecnorm(yw-y1)<4e-4
@test abs(vecnorm(x1)-vecnorm(yw))<1e-12
@test_approx_eq xtt x1

# 2-d
wf = POfilter("db3")
L = nscales(n2)
yw = fwt(x2, L, wf)
xtt = iwt(yw,L,wf)

@test vecnorm(yw-y2)<4e-4
@test abs(vecnorm(x2)-vecnorm(yw))<1e-12
@test_approx_eq xtt x2

# 2D lifting and filtering comparison
for WT in ("db1","db2")
	wf = POfilter(WT)
	wls = GPLS(WT)
	n = 64
	x = randn(n)
	
	x2 = copy(x)
	tmp = zeros(n)
	tmpsub = zeros(n)
	
	for L in (nscales(n),0,1,2)
		yf = fwt(x, L, wf)
		yls = fwt(x, L, wls)
		dwt!(x2, L, wls, true, tmp, oopc=true, oopv=tmpsub)
		
		@test_approx_eq yf yls
		@test_approx_eq yf x2
		
		ytf = iwt(yf, L, wf)
		ytls = iwt(yls, L, wls)
		dwt!(x2, L, wls, false, tmp, oopc=true, oopv=tmpsub)
		
		@test_approx_eq ytf x
		@test_approx_eq ytls x
		@test_approx_eq x2 x
		
	end
end

# TYPES AND SIZES

# 1-d
n = 8
wf = POfilter("db2")
L = 2
ft = Float64; x = ft[1:n]; y=fwt(x,wf);
@test Array{ft,1} == typeof(y) && length(y) == n
ft = Float32; x = ft[1:n]; y=fwt(x,wf);
@test Array{ft,1} == typeof(y) && length(y) == n
ft = Int64; x = ft[1:n]; y=fwt(x,wf);
@test Array{typeof(float(x[1])),1} == typeof(y) && length(y) == n
ft = Int32; x = ft[1:n]; y=fwt(x,wf);
@test Array{typeof(float(x[1])),1} == typeof(y) && length(y) == n

# 2-d
n = 8
wf = POfilter("db2")
L = 2
ft = Float64; x = ft[1:n]; y=fwt(x*x',wf);
@test Array{ft,2} == typeof(y) && length(y) == n*n
ft = Float32; x = ft[1:n]; y=fwt(x*x',wf);
@test Array{ft,2} == typeof(y) && length(y) == n*n
ft = Int64; x = ft[1:n]; y=fwt(int(x*x'),wf);
@test Array{typeof(float(x[1])),2} == typeof(y) && length(y) == n*n
ft = Int32; x = ft[1:n]; y=fwt(int(x*x'),wf);
@test Array{typeof(float(x[1])),2} == typeof(y) && length(y) == n*n




