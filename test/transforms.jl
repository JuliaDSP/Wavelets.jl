

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

# 1-d and 2-d lifting and filtering comparison
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
    
    x = randn(n,n)
    x2 = copy(x)
    for L in (nscales(n),0,1,2)
        yf = fwt(x, L, wf)
        yls = fwt(x, L, wls)
        
        @test_approx_eq yf yls
        
        ytf = iwt(yf, L, wf)
        ytls = iwt(yls, L, wls)
        
        @test_approx_eq ytf x
        @test_approx_eq ytls x
        
    end
end

# ============= transform functionality ================

# column-wise 1-d
wf = POfilter("db2")
x = randn(16,2)
y = copy(x)
y[:,1] = fwt(vec(x[:,1]),wf)
y[:,2] = fwt(vec(x[:,2]),wf)
@test_approx_eq fwtc(x,wf) y

# column-wise 2-d
n = 16
x = randn(n,n,2)
y = copy(x)
y[:,:,1] = fwt(reshape(x[:,:,1],n,n),wf)
y[:,:,2] = fwt(reshape(x[:,:,2],n,n),wf)
@test_approx_eq fwtc(x,wf) y

# "inplace" for filter
wf = POfilter("db2")
x = randn(16)
@test_approx_eq fwt(x,2,wf) dwt!(copy(x),2,wf,true)


# ============= types and sizes ================

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

ft = Int64
x = ft[1:n]
@test_approx_eq fwt(x,wf,2) fwt(x,2,wf)
@test_approx_eq fwt(x,wf) fwt(float(x),wf)
ft = Float64
x = ft[1:n]
@test_approx_eq fwt(x,wf,2) fwt(x,2,wf)


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

ft = Int64
x = ft[1:n]
@test_approx_eq fwt(x*x',wf,2) fwt(x*x',2,wf)
@test_approx_eq fwt(x*x',wf) fwt(float(x*x'),wf)
ft = Float64
x = ft[1:n]
@test_approx_eq fwt(x*x',wf,2) fwt(x*x',2,wf)



# ============= error tests ================

type wunknownt <: WaveletType
    a::Int
end
uwt = wunknownt(1)

if VERSION < v"0.3.0-"  # v0.2 support
    tests=quote
    @test_throws fwt(randn(4),uwt)
    @test_throws fwt(randn(4,4),uwt)
    end
    eval(tests)
else
    EE = ErrorException
    tests=quote
    @test_throws EE fwt(randn(4),uwt)
    @test_throws EE fwt(randn(4,4),uwt)
    end
    eval(tests)
end


