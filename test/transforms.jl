
# ============= accuracy tests ================

# test against data made by data/make_filter_data.m in Octave
wname = {"Daubechies","Coiflet","Haar","Symmlet","Battle","Vaidyanathan","Beylkin"};
wnum = {[4:2:20],[2:5],[0],[4:10],[1,3,5],[0],[0]};
wvm = {[2:1:10],[4:2:10],[0],[4:10],[2,4,6],[0],[0]};

name = "filter"
data = vec(readdlm(joinpath(dirname(@__FILE__), "data", "filter1d_data.txt"),'\t'))
data2 = readdlm(joinpath(dirname(@__FILE__), "data", "filter2d_data.txt"),'\t')

for i = 1:length(wname)
    for num = 1:length(wnum[i])
        
        ye = vec(readdlm(joinpath(dirname(@__FILE__), "data", string(name,"1d_",wname[i],wnum[i][num],".txt")),'\t'))
        ye2 = readdlm(joinpath(dirname(@__FILE__), "data", string(name,"2d_",wname[i],wnum[i][num],".txt")),'\t')
        
        wn = wnum[i][num]
        wn!=0 && (wt = POfilter(wname[i],wvm[i][num]))
        wn==0 && (wt = POfilter(lowercase(wname[i][1:4])))
        y = fwt(data, wt)
        y2 = fwt(data2, wt)
        #println(wname[i],wn,' ',wvm[i][num], ' ',vecnorm(y - ye))
        #println(vecnorm(iwt(y,wt) - data))
        #println(wname[i],wn,' ',wvm[i][num], ' ',vecnorm(y2 - ye2))
        #println(vecnorm(iwt(y2,wt) - data2))
        #println(abs(vecnorm(data)-vecnorm(y)))
        #println(abs(vecnorm(data2)-vecnorm(y2)))
        @test vecnorm(y - ye) < 5e-9
        @test vecnorm(y2 - ye2) < 5e-9
        # a few bad cases have quite large norms
        if wname[i] != "Battle" && (wname[i] != "Coiflet" && wvm[i][num]==10)
            @test abs(vecnorm(data)-vecnorm(y)) < 1e-9
            @test abs(vecnorm(data2)-vecnorm(y2)) < 1e-9
            @test vecnorm(iwt(y,wt) - data) < 5e-7
            @test vecnorm(iwt(y,wt) - data) < 5e-7
        end
    end
end

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

# "out of place" for LS
wt = GPLS("db2")
x = randn(16)
@test_approx_eq fwt(x,2,wt) dwt!(similar(x),x,2,wt,true)


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


