
# ============= accuracy tests ================

# test against data made by data/make_filter_data.m in Octave
wname = {"Daubechies","Coiflet","Haar","Symmlet","Battle","Vaidyanathan","Beylkin"};
wnum = {[4:2:20],[2:5],[0],[4:10],[1,3,5],[0],[0]};
wvm = {[2:1:10],[4:2:10],[0],[4:10],[2,4,6],[0],[0]};

name = "filter"
data = vec(readdlm(joinpath(dirname(@__FILE__), "data", "filter1d_data.txt"),'\t'))
data2 = readdlm(joinpath(dirname(@__FILE__), "data", "filter2d_data.txt"),'\t')

stderr = 1e-9*sqrt(length(data))
stderr2 = 1e-9*sqrt(length(data2))

for i = 1:length(wname)
    for num = 1:length(wnum[i])
        # expected transforms
        ye = vec(readdlm(joinpath(dirname(@__FILE__), "data", string(name,"1d_",wname[i],wnum[i][num],".txt")),'\t'))
        ye2 = readdlm(joinpath(dirname(@__FILE__), "data", string(name,"2d_",wname[i],wnum[i][num],".txt")),'\t')
        
        wn = wnum[i][num]
        wn!=0 && (wt = POfilter(wname[i],wvm[i][num]))
        wn==0 && (wt = POfilter(lowercase(wname[i][1:4])))
        # transform data
        y = fwt(data, wt)
        y2 = fwt(data2, wt)
        
        @test_vecnorm_eq_eps y ye stderr
        @test_vecnorm_eq_eps y2 ye2 stderr2
        # a few bad cases have quite large norms
        if wname[i] != "Battle" && (wname[i] != "Coiflet" && wvm[i][num]==10)
            @test abs(vecnorm(data)-vecnorm(y)) < 1e-9
            @test abs(vecnorm(data2)-vecnorm(y2)) < 1e-9
            @test_vecnorm_eq_eps iwt(y,wt) data stderr*100
            @test_vecnorm_eq_eps iwt(y2,wt) data2 stderr2*100
        end
    end
end

# 1-d and 2-d lifting and filtering comparison
for WT in ("db1","db2")
    wf = POfilter(WT)
    wls = GPLS(WT)
    n = 64
    x = randn(n)
    stderr = 1e-10*sqrt(n)
    stderr2 = 1e-10*sqrt(n*n)
    
    x2 = copy(x)
    tmp = zeros(n)
    tmpsub = zeros(n)
    
    for L in (nscales(n),0,1,2)
        yf = fwt(x, L, wf)
        yls = fwt(x, L, wls)
        dwt!(x2, L, wls, true, tmp, oopc=true, oopv=tmpsub)
        
        @test_vecnorm_eq_eps yf yls stderr
        @test_vecnorm_eq_eps yf x2 stderr
        
        ytf = iwt(yf, L, wf)
        ytls = iwt(yls, L, wls)
        dwt!(x2, L, wls, false, tmp, oopc=true, oopv=tmpsub)
        
        @test_vecnorm_eq_eps ytf x stderr
        @test_vecnorm_eq_eps ytls x stderr
        @test_vecnorm_eq_eps x2 x stderr
        
    end
    
    x = randn(n,n)
    x2 = copy(x)
    for L in (nscales(n),0,1,2)
        yf = fwt(x, L, wf)
        yls = fwt(x, L, wls)
        
        @test_vecnorm_eq_eps yf yls stderr2
        @test_approx_eq yf yls
        
        ytf = iwt(yf, L, wf)
        ytls = iwt(yls, L, wls)
        
        @test_vecnorm_eq_eps ytf x stderr2
        @test_vecnorm_eq_eps ytls x stderr2
        
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

function makefwt(ft::Type, n, wf, L)
	x0 = rand(-5:5, n)
	x = zeros(ft, n)
	copy!(x,x0)
	return (x, fwt(x, wf, L))
end

# 1-d
n = 8; wf = POfilter("db2"); L = 2
sett = (n,wf,L)

ft = Float64; x, y = makefwt(ft, sett...)
@test Array{ft,1} == typeof(y) && length(y) == n
ft = Float32; x, y = makefwt(ft, sett...)
@test Array{ft,1} == typeof(y) && length(y) == n
ft = Int64; x, y = makefwt(ft, sett...)
@test Array{typeof(float(x[1])),1} == typeof(y) && length(y) == n
ft = Int32; x, y = makefwt(ft, sett...)
@test Array{typeof(float(x[1])),1} == typeof(y) && length(y) == n
@test_approx_eq fwt(x,wf) fwt(float(x),wf)
@test_approx_eq fwt(x,wf,2) fwt(x,2,wf)


# 2-d
sett = ((n,n),wf,L)
ft = Float64; x, y = makefwt(ft, sett...)
@test Array{ft,2} == typeof(y) && length(y) == n*n
ft = Float32; x, y = makefwt(ft, sett...)
@test Array{ft,2} == typeof(y) && length(y) == n*n
ft = Int64; x, y = makefwt(ft, sett...)
@test Array{typeof(float(x[1])),2} == typeof(y) && length(y) == n*n
ft = Int32; x, y = makefwt(ft, sett...)
@test Array{typeof(float(x[1])),2} == typeof(y) && length(y) == n*n
@test_approx_eq fwt(x,wf) fwt(float(x),wf)
@test_approx_eq fwt(x,wf,2) fwt(x,2,wf)


# ============= error tests ================

type wunknownt <: WaveletType end
uwt = wunknownt()
EE = ErrorException
@test_throws EE fwt(randn(4),uwt)
@test_throws EE fwt(randn(4,4),uwt)




