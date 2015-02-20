
# ============= accuracy tests ================
print("transforms: accuracy tests ...\n")

# test against data made by data/make_filter_data.m in Octave
wname = Any["Daubechies","Coiflet","Haar","Symmlet","Battle","Vaidyanathan","Beylkin"];
wtype = Any[WT.Daubechies, WT.Coiflet, WT.Haar, WT.Symlet, WT.Battle, WT.Vaidyanathan, WT.Beylkin];
wnum = Any[[4:2:20],[2:5],[0],[4:10],[1,3,5],[0],[0]];
wvm = Any[[2:1:10],[4:2:10],[0],[4:10],[2,4,6],[0],[0]];

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
        
        wtc = wtype[i]
        if wnum[i][num] != 0 
            class = wtc{wvm[i][num]}()
        else
            class = wtc()
        end
        wt = waveletfilter(class)

        # transform data
        y = dwt(data, wt)
        y2 = dwt(data2, wt)
        
        @test_vecnorm_eq_eps y ye stderr
        @test_vecnorm_eq_eps y2 ye2 stderr2
        # a few bad cases have quite large norms
        if wname[i] != "Battle" && (wname[i] != "Coiflet" && wvm[i][num]==10)
            @test abs(vecnorm(data)-vecnorm(y)) < 1e-9
            @test abs(vecnorm(data2)-vecnorm(y2)) < 1e-9
            @test_vecnorm_eq_eps idwt(y,wt) data stderr*100
            @test_vecnorm_eq_eps idwt(y2,wt) data2 stderr2*100
        end
    end
end

# 1-d and 2-d lifting and filtering comparison
for wclass in (WT.db1, WT.db2)
    wf = waveletfilter(wclass)
    wls = GLS(wclass)
    n = 64
    x = randn(n)
    stderr = 1e-10*sqrt(n)
    stderr2 = 1e-10*sqrt(n*n)
    
    x2 = copy(x)
    tmp = zeros(n)
    
    for L in (ndyadicscales(n),0,1,2)
        yf = dwt(x, wf, L)
        yls = dwt(x, wls, L)
        dwt!(x2, wls, L, true, tmp)
        
        @test_vecnorm_eq_eps yf yls stderr
        @test_vecnorm_eq_eps yf x2 stderr
        
        ytf = idwt(yf, wf, L)
        ytls = idwt(yls, wls, L)
        dwt!(x2, wls, L, false, tmp)
        
        @test_vecnorm_eq_eps ytf x stderr
        @test_vecnorm_eq_eps ytls x stderr
        @test_vecnorm_eq_eps x2 x stderr
        
    end
    
    x = randn(n,n)
    x2 = copy(x)
    for L in (ndyadicscales(n),0,1,2)
        yf = dwt(x, wf, L)
        yls = dwt(x, wls, L)
        
        @test_vecnorm_eq_eps yf yls stderr2
        @test_approx_eq yf yls
        
        ytf = idwt(yf, wf, L)
        ytls = idwt(yls, wls, L)
        
        @test_vecnorm_eq_eps ytf x stderr2
        @test_vecnorm_eq_eps ytls x stderr2
        
    end
end

# ============= transform functionality ================
print("transforms: transform functionality ...\n")

# column-wise 1-d
wf = waveletfilter(WT.db2)
x = randn(16,2)
y = copy(x)
y[:,1] = dwt(vec(x[:,1]),wf)
y[:,2] = dwt(vec(x[:,2]),wf)
@test_approx_eq dwtc(x,wf) y
@test_approx_eq dwtc(x',wf,ndyadicscales(16),1)' y

# column-wise 2-d
n = 16
x = randn(n,n,2)
y = copy(x)
y[:,:,1] = dwt(reshape(x[:,:,1],n,n),wf)
y[:,:,2] = dwt(reshape(x[:,:,2],n,n),wf)
@test_approx_eq dwtc(x,wf) y
x = randn(n,2,n)
y = copy(x)
y[:,1,:] = dwt(reshape(x[:,1,:],n,n),wf)
y[:,2,:] = dwt(reshape(x[:,2,:],n,n),wf)
@test_approx_eq dwtc(x,wf,ndyadicscales(16),2) y

# "inplace" for filter
wf = waveletfilter(WT.db2)
x = randn(16)
@test_approx_eq dwt(x,wf,2) dwt!(copy(x),wf,2,true)

# "out of place" for LS
wt = GLS(WT.db2)
x = randn(16)
@test_approx_eq dwt(x,wt,2) dwt!(similar(x),x,wt,2,true)


# ============= types and sizes ================
print("transforms: types and sizes ...\n")

function makedwt(ft::Type, n, wf, L)
	x0 = rand(-5:5, n)
	x = zeros(ft, n)
	copy!(x,x0)
	return (x, dwt(x, wf, L))
end

# 1-d
n = 8; wf = waveletfilter(WT.db2); L = 2
sett = (n,wf,L)

ft = Float64; x, y = makedwt(ft, sett...)
@test Array{ft,1} == typeof(y) && length(y) == n
ft = Float32; x, y = makedwt(ft, sett...)
@test Array{ft,1} == typeof(y) && length(y) == n
ft = Int64; x, y = makedwt(ft, sett...)
@test Array{typeof(float(x[1])),1} == typeof(y) && length(y) == n
ft = Int32; x, y = makedwt(ft, sett...)
@test Array{typeof(float(x[1])),1} == typeof(y) && length(y) == n
@test_approx_eq dwt(x,wf) dwt(float(x),wf)


# 2-d
sett = ((n,n),wf,L)
ft = Float64; x, y = makedwt(ft, sett...)
@test Array{ft,2} == typeof(y) && length(y) == n*n
ft = Float32; x, y = makedwt(ft, sett...)
@test Array{ft,2} == typeof(y) && length(y) == n*n
ft = Int64; x, y = makedwt(ft, sett...)
@test Array{typeof(float(x[1])),2} == typeof(y) && length(y) == n*n
ft = Int32; x, y = makedwt(ft, sett...)
@test Array{typeof(float(x[1])),2} == typeof(y) && length(y) == n*n
@test_approx_eq dwt(x,wf) dwt(float(x),wf)

# non-Array type
wt = waveletls(WT.db2)
x = randn(16, 16)
xs = sub(copy(x), 1:16, 1:16)
@test_approx_eq dwt(x,wt,2) dwt(xs,wt,2)

#util functions
for class in (WT.haar, WT.db2, WT.cdf97)
    WT.class(class)
    WT.name(class)
    WT.vanishingmoments(class)
end
class = WT.db1
wt = waveletfilter(class)
@test length(wt) == 2
@test_approx_eq wt.qmf*0.7 scale(wt, 0.7).qmf
@test OrthoFilter{PerBoundary}(wt.qmf, "db1").qmf == waveletfilter(class).qmf

# ============= error tests ================
print("transforms: error tests ...\n")

type wunknownt <: DiscreteWavelet end
uwt = wunknownt()
EE = Exception
@test_throws EE dwt(randn(4),uwt)
@test_throws EE dwt(randn(4,4),uwt)
@test_throws EE waveletfilter(WT.Coiflet{33}())
@test_throws EE waveletfilter("db2asdsad")
@test_throws EE waveletfilter("db2", "ppppp")

# ============= WPT ================
print("transforms: WPT ...\n")

wf = waveletfilter(WT.db2)
x = randn(16)

L = 1
wp = wpt(x,wf,L)
dw = dwt(x,wf,L)
@test_approx_eq wp dw
@test_approx_eq iwpt(wp,wf,L) x

L = 2
wp = wpt(x,wf,L)
dw = dwt(x,wf,L)
dw2 = copy(dw)
dw2[9:end] = dwt(dw[9:end],wf,1)
@test_approx_eq dw[1:8] wp[1:8]
@test_approx_eq dw2 wp
@test_approx_eq iwpt(wp,wf,L) x

L = 3
wp = wpt(x,wf,L)
dw = dwt(x,wf,L)
@test_approx_eq dw[1:4] wp[1:4]
@test_approx_eq dwt(dw2[5:8],wf,1) wp[5:8]
@test_approx_eq dwt(dw2[9:12],wf,1) wp[9:12]
@test_approx_eq dwt(dw2[13:16],wf,1) wp[13:16]
@test_approx_eq iwpt(wp,wf,L) x

wl = waveletls(WT.db2)
x = randn(128)
@test_approx_eq iwpt(wpt(x,wf),wf) x
@test_approx_eq iwpt(wpt(x,wl),wl) x
@test_approx_eq wpt(x,wl,1) wpt(x,wf,1)
@test_approx_eq wpt(x,wl,2) wpt(x,wf,2)
@test_approx_eq wpt(x,wl,4) wpt(x,wf,4)
@test_approx_eq wpt(x,wl) wpt(x,wf)

wl = waveletls(WT.db2)
n = 128
x = randn(n)
for L = 0:ndyadicscales(n)
    @test_approx_eq wpt(x, wl, maketree(n, L, :dwt)) dwt(x, wl, L)
    @test_approx_eq iwpt(x, wl, maketree(n, L, :dwt)) idwt(x, wl, L)
    @test_approx_eq wpt(x, wf, maketree(n, L, :dwt)) dwt(x, wf, L)
    @test_approx_eq iwpt(x, wf, maketree(n, L, :dwt)) idwt(x, wf, L)
end

# ============= tranform low level functions ================

#...



