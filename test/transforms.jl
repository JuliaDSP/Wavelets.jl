
import Compat.view

# ============= accuracy tests ================
print("transforms: accuracy tests ...\n")

# test against data made by data/make_filter_data.m in Octave
wname = Any["Daubechies","Coiflet","Haar","Symmlet","Battle","Vaidyanathan","Beylkin"];
wtype = Any[WT.Daubechies, WT.Coiflet, WT.Haar, WT.Symlet, WT.Battle, WT.Vaidyanathan, WT.Beylkin];
wnum = Any[collect(4:2:20),collect(2:5),[0],collect(4:10),[1,3,5],[0],[0]];
wvm = Any[collect(2:1:10),collect(4:2:10),[0],collect(4:10),[2,4,6],[0],[0]];

name = "filter"
data = vec(readdlm(joinpath(dirname(@__FILE__), "data", "filter1d_data.txt"),'\t'))
data2 = readdlm(joinpath(dirname(@__FILE__), "data", "filter2d_data.txt"),'\t')

stderr = 1e-9*sqrt(length(data))
stderr2 = 1e-9*sqrt(length(data2))

for i in eachindex(wname)
    for num in eachindex(wnum[i])
        # expected transforms
        ye = vec(readdlm(joinpath(dirname(@__FILE__), "data", string(name,"1d_",wname[i],wnum[i][num],".txt")),'\t'))
        ye2 = readdlm(joinpath(dirname(@__FILE__), "data", string(name,"2d_",wname[i],wnum[i][num],".txt")),'\t')

        wtc = wtype[i]
        if wnum[i][num] != 0
            class = wtc{wvm[i][num]}()
        else
            class = wtc()
        end
        wt = wavelet(class, WT.Filter)

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

# 1-D, 2-D, 3-D lifting vs filtering / inverse vs original tests
for wclass in (WT.db1, WT.db2)
    wf = wavelet(wclass, WT.Filter)
    wls = wavelet(wclass, WT.Lifting)
    n = 64
    x = randn(n)
    stderr = 1e-10*sqrt(n)
    stderr2 = 1e-10*sqrt(n*n)

    x2 = copy(x)

    for L in (ndyadicscales(n),0,1,2)
        yf = dwt(x, wf, L)
        yls = dwt(x, wls, L)

        # filter vs lifting
        @test_vecnorm_eq_eps yf yls stderr

        ytf = idwt(yf, wf, L)
        ytls = idwt(yls, wls, L)

        # inverse vs original
        @test_vecnorm_eq_eps ytf x stderr
        @test_vecnorm_eq_eps ytls x stderr
    end

    x = randn(n,n)
    x2 = copy(x)
    for L in (ndyadicscales(n),0,1,2)
        yf = dwt(x, wf, L)
        yls = dwt(x, wls, L)

        # filter vs lifting
        @test_vecnorm_eq_eps yf yls stderr2

        ytf = idwt(yf, wf, L)
        ytls = idwt(yls, wls, L)

        # inverse vs original
        @test_vecnorm_eq_eps ytf x stderr2
        @test_vecnorm_eq_eps ytls x stderr2
    end

    x = randn(n,n,n)
    x2 = copy(x)
    for L in (ndyadicscales(n),0,1,2)
        yf = dwt(x, wf, L)
        yls = dwt(x, wls, L)

        # filter vs lifting
        @test_vecnorm_eq_eps yf yls stderr2

        ytf = idwt(yf, wf, L)
        ytls = idwt(yls, wls, L)

        # inverse vs original
        @test_vecnorm_eq_eps ytf x stderr2
        @test_vecnorm_eq_eps ytls x stderr2

    end
end

# ============= transform functionality ================
print("transforms: transform functionality ...\n")

#=
# "inplace" for filter
wf = wavelet(WT.db2, WT.Filter)
x = randn(16)
@test_approx_eq dwt(x,wf,2) dwt!(copy(x),wf,2)

# "out of place" for LS
wt = GLS(WT.db2)
x = randn(16)
@test_approx_eq dwt(x,wt,2) dwt!(similar(x),x,wt,2)
=#

# ============= types and sizes ================
print("transforms: types and sizes ...\n")

function makedwt(ft::Type, n, wf, L)
	x0 = rand(-5:5, n)
	x = zeros(ft, n)
	copy!(x,x0)
	return (x, dwt(x, wf, L))
end

# 1-d
n = 8; wf = wavelet(WT.db2); L = 2
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

# Complex types
for wfc in (wavelet(WT.db2), wavelet(WT.db2, WT.Lifting))
    xc = zeros(Complex128, n)
    map!(i->rand(Complex128), xc)
    yc = dwt(xc, wfc, L)
    @test Array{Complex128,1} == typeof(yc) && length(yc) == n
    @test_approx_eq xc idwt(dwt(xc, wfc), wfc)
end

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
wt = wavelet(WT.db2, WT.Lifting)
x = randn(16, 16)
xs = view(copy(x), 1:16, 1:16)
@test_approx_eq dwt(x,wt,2) dwt(xs,wt,2)

#util functions
for class in (WT.haar, WT.db2, WT.cdf97)
    WT.class(class)
    WT.name(class)
    WT.vanishingmoments(class)
end
class = WT.db1
wt = wavelet(class, WT.Filter)
@test length(wt) == 2
@test_approx_eq wt.qmf*0.7 WT.scale(wt, 0.7).qmf

# ============= error tests ================
print("transforms: error tests ...\n")

type wunknownt <: DiscreteWavelet end
uwt = wunknownt()
EE = Exception
@test_throws EE dwt(randn(4),uwt)
@test_throws EE dwt(randn(4,4),uwt)
@test_throws EE wavelet(WT.Coiflet{33}())
@test_throws EE wavelet("db2asdsad")
@test_throws EE wavelet("db2", "ppppp")

# ============= WPT ================
print("transforms: WPT ...\n")

wf = wavelet(WT.db2, WT.Filter)
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

wl = wavelet(WT.db2, WT.Lifting)
x = randn(128)
@test_approx_eq iwpt(wpt(x,wf),wf) x
@test_approx_eq iwpt(wpt(x,wl),wl) x
@test_approx_eq wpt(x,wl,1) wpt(x,wf,1)
@test_approx_eq wpt(x,wl,2) wpt(x,wf,2)
@test_approx_eq wpt(x,wl,4) wpt(x,wf,4)
@test_approx_eq wpt(x,wl) wpt(x,wf)

wl = wavelet(WT.db2, WT.Lifting)
n = 128
x = randn(n)
for L = 0:ndyadicscales(n)
    @test_approx_eq wpt(x, wl, maketree(n, L, :dwt)) dwt(x, wl, L)
    @test_approx_eq iwpt(x, wl, maketree(n, L, :dwt)) idwt(x, wl, L)
    @test_approx_eq wpt(x, wf, maketree(n, L, :dwt)) dwt(x, wf, L)
    @test_approx_eq iwpt(x, wf, maketree(n, L, :dwt)) idwt(x, wf, L)
end

# non-dyadic tests
wt = wavelet(WT.db2, WT.Filter)
n = 5*8
x = randn(n)
for L in 0:maxtransformlevels(n)
    for wt in (wavelet(WT.db2, WT.Filter), wavelet(WT.db2, WT.Lifting))
        @test_approx_eq wpt(x, wt, maketree(n, L, :dwt)) dwt(x, wt, L)
        @test_approx_eq iwpt(x, wt, maketree(n, L, :dwt)) idwt(x, wt, L)
    end
end


# ============= tranform low level functions ================

#...


# ============= wavelet types ================
print("transforms: wavelet types ...\n")

wt = wavelet(WT.db2, WT.Filter)
@test length(wt) == 4
same = WT.name(wt) == "db2"
@test same

wt = wavelet(WT.db2, WT.Lifting)
@test WT.name(wt) == "db2"
