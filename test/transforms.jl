
@testset "Accuracy" begin
    # test against data made by data/make_filter_data.m in Octave
    wname = Any["Daubechies","Coiflet","Haar","Symmlet","Battle","Vaidyanathan","Beylkin"];
    wtype = Any[WT.Daubechies, WT.Coiflet, WT.Haar, WT.Symlet, WT.Battle, WT.Vaidyanathan, WT.Beylkin];
    wnum = Any[collect(4:2:20),collect(2:5),[0],collect(4:10),[1,3,5],[0],[0]];
    wvm = Any[collect(2:1:10),collect(4:2:10),[0],collect(4:10),[2,4,6],[0],[0]];

    #print("transforms: reading 1d and 2d data ...\n")
    name = "filter"
    data = vec(readdlm(joinpath(dirname(@__FILE__), "data", "filter1d_data.txt"),'\t'))
    data2 = readdlm(joinpath(dirname(@__FILE__), "data", "filter2d_data.txt"),'\t')

    #print("transforms: testing 1d and 2d data ...\n")
    stderr1 = 1e-9*sqrt(length(data))
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

            @test (@vecnorm_eq_eps y ye stderr1)
            @test (@vecnorm_eq_eps y2 ye2 stderr2)
            # a few bad cases have quite large norms
            if wname[i] != "Battle" && (wname[i] != "Coiflet" && wvm[i][num]==10)
                @test abs(norm(data)-norm(y)) < 1e-9
                @test abs(norm(data2)-norm(y2)) < 1e-9
                @test (@vecnorm_eq_eps idwt(y,wt) data stderr1*100)
                @test (@vecnorm_eq_eps idwt(y2,wt) data2 stderr2*100)
            end
        end
    end
end

@testset "Accuracy non-square" begin
    data2 = readdlm(joinpath(dirname(@__FILE__), "data", "filter2d_nonsquare_data.txt"))
    y2  = dwt(data2, wavelet(WT.haar), 1)
    ye2 = readdlm(joinpath(dirname(@__FILE__), "data", "filter2d_nonsquare_Haar0.txt"))
    stderr2 = 1e-9*sqrt(length(ye2))
    @test (@vecnorm_eq_eps y2 ye2 stderr2)
end

@testset "Lifting vs filter" begin
    n = 32
    stderr1 = 1e-10*sqrt(n)
    stderr2 = 1e-10*sqrt(n*n)
    # 1-D, 2-D, 3-D lifting vs filtering / inverse vs original tests
    @testset "wclass $wclass" for wclass in (WT.db1, WT.db2)
        wf = wavelet(wclass, WT.Filter)
        wls = wavelet(wclass, WT.Lifting)
        x = randn(n)
        x2 = copy(x)

        for L in (ndyadicscales(n),0,1,2)

            yf = dwt(x, wf, L)
            yls = dwt(x, wls, L)

            # filter vs lifting
            @test (@vecnorm_eq_eps yf yls stderr1)

            ytf = idwt(yf, wf, L)
            ytls = idwt(yls, wls, L)

            # inverse vs original
            @test (@vecnorm_eq_eps ytf x stderr1)
            @test (@vecnorm_eq_eps ytls x stderr1)
        end

        x = randn(n,n)
        x2 = copy(x)
        for L in (ndyadicscales(n),0,1,2)

            yf = dwt(x, wf, L)
            yls = dwt(x, wls, L)

            # filter vs lifting
            @test (@vecnorm_eq_eps yf yls stderr2)

            ytf = idwt(yf, wf, L)
            ytls = idwt(yls, wls, L)

            # inverse vs original
            @test (@vecnorm_eq_eps ytf x stderr2)
            @test (@vecnorm_eq_eps ytls x stderr2)
        end

        x = randn(n,n,n)
        x2 = copy(x)
        for L in (ndyadicscales(n),0,1,2)

            yf = dwt(x, wf, L)
            #if L == 5
            #    @profile dwt(x, wls, L)
            #    Profile.print()
            #end
            yls = dwt(x, wls, L)

            # filter vs lifting
            #@test (@vecnorm_eq_eps yf yls stderr2)
            @test yf ≈ yls rtol=stderr2

            ytf = idwt(yf, wf, L)
            ytls = idwt(yls, wls, L)

            # inverse vs original
            #@test (@vecnorm_eq_eps ytf x stderr2)
            #@test (@vecnorm_eq_eps ytls x stderr2)
            @test ytf ≈ x rtol=stderr2
            @test ytls ≈ x rtol=stderr2

        end
    end
end

@testset "Types and sizes" begin
    function makedwt(ft::Type, n, wf, L)
    	x0 = rand(-5:5, n)
    	x = zeros(ft, n)
    	copyto!(x,x0)
    	return (x, dwt(x, wf, L))
    end

    # 1-d
    n = 8
    wf = wavelet(WT.db2)
    L = 2

    x, y = makedwt(Float64, n, wf, L)
    @test Vector{Float64} == typeof(y) && length(y) == n
    x, y = makedwt(Float32, n, wf, L)
    @test Vector{Float32} == typeof(y) && length(y) == n
    x, y = makedwt(Int64, n, wf, L)
    @test Vector{typeof(float(x[1]))} == typeof(y) && length(y) == n
    x, y = makedwt(Int32, n, wf, L)
    @test Vector{typeof(float(x[1]))} == typeof(y) && length(y) == n
    @test dwt(x,wf) ≈ dwt(float(x),wf)

    # Complex types
    for wfc in (wavelet(WT.db2), wavelet(WT.db2, WT.Lifting))
        xc = zeros(ComplexF64, n)
        map!(i->rand(ComplexF64), xc, xc)
        yc = dwt(xc, wfc, L)
        @test Vector{ComplexF64} == typeof(yc) && length(yc) == n
        @test xc ≈ idwt(dwt(xc, wfc), wfc)
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
    @test dwt(x,wf) ≈ dwt(float(x),wf)

    # non-Array type
    wt = wavelet(WT.db2, WT.Lifting)
    x = randn(16, 16)
    xs = view(copy(x), 1:16, 1:16)
    @test dwt(x,wt,2) ≈ dwt(xs,wt,2)

    # util functions
    for class in (WT.haar, WT.db2, WT.cdf97)
        WT.class(class)
        WT.name(class)
        WT.vanishingmoments(class)
    end

    wt = wavelet(WT.db1, WT.Filter)
    @test length(wt) == 2
    @test wt.qmf*0.7 ≈ WT.scale(wt, 0.7).qmf

    # inplace methods
    wt = wavelet(WT.db1, WT.Filter)
    x = randn(8)
    y = similar(x)
    @test dwt!(y, x, wt) ≈ dwt(x, wt)

    wt = wavelet(WT.db1, WT.Lifting)
    x = randn(8)
    y = copy(x)
    @test dwt!(x, wt) ≈ dwt(y, wt)
end

struct wunknownt <: DiscreteWavelet{Float64} end
@testset "Errors" begin
    uwt = wunknownt()
    EE = Exception
    @test_throws EE dwt(randn(4),uwt)
    @test_throws EE dwt(randn(4,4),uwt)
    @test_throws EE wavelet(WT.Coiflet{33}())
    @test_throws EE wavelet("db2asdsad")
    @test_throws EE wavelet("db2", "ppppp")
end


# 2-d
function makedwt(ft::Type, n, wf, L)
    x0 = rand(-5:5, n)
    x = zeros(ft, n)
    copyto!(x,x0)
    return (x, dwt(x, wf, L))
end

n=8
wf = wavelet(WT.db2)
L = 2
sett = ((n,n),wf,L)
ft = Float64; x, y = makedwt(ft, sett...)
@test Array{ft,2} == typeof(y) && length(y) == n*n
ft = Float32; x, y = makedwt(ft, sett...)
@test Array{ft,2} == typeof(y) && length(y) == n*n
ft = Int64; x, y = makedwt(ft, sett...)
@test Array{typeof(float(x[1])),2} == typeof(y) && length(y) == n*n
ft = Int32; x, y = makedwt(ft, sett...)
@test Array{typeof(float(x[1])),2} == typeof(y) && length(y) == n*n
@test dwt(x,wf) ≈ dwt(float(x),wf)

# non-Array type
wt = wavelet(WT.db2, WT.Lifting)
x = randn(16, 16)
xs = view(copy(x), 1:16, 1:16)
@test dwt(x,wt,2) ≈ dwt(xs,wt,2)

# util functions
for class in (WT.haar, WT.db2, WT.cdf97)
    WT.class(class)
    WT.name(class)
    WT.vanishingmoments(class)
end
class = WT.db1
wt = wavelet(class, WT.Filter)
@test length(wt) == 2
@test wt.qmf*0.7 ≈ WT.scale(wt, 0.7).qmf

# inplace methods
class = WT.db1
wt = wavelet(class, WT.Filter)
x = randn(8)
y = similar(x)
@test dwt!(y, x, wt) ≈ dwt(x, wt)

wt = wavelet(class, WT.Lifting)
x = randn(8)
y = copy(x)
@test dwt!(x, wt) ≈ dwt(y, wt)

# continuous 1-d; different scalings should lead to different sizes, different boundary condtions shouldn't
@testset "Continuous Wavelet Transform" begin
    for xSize = (33, 67)
        for boundary = (WT.DEFAULT_BOUNDARY, WT.padded, WT.NaivePer)
            for s=1:2:8
                for wfc in (wavelet(WT.morl,s=s,boundary=boundary), wavelet(WT.dog0,s=s,boundary=boundary), wavelet(WT.paul4,s=s,boundary=boundary))
                    xc = rand(Float64,xSize)
                    yc = cwt(xc,wfc)
                    if typeof(wfc.waveType) <: Union{WT.Morlet, WT.Paul}
                        @test Array{ComplexF64,2}==typeof(yc)
                    else
                        @test Array{Float64,2}==typeof(yc)
                    end
                    nOctaves= log2(xSize) - wfc.averagingLength; 
                    nWaveletsInOctave = reverse([max(1, round(Int, s / x^(1))) for
                                                 x=1:round(Int, nOctaves)])
                    totalWavelets = round(Int, sum(nWaveletsInOctave))
                    @test size(yc) == (xSize, totalWavelets+1)
                end
            end
        end
    end
end
# TODO: test actual values
#       test averaging types
#            various extra dimensions


@testset "WPT" begin
    wf = wavelet(WT.db2, WT.Filter)
    x = randn(16)

    L = 1
    wp = wpt(x,wf,L)
    dw = dwt(x,wf,L)
    @test wp ≈ dw
    @test iwpt(wp,wf,L) ≈ x

    L = 2
    wp = wpt(x,wf,L)
    dw = dwt(x,wf,L)
    dw2 = copy(dw)
    dw2[9:end] = dwt(dw[9:end],wf,1)
    @test dw[1:8] ≈ wp[1:8]
    @test dw2 ≈ wp
    @test iwpt(wp,wf,L) ≈ x

    L = 3
    wp = wpt(x,wf,L)
    dw = dwt(x,wf,L)
    @test dw[1:4] ≈ wp[1:4]
    @test dwt(dw2[5:8],wf,1) ≈ wp[5:8]
    @test dwt(dw2[9:12],wf,1) ≈ wp[9:12]
    @test dwt(dw2[13:16],wf,1) ≈ wp[13:16]
    @test iwpt(wp,wf,L) ≈ x

    wl = wavelet(WT.db2, WT.Lifting)
    x = randn(128)
    @test iwpt(wpt(x,wf),wf) ≈ x
    @test iwpt(wpt(x,wl),wl) ≈ x
    @test wpt(x,wl,1) ≈ wpt(x,wf,1)
    @test wpt(x,wl,2) ≈ wpt(x,wf,2)
    @test wpt(x,wl,4) ≈ wpt(x,wf,4)
    @test wpt(x,wl) ≈ wpt(x,wf)

    wl = wavelet(WT.db2, WT.Lifting)
    n = 128
    x = randn(n)
    for L = 0:ndyadicscales(n)
        @test wpt(x, wl, maketree(n, L, :dwt)) ≈ dwt(x, wl, L)
        @test iwpt(x, wl, maketree(n, L, :dwt)) ≈ idwt(x, wl, L)
        @test wpt(x, wf, maketree(n, L, :dwt)) ≈dwt(x, wf, L)
        @test iwpt(x, wf, maketree(n, L, :dwt)) ≈ idwt(x, wf, L)
    end

    # non-dyadic tests
    wt = wavelet(WT.db2, WT.Filter)
    n = 5*8
    x = randn(n)
    for L in 0:maxtransformlevels(n)
        for wt in (wavelet(WT.db2, WT.Filter), wavelet(WT.db2, WT.Lifting))
            @test wpt(x, wt, maketree(n, L, :dwt)) ≈ dwt(x, wt, L)
            @test iwpt(x, wt, maketree(n, L, :dwt)) ≈ idwt(x, wt, L)
        end
    end
end

@testset "MODWT" begin
    wf = wavelet(WT.db4)
    # power-of-two
    x = randn(128)
    W = modwt(x, wf)
    x_back = imodwt(W, wf)
    @test x ≈ x_back

    # non power-of-two
    x = cumsum(randn(129))
    W = modwt(x, wf)
    x_back = imodwt(W, wf)
    @test x ≈ x_back
    @test size(W) == (length(x), maxmodwttransformlevels(x)+1)

    # less than max number of levels
    L = 4
    Wl = modwt(x, wf, L)
    @test W[:, 1:L-1] ≈ Wl[:, 1:L-1]
end

@testset "Wavelet types" begin
    wt = wavelet(WT.db2, WT.Filter)
    @test length(wt) == 4
    same = WT.name(wt) == "db2"
    @test same

    wt = wavelet(WT.db2, WT.Lifting)
    @test WT.name(wt) == "db2"
end
