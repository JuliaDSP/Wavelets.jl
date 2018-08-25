
@testset "Threshold" begin # TODO @test something
    x = randn(200)*2
    for th in (BiggestTH(), HardTH(), SoftTH(), SemiSoftTH(), SteinTH())
        threshold(x, th, 2)
    end
    for th in (PosTH(), NegTH())
        threshold(x, th)
    end
end

@testset "Denoising" begin
    vs = VisuShrink(10)
    n = 2^8
    x0 = testfunction(n, "Doppler")
    x = x0 + 0.05*randn(n)
    y = denoise(x, TI=true)
    y = denoise(x, TI=true, nspin = 8)
    y = denoise(x, TI=false)
    y = denoise(x, nothing)
    y = denoise(randn(32,32), TI=true)
end

@testset "Best basis" begin
    wt = wavelet(WT.db4)
    x = sin.(4*range(0, stop=2*pi-eps(), length=1024))
    tree = bestbasistree(x, wt)
    xtb = wpt(x, wt, tree)
    @test iwpt(xtb, wt, tree) ≈ x

    x = sin.(4*range(0, stop=2*pi-eps(), length=5*64)) # non-dyadic
    tree = bestbasistree(x, wt)
    xtb = wpt(x, wt, tree)
    @test iwpt(xtb, wt, tree) ≈ x
end

@testset "Matching pursuit" begin
    N = 128
    M = 64
    y=randn(N)
    Arr = randn(M,N)
    func(a::AbstractVector) = Arr*a
    funct(a::AbstractVector) = Arr'*a
    x = func(y)
    matchingpursuit(x, func, funct, 0.1)
end
