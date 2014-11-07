
include("setup_2d.jl")

println("fft (N=",N,"x",N,"), (FFTW)")
f(x0) = for i = 1:tn; fft(x0); end
f(x0);
@time f(x0);

