
include("setup_1d.jl")

println("fft (N=",N,"), (FFTW)")
f(x0) = for i = 1:tn; fft(x0); end
f(x0);
@time f(x0);

