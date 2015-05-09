
include("setup_1d.jl")

x0c = complex(x0)

println("fft! (N=",N,"), (FFTW)")
f(x0c) = for i = 1:tn; fft!(x0c); end
f(x0c);
@time f(x0c);

