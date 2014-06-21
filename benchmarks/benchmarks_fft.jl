



N = 1024*32;
x0 = rand(N);

tn = 1000

println("fft (N=",N,"), (FFTW)")
f(x0) = for i = 1:tn; fft(x0); end
f(x0);
@time f(x0);



