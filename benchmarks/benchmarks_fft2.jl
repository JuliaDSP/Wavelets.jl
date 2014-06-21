



N = 1024;
x0 = rand(N,N);

tn = 10

println("fft (N=",N,"x",N,"), (FFTW)")
f(x0) = for i = 1:tn; fft(x0); end
f(x0);
@time f(x0);



