Wavelets
=========

[![Build Status](https://travis-ci.org/gummif/Wavelets.jl.svg?branch=master)](https://travis-ci.org/gummif/Wavelets.jl)

A [Julia](https://github.com/JuliaLang/julia) module for very fast wavelet transforms (orthogonal, periodic, 1D, 2D).

Rouchly 20x speedup and 50x less memory usage than [this](https://github.com/tomaskrehlik/Wavelets) implementation of `dwt`. Loosely inspired by [this](https://github.com/tomaskrehlik/Wavelets) and [this](http://statweb.stanford.edu/~wavelab). See benchmarks and a todo list below.

Written by Gudmundur Adalsteinsson. See license in LICENSE.md.

Usage
---------

A rough idea of the API:

```julia
# set up filter (Periodic, Orthogonal, 4 vanishing moments)
wf = POfilter("Coiflet", 4)
# or do it this way
wf = POfilter("coif4")
# xt is a 5 level transform of vector x
xt = fwt(x, 5, wf)
# x2 is approx. equal to x
x2 = iwt(xt, 5, wf)

# scaling filters is easy
wf1 = scale(POfilter("haar"), 1/sqrt(2))
# signals can be transformed inplace with a low-level command
# requiring very little memory allocation (especially for L=1)
dwt!(xt, x, L, wf, true)
```


Benchmarks
---------

Timing of `fwt`, `iwt` (using `db6` filter of length 12) and `fft`. The wavelet transforms are faster and use less memory than `fft` in all cases.

```julia
# 100 iterations
fwt (N=1048576), 18 levels
elapsed time: 17.024847664 seconds (1258560248 bytes allocated, 1.34% gc time)
iwt (N=1048576), 18 levels
elapsed time: 11.731993823 seconds (1258562648 bytes allocated, 1.93% gc time)
fft (N=1048576), (FFTW)
elapsed time: 29.745145114 seconds (5872236248 bytes allocated, 4.44% gc time)

# 1000 iterations
fwt (N=32768), 13 levels
elapsed time: 5.30228259 seconds (395317512 bytes allocated, 1.37% gc time)
iwt (N=32768), 13 levels
elapsed time: 3.581545857 seconds (395341512 bytes allocated, 2.01% gc time)
fft (N=32768), (FFTW)
elapsed time: 8.85191674 seconds (1837005496 bytes allocated, 5.72% gc time)

# 10000 iterations
fwt (N=512), 7 levels
elapsed time: 0.932289745 seconds (77037512 bytes allocated, 1.15% gc time)
iwt (N=512), 7 levels
elapsed time: 0.642051702 seconds (77277512 bytes allocated, 1.69% gc time)
fft (N=512), (FFTW)
elapsed time: 1.457472828 seconds (306860616 bytes allocated, 6.14% gc time)
```

For 2D transforms:
```julia
# 10 iterations
fwt (N=1024 x 1024), 8 levels
elapsed time: 2.422658432 seconds (98389208 bytes allocated, 0.79% gc time)
iwt (N=1024 x 1024), 8 levels
elapsed time: 2.048636659 seconds (95275608 bytes allocated, 0.98% gc time)
fft (N=1024 x 1024), (FFTW)
elapsed time: 2.945895417 seconds (587236728 bytes allocated, 4.69% gc time)
```

By using the low-level function `dwt!` and pre-allocating temporary arrays, significant gains can be made in terms of memory usage (and a little speedup). This is useful when transforming multiple signals.
```julia
# 1000 iterations
dwt! (N=32768), 13 levels
elapsed time: 5.100670485 seconds (132013512 bytes allocated, 0.62% gc time)
fwt (N=32768), 13 levels
elapsed time: 5.151621451 seconds (395317512 bytes allocated, 1.62% gc time)
```



TODO
---------

* Boundary orthogonal wavelets
* Biorthogonal wavelets
* Lifting scheme wavelets
* Interpolating wavelets
* Thresholding/denoising functions
* Best M-term approx. and sparsity utilities
* Wavelet packets
* Continuous wavelets
* Visualization functions
* Wavelet scalogram



