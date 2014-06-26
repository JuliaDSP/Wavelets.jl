Wavelets
=========

[![Build Status](https://travis-ci.org/gummif/Wavelets.jl.svg?branch=master)](https://travis-ci.org/gummif/Wavelets.jl)

A [Julia](https://github.com/JuliaLang/julia) package for very fast wavelet transforms (1D, 2D, by filtering or lifting).

Rouchly 20x speedup and 50x less memory usage than [this](https://github.com/tomaskrehlik/Wavelets) implementation of `dwt`. Loosely inspired by [this](https://github.com/tomaskrehlik/Wavelets) and [this](http://statweb.stanford.edu/~wavelab). See benchmarks and a todo list below.

* 1st generation wavelets using filter banks (periodic and orthogonal). Filters are included for the following types: Haar, Daubechies, Coiflet, Symmlet, Battle-Lemarie, Beylkin, Vaidyanathan.

* 2nd generation wavelets by lifting (periodic and general type including orthogonal and biorthogonal). Included lifting schemes are currently only for Haar and Daubechies (under development). A new lifting scheme can be easily constructed by users. The current implementation of the lifting transforms is 10x faster than the filter transforms.

Written by Gudmundur Adalsteinsson (c) 2014. See license (MIT) in LICENSE.md.

Usage
---------

A rough idea of the API:

```julia
# set up filter (Periodic, Orthogonal, 4 vanishing moments)
wt = POfilter("Coiflet", 4)
# or do it this way
wt = POfilter("coif4")
# the input object determines the transform type 
# wt now contains instructions for a periodic lifting scheme
wt = GPLS("coif4")
# xt is a 5 level transform of vector x
xt = fwt(x, 5, wt)
# x2 is approx. equal to x
x2 = iwt(xt, 5, wt)

# scaling filters is easy
wt1 = scale(POfilter("haar"), 1/sqrt(2))
# signals can be transformed inplace with a low-level command
# requiring very little memory allocation (especially for L=1)
dwt!(x, L, wt, true)      # inplace by lifting
dwt!(xt, x, L, wt, true)  # write to xt by filtering
```


Benchmarks
---------

Timing of `fwt` (using `db2` filter of length 4) and `fft`. The wavelet transforms are faster and use less memory than `fft` in all cases. `fwt` by lifting is currently an order of magnitude faster than by filtering.

```julia
# 10 iterations
fwt by filtering (N=1048576), 18 levels
elapsed time: 1.262672396 seconds (125866088 bytes allocated, 2.35% gc time)
fwt by lifting (N=1048576), 18 levels
elapsed time: 0.153162937 seconds (104927608 bytes allocated, 13.18% gc time)
fft (N=1048576), (FFTW)
elapsed time: 2.836224168 seconds (587236088 bytes allocated, 4.83% gc time)
```

For 2D transforms by filtering (using `db6` of length 12):
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
elapsed time: 5.081105993 seconds (701512 bytes allocated)
fwt (N=32768), 13 levels
elapsed time: 5.151621451 seconds (395317512 bytes allocated, 1.62% gc time)
```



TODO
---------

* Boundary orthogonal wavelets
* Biorthogonal wavelets
* Thresholding/denoising functions
* Best M-term approx. and sparsity utilities
* Wavelet packets
* Continuous wavelets
* Visualization functions
* Wavelet scalogram



