<img src="wavelets.png" alt="Wavelets">

---------

[![Build Status](https://travis-ci.org/JuliaDSP/Wavelets.jl.svg?branch=master)](https://travis-ci.org/JuliaDSP/Wavelets.jl)
[![Coverage Status](https://coveralls.io/repos/JuliaDSP/Wavelets.jl/badge.png?branch=master)](https://coveralls.io/r/JuliaDSP/Wavelets.jl?branch=master)

A [Julia](https://github.com/JuliaLang/julia) package for fast wavelet transforms (1-D, 2-D, by filtering or lifting).

* 1st generation wavelets using filter banks (periodic and orthogonal). Filters are included for the following types: Haar, Daubechies, Coiflet, Symmlet, Battle-Lemarie, Beylkin, Vaidyanathan.

* 2nd generation wavelets by lifting (periodic and general type including orthogonal and biorthogonal). Included lifting schemes are currently only for Haar and Daubechies (under development). A new lifting scheme can be easily constructed by users. The current implementation of the lifting transforms is 10x faster than the filter transforms.

* Denoising and thresholding functions and utilities, e.g. TI denoising by cycle spinning, noise estimation, matching pursuit. See example code and image below.

* Wavelet utilities e.g. indexing and size calculation, scaling and wavelet functions computation, test functions, up and down sampling, filter mirrors, coefficient counting, inplace circshifts, and more.

* Plotting/visualization utilities for 1-D and 2-D signals.

Roughly 20x speedup and 50x less memory usage than [this](https://github.com/tomaskrehlik/Wavelets) implementation of `dwt`. Loosely inspired by [this](https://github.com/tomaskrehlik/Wavelets) and [this](http://statweb.stanford.edu/~wavelab). 

See license (MIT) in LICENSE.md.


Usage
---------

Install via the package manager and load with `using`

```julia
julia> Pkg.add("Wavelets")
julia> using Wavelets
```

A rough idea of the API:

```julia
# the simplest way to transform a signal x is
xt = dwt(x, wavelet("db2"))

# the transform type can be more explicitly specified
# set up wavelet type (filter, Periodic, Orthogonal, 4 vanishing moments)
wt = wavelet("Coiflet", 4)  # or
wt = wavelet("Coiflet", 4, transform="filter", boundary="per")  # or
wt = waveletfilter("coif4")
# which is equivalent to 
wt = OrthoFilter("coif4")
# the object wt determines the transform type 
# wt now contains instructions for a periodic biorthogonal CDF 9/7 lifting scheme
wt = waveletls("cdf9/7")
# xt is a 5 level transform of vector x
xt = dwt(x, wt, 5)
# inverse tranform
xti = idwt(xt, wt, 5)
# a full transform
xt = dwt(x, wt)

# scaling filters is easy
wt = scale(wt, 1/sqrt(2))
# signals can be transformed inplace with a low-level command
# requiring very little memory allocation (especially for L=1 for filters)
dwt!(x, wt, L, true)      # inplace (lifting)
dwt!(xt, x, wt, L, true)  # write to xt (filter)

# denoising with default parameters (VisuShrink hard thresholding)
x0 = testfunction(n, "HeaviSine")
x = x0 + 0.3*randn(n)
y = denoise(x)

# plotting utilities 1-d (see images and code in /example)
x = testfunction(n, "Bumps")
y = dwt(x, waveletls("cdf9/7"))
d,l = wplotdots(y, 0.1, n)
A = wplotim(y)
# plotting utilities 2-d
img = imread("lena.png")
x = permutedims(img.data, [ndims(img.data):-1:1])
L = 2
xts = wplotim(x, L, waveletfilter("db3"))
```

![Bumps](/example/transform1d_bumps.png)

![Lena](/example/transform2d_lena.jpg)

#### Wavelet class information

| Long name | Short | Type | Numbers |
|:------- |:------ |:----- |:----- |
| `Haar` | `haar` | Ortho |   |
| `Coiflet` | `coif` | Ortho | 2:2:8 |
| `Daubechies` | `db` | Ortho | 1:10 |
| `Symmlet`/`Symlet` | `sym` | Ortho | 4:10 |
| `Battle` | `batt` | Ortho | 2:2:6
| `Beylkin` | `beyl` | Ortho |  |
| `Vaidyanathan` | `vaid` | Ortho |  |
| `CDF` | `cdf` | BiOrtho | "9/7" |


Benchmarks
---------

Timing of `dwt` (using `db2` filter of length 4) and `fft`. The wavelet transforms are faster and use less memory than `fft` in all cases. `dwt` by lifting is currently an order of magnitude faster than by filtering.

```julia
# 10 iterations
dwt by filtering (N=1048576), 18 levels
elapsed time: 1.262672396 seconds (125866088 bytes allocated, 2.35% gc time)
dwt by lifting (N=1048576), 18 levels
elapsed time: 0.153162937 seconds (104927608 bytes allocated, 13.18% gc time)
fft (N=1048576), (FFTW)
elapsed time: 2.836224168 seconds (587236088 bytes allocated, 4.83% gc time)
```

For 2D transforms (the `dwt` using a CDF 9/7 lifting scheme):
```julia
# 10 iterations
dwt by lifting (N=1024x1024), 10 levels
elapsed time: 0.952415385 seconds (109075368 bytes allocated, 2.14% gc time)
fft (N=1024 x 1024), (FFTW)
elapsed time: 2.945895417 seconds (587236728 bytes allocated, 4.69% gc time)
```

By using the low-level function `dwt!` and pre-allocating temporary arrays, significant gains can be made in terms of memory usage (and a little speedup). This is useful when transforming multiple signals.
```julia
# 1000 iterations
dwt! (N=32768), 13 levels
elapsed time: 5.081105993 seconds (701512 bytes allocated)
dwt (N=32768), 13 levels
elapsed time: 5.151621451 seconds (395317512 bytes allocated, 1.62% gc time)
```

Denoising and thesholding
---------

The `Wavelets.Threshold` module includes the following utilities

* denoising (VisuShrink, translation invariant (TI))
* noise estimator
* matching pursuit
* hard threshold
* soft threshold
* semisoft threshold
* stein threshold
* biggest m-term (best m-term) approximation
* positive and negative thresholds

Example of usage:
```julia
n = 2^11;
x0 = testfunction(n,"Doppler")
x = x0 + 0.05*randn(n)
y = denoise(x,TI=true)
```
![Doppler](/example/denoise_doppler.png)

To-do list
---------

* don't use sub arrays in 2d
* Boundary orthogonal wavelets
* Define more lifting schemes
* Redundant transforms and wavelet packets
* Continuous wavelets
* Wavelet scalogram



