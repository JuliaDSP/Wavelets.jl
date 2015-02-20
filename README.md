<img src="wavelets.png" alt="Wavelets">

---------

[![Build Status](https://travis-ci.org/JuliaDSP/Wavelets.jl.svg?branch=master)](https://travis-ci.org/JuliaDSP/Wavelets.jl)
[![Coverage Status](https://coveralls.io/repos/JuliaDSP/Wavelets.jl/badge.svg?branch=master)](https://coveralls.io/r/JuliaDSP/Wavelets.jl?branch=master)

A [Julia](https://github.com/JuliaLang/julia) package for fast wavelet transforms (1-D, 2-D, by filtering or lifting). The package includes discrete wavelet transforms, column-wise discrete wavelet transforms, and wavelet packet transforms.

* 1st generation wavelets using filter banks (periodic and orthogonal). Filters are included for the following types: Haar, Daubechies, Coiflet, Symmlet, Battle-Lemarie, Beylkin, Vaidyanathan.

* 2nd generation wavelets by lifting (periodic and general type including orthogonal and biorthogonal). Included lifting schemes are currently only for Haar and Daubechies (under development). A new lifting scheme can be easily constructed by users. The current implementation of the lifting transforms is 2x faster than the filter transforms.

* Thresholding, best basis and denoising functions, e.g. TI denoising by cycle spinning, best basis for WPT, noise estimation, matching pursuit. See example code and image below.

* Wavelet utilities e.g. indexing and size calculation, scaling and wavelet functions computation, test functions, up and down sampling, filter mirrors, coefficient counting, inplace circshifts, and more.

* Plotting/visualization utilities for 1-D and 2-D signals.

Loosely inspired by [this](https://github.com/tomaskrehlik/Wavelets) and [this](http://statweb.stanford.edu/~wavelab). 

See license (MIT) in LICENSE.md.


Usage
---------

Install via the package manager and load with `using`

```julia
julia> Pkg.add("Wavelets")
julia> using Wavelets
```


API
---------

#### Wavelet transforms and types
```julia
# Transform Type construction
wavelet(w::WT.WaveletClass, boundary::WaveletBoundary=Periodic)  # defaults to filter
waveletfilter(w::WT.WaveletClass, boundary::WaveletBoundary=Periodic)
waveletls(w::WT.WaveletClass, boundary::WaveletBoundary=Periodic)
# DWT (discrete wavelet transform)
dwt(x::AbstractArray, wt::DiscreteWavelet, L::Integer=nscales(x))
idwt(x::AbstractArray, wt::DiscreteWavelet, L::Integer=nscales(x))
dwt!(y::AbstractArray, x::AbstractArray, filter::OrthoFilter, L::Integer, fw::Bool)
dwt!(y::AbstractArray, scheme::GLS, L::Integer, fw::Bool)
# DWTC (discrete wavelet transform along dimension td (default is last dim.))
dwtc(x::AbstractArray, wt::DiscreteWavelet, L::Integer=nscales(x), td::Integer=ndims(x))
idwtc(x::AbstractArray, wt::DiscreteWavelet, L::Integer=nscales(x), td::Integer=ndims(x))
# WPT (wavelet packet transform)
# Ltree can be L::Integer=nscales(x) or tree::BitVector=maketree(length(y), L, :full)
wpt(x::AbstractArray, wt::DiscreteWavelet, Ltree)
iwpt(x::AbstractArray, wt::DiscreteWavelet, Ltree)
wpt!(y::AbstractArray, x::AbstractArray, filter::OrthoFilter, Ltree, fw::Bool)
wpt!(y::AbstractArray, scheme::GLS, Ltree, fw::Bool)
```

#### Wavelet classes

The module WT contains the types for wavelet classes. The module defines constants of the form e.g. `WT.coif4` as shortcuts for `WT.Coiflet{4}()`.
The numbers for orthogonal wavelets indicate the number vanishing moments of the wavelet function. The length of a `wt::OrthoFilter` can be obtained with `length(wt)`.

| Class Type | Namebase | Supertype | Numbers |
|:------- |:------ |:----- |:----- |
| `Haar` | haar | `OrthoWaveletClass` |   |
| `Coiflet` | coif | `OrthoWaveletClass` | 2:2:8 |
| `Daubechies` | db | `OrthoWaveletClass` | 1:Inf |
| `Symlet` | sym | `OrthoWaveletClass` | 4:10 |
| `Battle` | batt | `OrthoWaveletClass` | 2:2:6
| `Beylkin` | beyl | `OrthoWaveletClass` |  |
| `Vaidyanathan` | vaid | `OrthoWaveletClass` |  |
| `CDF` | cdf | `BiOrthoWaveletClass` | (9,7) |

Class information functions
```julia
class(::WaveletClass) ::ASCIIString         # class full name
name(::WaveletClass) ::ASCIIString          # type short name
vanishingmoments(::WaveletClass)            # vanishing moments of wavelet function
```

Examples
---------

```julia
# the simplest way to transform a signal x is
xt = dwt(x, wavelet(WT.db2))

# the transform type can be more explicitly specified
# set up wavelet type (filter, Periodic, Orthogonal, 4 vanishing moments)
wt = wavelet(WT.Coiflet{4}(), Periodic)  # or
wt = waveletfilter(WT.coif4)
# which is equivalent to 
wt = OrthoFilter(WT.coif4)
# the object wt determines the transform type 
# wt now contains instructions for a periodic biorthogonal CDF 9/7 lifting scheme
wt = waveletls(WT.cdf97)
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
y = dwt(x, waveletls(WT.cdf97))
d,l = wplotdots(y, 0.1, n)
A = wplotim(y)
# plotting utilities 2-d
img = imread("lena.png")
x = permutedims(img.data, [ndims(img.data):-1:1])
L = 2
xts = wplotim(x, L, waveletfilter(WT.db3))
```

See [Bumps](/example/transform1d_bumps.png) and [Lena](/example/transform2d_lena.jpg) for plot images.


Thesholding
---------

The `Wavelets.Threshold` module includes the following utilities

* denoising (VisuShrink, translation invariant (TI))
* best basis for WPT
* noise estimator
* matching pursuit
* threshold functions (see table for types)

#### API
```julia
# threshold types with parameter
threshold!(x::AbstractArray, TH::THType, t::Real)
threshold(x::AbstractArray, TH::THType, t::Real)
# without parameter (PosTH, NegTH)
threshold!(x::AbstractArray, TH::THType)
threshold(x::AbstractArray, TH::THType)
# denoising
denoise(x::AbstractArray,
        wt=DEFAULT_WAVELET;
        L::Int=min(maxtranformlevels(x),6),
        dnt=VisuShrink(size(x,1)),
        estnoise::Function=noisest, 
        TI=false,
        nspin=tuple([8 for i=1:ndims(x)]...) )
# entropy
coefentropy(x, et::Entropy, nrm)
# best basis for WPT limited to active inital tree nodes
bestbasistree(y::AbstractVector, wt::DiscreteWavelet, 
        L::Integer=nscales(y), et::Entropy=ShannonEntropy() )
bestbasistree(y::AbstractVector, wt::DiscreteWavelet, 
        tree::BitVector, et::Entropy=ShannonEntropy() )
```

| Type | Details |
|:------- |:------ |
| **Thresholding** | ` <: THType` |
| `HardTH` | hard thresholding |
| `SoftTH` | soft threshold |
| `SemiSoftTH` | semisoft thresholding |
| `SteinTH` | stein thresholding |
| `PosTH` | positive thresholding |
| `NegTH` | negative thresholding |
| `BiggestTH` | biggest m-term (best m-term) approx. |
| **Denoising** | ` <: DNFT` |
| `VisuShrink` | VisuShrink denoising |
| **Entropy** | ` <: Entropy` |
| `ShannonEntropy` | `-v^2*log(v^2)` (Coifman-Wickerhauser) |
| `LogEnergyEntropy` | `-log(v^2)` |



#### Examples
Find best basis tree for WPT:
```julia
wt = wavelet(WT.db4)
x = sin(4*linspace(0,2*pi-eps(),1024))
tree = bestbasistree(x, wt)
xtb = wpt(x, wt, tree)
xt = dwt(x, wt)
# estimate sparsity in ell_1 norm
vecnorm(xtb, 1)  # best basis wpt
vecnorm(xt, 1)   # regular dwt
```
Denoising:
```julia
n = 2^11;
x0 = testfunction(n,"Doppler")
x = x0 + 0.05*randn(n)
y = denoise(x,TI=true)
```
![Doppler](/example/denoise_doppler.png)


Benchmarks
---------

Timing of `dwt` (using `db2` filter of length 4) and `fft`. The lifting wavelet transforms are faster and use less memory than `fft` in 1-D and 2-D. `dwt` by lifting is currently 2x faster than by filtering.

```julia
# 10 iterations
dwt by filtering (N=1048576), 20 levels
elapsed time: 0.247907616 seconds (125861504 bytes allocated, 8.81% gc time)
dwt by lifting (N=1048576), 20 levels
elapsed time: 0.131240966 seconds (104898544 bytes allocated, 17.48% gc time)
fft (N=1048576), (FFTW)
elapsed time: 0.487693289 seconds (167805296 bytes allocated, 8.67% gc time)
```

For 2-D transforms (using a `db4` filter and CDF 9/7 lifting scheme):
```julia
# 10 iterations
dwt by filtering (N=1024x1024), 10 levels
elapsed time: 0.773281141 seconds (85813504 bytes allocated, 2.87% gc time)
dwt by lifting (N=1024x1024), 10 levels
elapsed time: 0.317705928 seconds (88765424 bytes allocated, 3.44% gc time)
fft (N=1024x1024), (FFTW)
elapsed time: 0.577537263 seconds (167805888 bytes allocated, 5.53% gc time)
```

By using the low-level function `dwt!` and pre-allocating temporary arrays, significant gains can be made in terms of memory usage (and some speedup). This is useful when transforming multiple times.

To-do list
---------

* Transforms for non-square 2-D signals
* Boundary extensions (other than periodic)
* Boundary orthogonal wavelets
* Define more lifting schemes
* WPT in 2-D
* Stationary transform
* Continuous wavelets
* Wavelet scalogram
* DWT in 3-D


