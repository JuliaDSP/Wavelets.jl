#!/usr/bin/env julia

using BenchmarkTools
using GPUEnv
using Printf
using Random
using Statistics
using Wavelets

GPUEnv.activate(persist = true, include_jlarrays = false, only_first = true)

const BACKENDS = gpu_backends()

if isempty(BACKENDS)
    println("No functional native GPU backend detected. Skipping GPU benchmark.")
    exit(0)
end

const GPU_BACKEND = first(BACKENDS)
Random.seed!(0)

function bench_time(f; samples=10, evals=1)
    trial = run(@benchmarkable begin
        local result = $f()
        synchronize_backend(GPU_BACKEND)
        result
    end samples=samples evals=evals)
    return median(trial.times) / 1e9
end

function print_header(title)
    println()
    println(repeat("=", 78))
    println(title)
    println(repeat("=", 78))
    @printf("%-16s %12s %12s %10s\n", "Size", "CPU (ms)", "GPU (ms)", "Speedup")
end

function print_row(label, cpu_t, gpu_t)
    speedup = cpu_t / gpu_t
    @printf("%-16s %12.3f %12.3f %10.2fx\n", label, cpu_t * 1e3, gpu_t * 1e3, speedup)
    return speedup
end

results = Dict{String, Vector{Float64}}()
record!(k, v) = push!(get!(results, k, Float64[]), v)

function benchmark_case(group, label, cpu_f, gpu_f)
    cpu_t = bench_time(cpu_f)
    gpu_t = bench_time(gpu_f)
    record!(group, print_row(label, cpu_t, gpu_t))
end

println("Wavelets.jl GPU benchmark")
println("Backend: ", GPU_BACKEND.name)

print_header("1D Filter DWT")
wtf = wavelet(WT.db4)
for p in (16, 18, 20)
    n = 2^p
    x = randn(Float32, n)
    xg = to_gpu(GPU_BACKEND, x)
    L = min(8, maxtransformlevels(n))
    yg = dwt(xg, wtf, L)
    yc = dwt(x, wtf, L)
    benchmark_case("1D filter dwt", "2^$p", () -> dwt(x, wtf, L), () -> dwt(xg, wtf, L))
end

print_header("1D Filter IDWT")
for p in (16, 18, 20)
    n = 2^p
    x = randn(Float32, n)
    xg = to_gpu(GPU_BACKEND, x)
    L = min(8, maxtransformlevels(n))
    yg = dwt(xg, wtf, L)
    yc = dwt(x, wtf, L)
    benchmark_case("1D filter idwt", "2^$p", () -> idwt(yc, wtf, L), () -> idwt(yg, wtf, L))
end

print_header("1D Filter WPT")
for p in (14, 16, 18)
    n = 2^p
    x = randn(Float32, n)
    xg = to_gpu(GPU_BACKEND, x)
    benchmark_case("1D filter wpt", "2^$p", () -> wpt(x, wtf), () -> wpt(xg, wtf))
end

print_header("1D Filter IWPT")
for p in (14, 16, 18)
    n = 2^p
    x = randn(Float32, n)
    xg = to_gpu(GPU_BACKEND, x)
    wp_cpu = wpt(x, wtf)
    wp_gpu = wpt(xg, wtf)
    benchmark_case("1D filter iwpt", "2^$p", () -> iwpt(wp_cpu, wtf), () -> iwpt(wp_gpu, wtf))
end

print_header("2D Filter DWT")
for n in (512, 1024, 2048)
    x = randn(Float32, n, n)
    xg = to_gpu(GPU_BACKEND, x)
    L = min(4, maxtransformlevels(n))
    yg = dwt(xg, wtf, L)
    yc = dwt(x, wtf, L)
    benchmark_case("2D filter dwt", "$(n)x$(n)", () -> dwt(x, wtf, L), () -> dwt(xg, wtf, L))
end

print_header("2D Filter IDWT")
for n in (512, 1024, 2048)
    x = randn(Float32, n, n)
    xg = to_gpu(GPU_BACKEND, x)
    L = min(4, maxtransformlevels(n))
    yg = dwt(xg, wtf, L)
    yc = dwt(x, wtf, L)
    benchmark_case("2D filter idwt", "$(n)x$(n)", () -> idwt(yc, wtf, L), () -> idwt(yg, wtf, L))
end

print_header("3D Filter DWT")
for n in (32, 64, 128)
    x = randn(Float32, n, n, n)
    xg = to_gpu(GPU_BACKEND, x)
    L = min(3, maxtransformlevels(n))
    yg = dwt(xg, wtf, L)
    yc = dwt(x, wtf, L)
    benchmark_case("3D filter dwt", "$(n)^3", () -> dwt(x, wtf, L), () -> dwt(xg, wtf, L))
end

print_header("3D Filter IDWT")
for n in (32, 64, 128)
    x = randn(Float32, n, n, n)
    xg = to_gpu(GPU_BACKEND, x)
    L = min(3, maxtransformlevels(n))
    yg = dwt(xg, wtf, L)
    yc = dwt(x, wtf, L)
    benchmark_case("3D filter idwt", "$(n)^3", () -> idwt(yc, wtf, L), () -> idwt(yg, wtf, L))
end

print_header("1D Lifting DWT")
wtl = wavelet(WT.cdf97, WT.Lifting)
for p in (16, 18, 20)
    n = 2^p
    x = randn(Float32, n)
    xg = to_gpu(GPU_BACKEND, x)
    L = min(8, maxtransformlevels(n))
    yg = dwt(xg, wtl, L)
    yc = dwt(x, wtl, L)
    benchmark_case("1D lifting dwt", "2^$p", () -> dwt(x, wtl, L), () -> dwt(xg, wtl, L))
end

print_header("1D Lifting IDWT")
wtl = wavelet(WT.cdf97, WT.Lifting)
for p in (16, 18, 20)
    n = 2^p
    x = randn(Float32, n)
    xg = to_gpu(GPU_BACKEND, x)
    L = min(8, maxtransformlevels(n))
    yg = dwt(xg, wtl, L)
    yc = dwt(x, wtl, L)
    benchmark_case("1D lifting idwt", "2^$p", () -> idwt(yc, wtl, L), () -> idwt(yg, wtl, L))
end

print_header("1D Lifting WPT")
wtl = wavelet(WT.cdf97, WT.Lifting)
for p in (16, 18, 20)
    n = 2^p
    x = randn(Float32, n)
    xg = to_gpu(GPU_BACKEND, x)
    L = min(8, maxtransformlevels(n))
    yg = dwt(xg, wtl, L)
    yc = dwt(x, wtl, L)
    benchmark_case("1D lifting wpt", "2^$p", () -> wpt(x, wtl), () -> wpt(xg, wtl))
end

print_header("1D Lifting IWPT")
wtl = wavelet(WT.cdf97, WT.Lifting)
for p in (16, 18, 20)
    n = 2^p
    x = randn(Float32, n)
    xg = to_gpu(GPU_BACKEND, x)
    wp_cpu = wpt(x, wtl)
    wp_gpu = wpt(xg, wtl)
    benchmark_case("1D lifting iwpt", "2^$p", () -> iwpt(wp_cpu, wtl), () -> iwpt(wp_gpu, wtl))
end

print_header("1D MODWT")
wtm = wavelet(WT.db4)
for p in (14, 16, 18)
    n = 2^p
    x = randn(Float32, n)
    xg = to_gpu(GPU_BACKEND, x)
    L = min(6, maxmodwttransformlevels(n))
    benchmark_case("1D modwt", "2^$p", () -> modwt(x, wtm, L), () -> modwt(xg, wtm, L))
end

print_header("1D IMODWT")
for p in (14, 16, 18)
    n = 2^p
    x = randn(Float32, n)
    xg = to_gpu(GPU_BACKEND, x)
    L = min(6, maxmodwttransformlevels(n))
    mw_cpu = modwt(x, wtm, L)
    mw_gpu = modwt(xg, wtm, L)
    benchmark_case("1D imodwt", "2^$p", () -> imodwt(mw_cpu, wtm), () -> imodwt(mw_gpu, wtm))
end

print_header("1D MODWT Step")
let
    g, h = WT.makereverseqmfpair(wtm)
    g ./= sqrt(2)
    h ./= sqrt(2)
    j = 2
    for p in (14, 16, 18)
        n = 2^p
        x = randn(Float32, n)
        xg = to_gpu(GPU_BACKEND, x)
        benchmark_case(
            "1D modwt step",
            "2^$p",
            () -> Wavelets.Transforms.modwt_step(x, j, h, g),
            () -> Wavelets.Transforms.modwt_step(xg, j, h, g)
        )
    end
end

print_header("1D IMODWT Step")
let
    g, h = WT.makereverseqmfpair(wtm)
    g ./= sqrt(2)
    h ./= sqrt(2)
    j = 2
    for p in (14, 16, 18)
        n = 2^p
        x = randn(Float32, n)
        xg = to_gpu(GPU_BACKEND, x)
        v_cpu, w_cpu = Wavelets.Transforms.modwt_step(x, j, h, g)
        v_gpu, w_gpu = Wavelets.Transforms.modwt_step(xg, j, h, g)
        benchmark_case(
            "1D imodwt step",
            "2^$p",
            () -> Wavelets.Transforms.imodwt_step(v_cpu, w_cpu, j, h, g),
            () -> Wavelets.Transforms.imodwt_step(v_gpu, w_gpu, j, h, g)
        )
    end
end

print_header("2D Lifting DWT")
wtl2 = wavelet(WT.db2, WT.Lifting)
for n in (128, 1024, 4096)
    x = randn(Float32, n, n)
    xg = to_gpu(GPU_BACKEND, x)
    L = min(4, maxtransformlevels(n))
    yg = dwt(xg, wtl2, L)
    yc = dwt(x, wtl2, L)
    benchmark_case("2D lifting dwt", "$(n)x$(n)", () -> dwt(x, wtl2, L), () -> dwt(xg, wtl2, L))
end

print_header("2D Lifting IDWT")
wtl2 = wavelet(WT.db2, WT.Lifting)
for n in (128, 1024, 4096)
    x = randn(Float32, n, n)
    xg = to_gpu(GPU_BACKEND, x)
    L = min(4, maxtransformlevels(n))
    yg = dwt(xg, wtl2, L)
    yc = dwt(x, wtl2, L)
    benchmark_case("2D lifting idwt", "$(n)x$(n)", () -> idwt(yc, wtl2, L), () -> idwt(yg, wtl2, L))
end

print_header("3D Lifting DWT")
wtl3 = wavelet(WT.haar, WT.Lifting)
for n in (32, 128, 512)
    x = randn(Float32, n, n, n)
    xg = to_gpu(GPU_BACKEND, x)
    L = min(3, maxtransformlevels(n))
    yg = dwt(xg, wtl3, L)
    yc = dwt(x, wtl3, L)
    benchmark_case("3D lifting dwt", "$(n)^3", () -> dwt(x, wtl3, L), () -> dwt(xg, wtl3, L))
end

print_header("3D Lifting IDWT")
wtl3 = wavelet(WT.haar, WT.Lifting)
for n in (32, 128, 512)
    x = randn(Float32, n, n, n)
    xg = to_gpu(GPU_BACKEND, x)
    L = min(3, maxtransformlevels(n))
    yg = dwt(xg, wtl3, L)
    yc = dwt(x, wtl3, L)
    benchmark_case("3D lifting idwt", "$(n)^3", () -> idwt(yc, wtl3, L), () -> idwt(yg, wtl3, L))
end

println()
println(repeat("=", 78))
println("Summary")
println(repeat("=", 78))
for (name, vals) in sort(collect(results), by=first)
    @printf("%-18s\t avg: %6.2fx, min: %6.2fx, max: %6.2fx\n", name, mean(vals), minimum(vals), maximum(vals))
end