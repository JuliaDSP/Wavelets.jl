# GPU tests for Wavelets.jl extension.
# The test runner activates a GPUEnv overlay before this file is included.

using GPUEnv
using Test
using Wavelets
using JLArrays

gpu_allclose(x, y; atol) = maximum(abs.(x .- Array(y))) < atol

@testset "GPU Wavelets" begin

@testset "$(backend.name)" for backend in gpu_backends()
    @testset "1D Filter DWT" begin
        for wname in (WT.haar, WT.db2, WT.db4)
            wt = wavelet(wname)
            for n in (16, 64)
                x_cpu = randn(Float32, n)
                x_gpu = to_gpu(backend, x_cpu)
                Lmax = min(3, maxtransformlevels(n))
                for L in 1:Lmax
                    y_cpu = dwt(x_cpu, wt, L)
                    y_gpu = dwt(x_gpu, wt, L)
                    @test maximum(abs.(y_cpu .- Array(y_gpu))) < 1e-6

                    r_cpu = idwt(y_cpu, wt, L)
                    r_gpu = idwt(y_gpu, wt, L)
                    @test maximum(abs.(r_cpu .- Array(r_gpu))) < 1e-6
                end
            end
        end
    end

    @testset "1D Lifting DWT" begin
        for wname in (WT.haar, WT.db2, WT.cdf97)
            wt = wavelet(wname, WT.Lifting)
            for n in (16, 64)
                x_cpu = randn(Float32, n)
                x_gpu = to_gpu(backend, x_cpu)
                Lmax = min(3, maxtransformlevels(n))
                for L in 1:Lmax
                    y_cpu = dwt(x_cpu, wt, L)
                    y_gpu = dwt(x_gpu, wt, L)
                    @test maximum(abs.(y_cpu .- Array(y_gpu))) < 1e-5

                    r_cpu = idwt(y_cpu, wt, L)
                    r_gpu = idwt(y_gpu, wt, L)
                    @test maximum(abs.(r_cpu .- Array(r_gpu))) < 1e-5
                end
            end
        end
    end

    @testset "2D Filter DWT" begin
        wt = wavelet(WT.haar)
        for n in (8, 16)
            x_cpu = randn(Float32, n, n)
            x_gpu = to_gpu(backend, x_cpu)
            Lmax = min(2, maxtransformlevels(n))
            for L in 1:Lmax
                y_cpu = dwt(x_cpu, wt, L)
                y_gpu = dwt(x_gpu, wt, L)
                @test maximum(abs.(y_cpu .- Array(y_gpu))) < 1e-5

                r_cpu = idwt(y_cpu, wt, L)
                r_gpu = idwt(y_gpu, wt, L)
                @test maximum(abs.(r_cpu .- Array(r_gpu))) < 1e-5
            end
        end
    end

    @testset "3D Filter DWT" begin
        wt = wavelet(WT.haar)
        n = 8
        x_cpu = randn(Float32, n, n, n)
        x_gpu = to_gpu(backend, x_cpu)
        for L in 1:2
            y_cpu = dwt(x_cpu, wt, L)
            y_gpu = dwt(x_gpu, wt, L)
            @test gpu_allclose(y_cpu, y_gpu; atol = 1e-5)

            r_cpu = idwt(y_cpu, wt, L)
            r_gpu = idwt(y_gpu, wt, L)
            @test gpu_allclose(r_cpu, r_gpu; atol = 1e-5)
        end
    end

    @testset "2D Lifting DWT" begin
        for wname in (WT.haar, WT.db2)
            wt = wavelet(wname, WT.Lifting)
            for n in (8, 16)
                x_cpu = randn(Float32, n, n)
                x_gpu = to_gpu(backend, x_cpu)
                Lmax = min(2, maxtransformlevels(n))
                for L in 1:Lmax
                    y_cpu = dwt(x_cpu, wt, L)
                    y_gpu = dwt(x_gpu, wt, L)
                    @test gpu_allclose(y_cpu, y_gpu; atol = 1e-5)

                    r_cpu = idwt(y_cpu, wt, L)
                    r_gpu = idwt(y_gpu, wt, L)
                    @test gpu_allclose(r_cpu, r_gpu; atol = 1e-5)
                end
            end
        end
    end

    @testset "WPT Filter" begin
        wt = wavelet(WT.haar)
        for n in (16, 64)
            x_cpu = randn(Float32, n)
            x_gpu = to_gpu(backend, x_cpu)

            y_cpu = wpt(x_cpu, wt)
            y_gpu = wpt(x_gpu, wt)
            @test maximum(abs.(y_cpu .- Array(y_gpu))) < 1e-5

            r_cpu = iwpt(y_cpu, wt)
            r_gpu = iwpt(y_gpu, wt)
            @test maximum(abs.(r_cpu .- Array(r_gpu))) < 1e-5
        end
    end

    @testset "WPT Lifting" begin
        for wname in (WT.haar, WT.db2, WT.cdf97)
            wt = wavelet(wname, WT.Lifting)
            for n in (16, 64)
                x_cpu = randn(Float32, n)
                x_gpu = to_gpu(backend, x_cpu)

                y_cpu = wpt(x_cpu, wt)
                y_gpu = wpt(x_gpu, wt)
                @test gpu_allclose(y_cpu, y_gpu; atol = 1e-5)

                r_cpu = iwpt(y_cpu, wt)
                r_gpu = iwpt(y_gpu, wt)
                @test gpu_allclose(r_cpu, r_gpu; atol = 1e-5)
            end
        end
    end

    @testset "MODWT" begin
        wt = wavelet(WT.db4)
        x_cpu = randn(Float32, 64)
        x_gpu = to_gpu(backend, x_cpu)
        L = min(4, maxmodwttransformlevels(length(x_cpu)))

        W_cpu = modwt(x_cpu, wt, L)
        W_gpu = modwt(x_gpu, wt, L)
        @test gpu_allclose(W_cpu, W_gpu; atol = 1e-5)

        xr_cpu = imodwt(W_cpu, wt)
        xr_gpu = imodwt(W_gpu, wt)
        @test gpu_allclose(xr_cpu, xr_gpu; atol = 1e-5)

        using Wavelets.Transforms: modwt_step, imodwt_step

        g, h = WT.makereverseqmfpair(wt, true, Float32)
        g ./= sqrt(2)
        h ./= sqrt(2)
        v1_cpu, w1_cpu = modwt_step(x_cpu, 2, h, g)
        v1_gpu, w1_gpu = modwt_step(x_gpu, 2, h, g)
        @test gpu_allclose(v1_cpu, v1_gpu; atol = 1e-5)
        @test gpu_allclose(w1_cpu, w1_gpu; atol = 1e-5)

        v0_cpu = imodwt_step(v1_cpu, w1_cpu, 2, h, g)
        v0_gpu = imodwt_step(v1_gpu, w1_gpu, 2, h, g)
        @test gpu_allclose(v0_cpu, v0_gpu; atol = 1e-5)
    end

    @testset "3D Lifting DWT" begin
        wt = wavelet(WT.haar, WT.Lifting)
        n = 8
        x_cpu = randn(Float32, n, n, n)
        x_gpu = to_gpu(backend, x_cpu)
        for L in 1:2
            y_cpu = dwt(x_cpu, wt, L)
            y_gpu = dwt(x_gpu, wt, L)
            @test gpu_allclose(y_cpu, y_gpu; atol = 1e-5)

            r_cpu = idwt(y_cpu, wt, L)
            r_gpu = idwt(y_gpu, wt, L)
            @test gpu_allclose(r_cpu, r_gpu; atol = 1e-5)
        end
    end

    @testset "Type preservation" begin
        wt = wavelet(WT.haar)
        x32_gpu = to_gpu(backend, randn(Float32, 64))
        y = dwt(x32_gpu, wt)
        @test eltype(y) == Float32
    end
end

end
