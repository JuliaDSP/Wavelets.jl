using BenchmarkTools
using Random
using Wavelets

Random.seed!(42)

const SUITE = BenchmarkGroup()

# ---------------------------------------------------------------------------
# Wavelet instances
# ---------------------------------------------------------------------------
const WTF  = wavelet(WT.db4,    WT.Filter)   # filter
const WTFL = wavelet(WT.haar,   WT.Filter)   # filter, haar (smallest kernel)
const WTL  = wavelet(WT.haar,   WT.Lifting)  # lifting, haar
const WTL2 = wavelet(WT.db2,    WT.Lifting)  # lifting, db2

# ---------------------------------------------------------------------------
# 1-D transforms
# ---------------------------------------------------------------------------
SUITE["1d"] = BenchmarkGroup(["1d"])

let g = SUITE["1d"]
    for (wname, wt) in (("filter_db4", WTF), ("lifting_haar", WTL), ("lifting_db2", WTL2))
        g[wname] = BenchmarkGroup()
        for p in (10, 15, 20)
            n  = 2^p
            x  = randn(Float32, n)
            L  = min(6, maxtransformlevels(n))
            y  = dwt(x, wt, L)

            g[wname]["dwt/2^$p"]   = @benchmarkable dwt($x,       $wt, $L) evals=1
            g[wname]["idwt/2^$p"]  = @benchmarkable idwt($y,      $wt, $L) evals=1
            g[wname]["wpt/2^$p"]   = @benchmarkable wpt($x,       $wt)     evals=1
            yw = wpt(x, wt)
            g[wname]["iwpt/2^$p"]  = @benchmarkable iwpt($yw,     $wt)     evals=1
        end
    end

    # MODWT — filter only, 1-D only
    # modwt returns an (N, L+1) matrix; imodwt takes that matrix (no L arg)
    let wname = "filter_haar_modwt"
        g[wname] = BenchmarkGroup()
        wt = WTFL
        for p in (10, 15, 20)
            n = 2^p
            x = randn(Float32, n)
            L = min(6, maxmodwttransformlevels(n))
            y = modwt(x, wt, L)
            g[wname]["modwt/2^$p"]  = @benchmarkable modwt($x, $wt, $L) evals=1
            g[wname]["imodwt/2^$p"] = @benchmarkable imodwt($y, $wt)    evals=1
        end
    end
end

# ---------------------------------------------------------------------------
# 2-D transforms
# ---------------------------------------------------------------------------
SUITE["2d"] = BenchmarkGroup(["2d"])

let g = SUITE["2d"]
    for (wname, wt) in (("filter_db4", WTF), ("lifting_haar", WTL), ("lifting_db2", WTL2))
        g[wname] = BenchmarkGroup()
        for n in (128, 512, 2048)
            x = randn(Float32, n, n)
            L = min(4, maxtransformlevels(n))
            y = dwt(x, wt, L)
            g[wname]["dwt/$(n)x$(n)"]  = @benchmarkable dwt($x,  $wt, $L) evals=1
            g[wname]["idwt/$(n)x$(n)"] = @benchmarkable idwt($y, $wt, $L) evals=1
        end
    end
end

# ---------------------------------------------------------------------------
# 3-D transforms
# ---------------------------------------------------------------------------
SUITE["3d"] = BenchmarkGroup(["3d"])

let g = SUITE["3d"]
    for (wname, wt) in (("filter_db4", WTF), ("lifting_haar", WTL), ("lifting_db2", WTL2))
        g[wname] = BenchmarkGroup()
        for n in (32, 128, 256)
            x = randn(Float32, n, n, n)
            L = min(3, maxtransformlevels(n))
            y = dwt(x, wt, L)
            g[wname]["dwt/$(n)^3"]  = @benchmarkable dwt($x,  $wt, $L) evals=1
            g[wname]["idwt/$(n)^3"] = @benchmarkable idwt($y, $wt, $L) evals=1
        end
    end
end
