# FHist.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://moelf.github.io/FHist.jl/dev/)
[![Build Status](https://github.com/Moelf/FHist.jl/workflows/CI/badge.svg)](https://github.com/Moelf/FHist.jl/actions)
[![Codecov](https://codecov.io/gh/Moelf/FHist.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/Moelf/FHist.jl)

Fast, error-aware, and thread-safe 1D/2D/3D histograms that are also compatible with `StatsBase.Histogram`

## Changelog
- 0.10
  - due to huge latency issue, UnicodePlots.jl (for text based terminal display) has been dropped.

## Quick Start
```julia
julia> using FHist

julia> a = randn(1000);

julia> h1 = Hist1D(a, -4:1.0:4)
edges: -4.0:1.0:4.0
bin counts: [0, 23, 123, 377, 320, 134, 21, 2]
total count: 1000

julia> h2 = Hist1D(Int; bins=-4.0:1.0:4.0);

julia> push!.(h2, a);

julia> Threads.@threads for v in a
           atomic_push!(h2, v) # thread-safe
       end

julia> h1+h1 == h2
true
```

Additionally, one can specify `overflow=true` when creating a histogram to clamp out-of-bounds values into 
the edge bins.
```julia
julia> Hist1D(rand(1000),0:0.2:0.9; overflow=true)
edges: 0.0:0.2:0.8
bin counts: [218, 183, 198, 401]
total count: 1000
```

## Speed

Single-threaded filling happens at almost 400 MHz
```julia
julia> a = randn(10^6);

julia> @benchmark Hist1D(a, -3:0.01:3)
 Range (min … max):  2.591 ms …   4.746 ms  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     2.802 ms               ┊ GC (median):    0.00%
 Time  (mean ± σ):   2.864 ms ± 195.931 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%
```

## Features

### 1D

```julia
julia> using FHist, Statistics, Random

julia> Random.seed!(42)
MersenneTwister(42)

julia> h1 = Hist1D(randn(10^4).+2, -5:0.5:5);

julia> h1 = (h1 + h1*2)
edges: -5.0:0.5:5.0
bin counts: [0, 0, 0, 0, 0, 3, 3, 39, 153, 540, 1356, 2640, 4551, 5757, 5727, 4494, 2652, 1326, 546, 168]
total count: 29955

julia> bincenters(h1)
-4.75:0.5:4.75

julia> println(bincounts(h1))
[0, 0, 0, 0, 0, 3, 3, 39, 153, 540, 1356, 2640, 4551, 5757, 5727, 4494, 2652, 1326, 546, 168]

julia> println(h1.sumw2);
[0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 5.0, 65.0, 255.0, 900.0, 2260.0, 4400.0, 7585.0, 9595.0, 9545.0, 7490.0, 4420.0, 2210.0, 910.0, 280.0]

julia> nbins(h1), integral(h1)
(20, 29955)

julia> mean(h1), std(h1)
(1.993865798698047, 1.0126978885814524)

julia> median(h1), quantile(h1, 0.5)
(1.7445284002084418, 1.7445284002084418)

julia> lookup.(h1, [-1.5,0,2.5]) # find bin counts for given x-axis values
3-element Vector{Int64}:
   39
 1356
 4494

julia> Hist1D(sample(h1; n=100))
edges: -1.0:1.0:5.0
bin counts: [4, 15, 41, 28, 10, 2]
total count: 100

julia> h1 |> normalize |> integral
1.0
```

### 2D

```julia
julia> h2 = Hist2D((randn(10^4),randn(10^4)), (-5:5,-5:5))

edges: (-5:5, -5:5)
bin counts: [0 0 … 0 0; 0 0 … 0 0; … ; 0 0 … 0 0; 0 0 … 0 0]
total count: 10000

julia> nbins(h2)
(10, 10)

julia> project(h2, :x)
edges: -5:5
bin counts: [0, 14, 228, 1345, 3399, 3462, 1300, 237, 15, 0]
total count: 10000

julia> h2 |> profile(:x) # or pipe

edges: -5:5
bin counts: [0.0, -0.21428571428571427, 0.0, 0.03382899628252788, 0.0025007355104442485, 0.012709416522241479, 0.018461538461538463, 0.035864978902953586, 0.1, 0.0]
total count: -0.010920048606008592

julia> h2 |> rebin(10,5)

edges: (-5.0:10.0:5.0, -5.0:5.0:5.0)
bin counts: [4953 5047]
total count: 10000
```

### 3D

```julia
julia> h3 = Hist3D((randn(10^4),randn(10^4),randn(10^4)), (-5:5,-5:5,-5:5))
Hist3D{Int64}, edges=(-5:5, -5:5, -5:5), integral=10000

julia> h3 |> project(:x) |> project(:x) |> std
1.051590599996025
```
