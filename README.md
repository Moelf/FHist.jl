# FHist.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://moelf.github.io/FHist.jl/dev/)
[![Build Status](https://github.com/Moelf/FHist.jl/workflows/CI/badge.svg)](https://github.com/Moelf/FHist.jl/actions)
[![Codecov](https://codecov.io/gh/Moelf/FHist.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/Moelf/FHist.jl)

Fast, error-aware, and thread-safe 1&2D histograms that are also compatible with `StatsBase.Histogram`

## Quick Start
```julia
julia> using FHist

julia> a = randn(1000);

julia> h1 = Hist1D(a)
                ┌                              ┐
   [-4.0, -3.0) ┤ 1
   [-3.0, -2.0) ┤▇▇ 25
   [-2.0, -1.0) ┤▇▇▇▇▇▇▇▇▇▇ 139
   [-1.0,  0.0) ┤▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇ 341
   [ 0.0,  1.0) ┤▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇ 347
   [ 1.0,  2.0) ┤▇▇▇▇▇▇▇▇▇ 123
   [ 2.0,  3.0) ┤▇▇ 23
   [ 3.0,  4.0) ┤ 1
                └                              ┘
edges: -4.0:1.0:4.0
bin counts: [1, 25, 139, 341, 347, 123, 23, 1]
total count: 1000

julia> h2 = Hist1D(Int; bins=-3:0.5:3)

julia> Threads.@threads for i in a
           push!(h2, i)
       end

julia> h1 == h2
true
```

## Features
```julia
julia> using FHist, Statistics, Random

julia> Random.seed!(42)
MersenneTwister(42)

julia> h1 = Hist1D(randn(10^4).+2, -5:0.5:5);

julia> h1 = (h1 + h1*2)
                ┌                              ┐
   [-5.0, -4.5) ┤ 0
   [-4.5, -4.0) ┤ 0
   [-4.0, -3.5) ┤ 0
   [-3.5, -3.0) ┤ 0
   [-3.0, -2.5) ┤ 0
   [-2.5, -2.0) ┤ 3
   [-2.0, -1.5) ┤ 3
   [-1.5, -1.0) ┤ 39
   [-1.0, -0.5) ┤▇ 153
   [-0.5,  0.0) ┤▇▇ 540
   [ 0.0,  0.5) ┤▇▇▇▇▇▇ 1356
   [ 0.5,  1.0) ┤▇▇▇▇▇▇▇▇▇▇▇ 2640
   [ 1.0,  1.5) ┤▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇ 4551
   [ 1.5,  2.0) ┤▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇ 5757
   [ 2.0,  2.5) ┤▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇ 5727
   [ 2.5,  3.0) ┤▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇ 4494
   [ 3.0,  3.5) ┤▇▇▇▇▇▇▇▇▇▇▇ 2652
   [ 3.5,  4.0) ┤▇▇▇▇▇▇ 1326
   [ 4.0,  4.5) ┤▇▇ 546
   [ 4.5,  5.0) ┤▇ 168
                └                              ┘
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

julia> lookup.(Ref(h1), [-1.5,0,2.5]) # find bin counts for given x-axis values
3-element Vector{Int64}:
   39
 1356
 4494

julia> Hist1D(sample(h1; n=100))
                ┌                              ┐
   [-1.0,  0.0) ┤▇▇▇ 4
   [ 0.0,  1.0) ┤▇▇▇▇▇▇▇▇▇▇ 15
   [ 1.0,  2.0) ┤▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇ 41
   [ 2.0,  3.0) ┤▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇ 28
   [ 3.0,  4.0) ┤▇▇▇▇▇▇ 10
   [ 4.0,  5.0) ┤▇ 2
                └                              ┘
edges: -1.0:1.0:5.0
bin counts: [4, 15, 41, 28, 10, 2]
total count: 100

julia> h1 |> normalize |> integral
1.0
```
