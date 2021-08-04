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
