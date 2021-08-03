# FHist.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://moelf.github.io/FHist.jl/dev/)
[![Build Status](https://github.com/Moelf/FHist.jl/workflows/CI/badge.svg)](https://github.com/Moelf/FHist.jl/actions)
[![Codecov](https://codecov.io/gh/Moelf/FHist.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/Moelf/FHist.jl)

Fast, error-aware, and thread-safe 1&2D histograms that are also compatible with `StatsBase.Histogram`

## Quick Start
```julia
julia> a = randn(1000);


julia> h1 = Hist1D(a);

julia> h2 = Hist1D(Int; bins=-3:0.5:3)

julia> Threads.@threads for i in a
           push!(h2, i)
       end

julia> h1 == h2
true
```
