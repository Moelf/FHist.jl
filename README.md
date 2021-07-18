# FHist.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://moelf.github.io/FHist.jl/dev/)
[![Build Status](https://github.com/Moelf/FHist.jl/workflows/CI/badge.svg)](https://github.com/Moelf/FHist.jl/actions)
[![Codecov](https://codecov.io/gh/Moelf/FHist.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/Moelf/FHist.jl)

Fast, error-aware, and thread-safe 1&2D histograms that are also compatible with `StatsBase.Histogram`

## Quick Start
```julia
julia> a = rand(1000);

julia> h1 = Hist1D(a)
StatsBase.Histogram{Int64, 1, Tuple{StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}}}}
edges:
  0.0:0.1:1.0
weights: [81, 106, 108, 99, 88, 113, 105, 96, 91, 113]
...
...

julia> h2 = Hist1D(Int; bins=0:0.1:1);

julia> Threads.@threads for i in a
           push!(h2, i)
       end

julia> update_error!(h2);

julia> h1 == h2
true
```
