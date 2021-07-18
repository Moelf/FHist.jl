# FHist.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://moelf.github.io/FHist.jl/dev/)
[![Build Status](https://github.com/Moelf/FHist.jl/workflows/CI/badge.svg)](https://github.com/Moelf/FHist.jl/actions)
[![Codecov](https://codecov.io/gh/Moelf/FHist.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/Moelf/FHist.jl)

Fast, error-aware, and thread-safe 1&2D histograms that are also compatible with `StatsBase.Histogram`

## Quick Start
```julia
julia> a = randn(1000);


julia> h1 = Hist1D(a)
   [-3.0, -2.5) ┤▇ 8                             
   [-2.5, -2.0) ┤▇▇ 16                           
   [-2.0, -1.5) ┤▇▇▇▇▇▇ 51                       
   [-1.5, -1.0) ┤▇▇▇▇▇▇▇▇▇▇▇ 95                  
   [-1.0, -0.5) ┤▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇ 152          
   [-0.5,  0.0) ┤▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇ 207   
   [ 0.0,  0.5) ┤▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇ 194     
   [ 0.5,  1.0) ┤▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇ 133            
   [ 1.0,  1.5) ┤▇▇▇▇▇▇▇▇▇▇ 84                   
   [ 1.5,  2.0) ┤▇▇▇▇▇ 42                        
   [ 2.0,  2.5) ┤▇ 11                            
   [ 2.5,  3.0) ┤▇ 7                             
                └                              ┘ 
                           Frequency
edges: -3.0:0.5:3.0
bin counts: [8, 16, 51, 95, 152, 207, 194, 133, 84, 42, 11, 7]
errors: 
  up  : [2.83, 4.0, 7.14, 9.75, 12.3, 14.4, 13.9, 11.5, 9.17, 6.48, 3.32, 2.65]
  down: [2.83, 4.0, 7.14, 9.75, 12.3, 14.4, 13.9, 11.5, 9.17, 6.48, 3.32, 2.65]
error_mode: sqrt

julia> h2 = Hist1D(Int; bins=-3:0.5:3)

julia> Threads.@threads for i in a
           push!(h2, i)
       end

julia> update_error!(h2);

julia> h1 == h2
true
```
