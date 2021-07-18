# FHist.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://moelf.github.io/FHist.jl/dev/)
[![Build Status](https://github.com/Moelf/FHist.jl/workflows/CI/badge.svg)](https://github.com/Moelf/FHist.jl/actions)
[![Codecov](https://codecov.io/gh/Moelf/FHist.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/Moelf/FHist.jl)

Fast, error-aware, and thread-safe 1&2D histograms that are also compatible with `StatsBase.Histogram`

## Quick Start
```julia
julia> a = rand(1000);


julia> h1 = Hist1D(a)
   [-3.0, -2.5) ┤ 3                                        
   [-2.5, -2.0) ┤▇▇▇ 18                                    
   [-2.0, -1.5) ┤▇▇▇▇▇▇ 36                                 
   [-1.5, -1.0) ┤▇▇▇▇▇▇▇▇▇▇▇▇ 74                           
   [-1.0, -0.5) ┤▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇ 159            
   [-0.5,  0.0) ┤▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇ 193       
   [ 0.0,  0.5) ┤▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇ 217   
   [ 0.5,  1.0) ┤▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇ 151              
   [ 1.0,  1.5) ┤▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇ 101                      
   [ 1.5,  2.0) ┤▇▇▇▇▇▇ 35                                 
   [ 2.0,  2.5) ┤▇ 8                                       
   [ 2.5,  3.0) ┤▇ 5                                       
                └                                        ┘ 
                                Frequency
edges: -3.0:0.5:3.0
bin counts: [3, 18, 36, 74, 159, 193, 217, 151, 101, 35, 8, 5]
errors: 
  up  : [1.73, 4.24, 6.0, 8.6, 12.6, 13.9, 14.7, 12.3, 10.0, 5.92, 2.83, 2.24]
  down: [1.73, 4.24, 6.0, 8.6, 12.6, 13.9, 14.7, 12.3, 10.0, 5.92, 2.83, 2.24]
error_mode: sqrt

julia> h2 = Hist1D(Int; bins=0:0.1:1);

julia> Threads.@threads for i in a
           push!(h2, i)
       end

julia> update_error!(h2);

julia> h1 == h2
true
```
