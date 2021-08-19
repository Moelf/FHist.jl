var documenterSearchIndex = {"docs":
[{"location":"api/","page":"APIs","title":"APIs","text":"Modules = [FHist]","category":"page"},{"location":"api/#FHist.Hist1D-Tuple{AbstractVector{T} where T, AbstractRange}","page":"APIs","title":"FHist.Hist1D","text":"Hist1D(array, edges::AbstractRange)\nHist1D(array, edges::AbstractVector)\n\nCreate a Hist1D with given bin edges and vlaues from array. Weight for each value is assumed to be 1.\n\n\n\n\n\n","category":"method"},{"location":"api/#FHist.Hist1D-Tuple{Any, StatsBase.AbstractWeights, AbstractRange}","page":"APIs","title":"FHist.Hist1D","text":"Hist1D(array, wgts::AbstractWeights, edges::AbstractRange)\nHist1D(array, wgts::AbstractWeights, edges::AbstractVector)\n\nCreate a Hist1D with given bin edges and vlaues from array. wgts should have the same size as array.\n\n\n\n\n\n","category":"method"},{"location":"api/#FHist.Hist1D-Union{Tuple{AbstractVector{T}}, Tuple{T}} where T","page":"APIs","title":"FHist.Hist1D","text":"Hist1D(A::AbstractVector{T}; nbins::Integer=StatsBase.sturges(length(A))) where T\nHist1D(A::AbstractVector{T}, wgts::AbstractWeights; nbins::Integer=StatsBase.sturges(length(A))) where T\n\nAutomatically determine number of bins based on Sturges algo.\n\n\n\n\n\n","category":"method"},{"location":"api/#FHist.Hist1D-Union{Tuple{}, Tuple{Type{T}}, Tuple{T}} where T","page":"APIs","title":"FHist.Hist1D","text":"Hist1D(elT::Type{T}=Float64; binedges) where {T}\n\nInitialize an empty histogram with bin content typed as T and bin edges. To be used with push!\n\n\n\n\n\n","category":"method"},{"location":"api/#Base.empty!-Union{Tuple{Hist1D{T, E}}, Tuple{E}, Tuple{T}} where {T, E}","page":"APIs","title":"Base.empty!","text":"empty!(h::Hist1D)\n\nResets a histogram's bin counts and sumw2.\n\n\n\n\n\n","category":"method"},{"location":"api/#Base.push!-Union{Tuple{E}, Tuple{T}, Tuple{Hist1D{T, E}, Real}, Tuple{Hist1D{T, E}, Real, Real}} where {T, E}","page":"APIs","title":"Base.push!","text":"unsafe_push!(h::Hist1D, val::Real, wgt::Real=1)\npush!(h::Hist1D, val::Real, wgt::Real=1)\n\nAdding one value at a time into histogram.  sumw2 (sum of weights^2) accumulates wgt^2 with a default weight of 1. unsafe_push! is a faster version of push! that is not thread-safe.\n\nN.B. To append multiple values at once, use broadcasting via push!.(h, [-3.0, -2.9, -2.8]) or push!.(h, [-3.0, -2.9, -2.8], 2.0)\n\n\n\n\n\n","category":"method"},{"location":"api/#FHist.bincenters-Tuple{Hist1D}","page":"APIs","title":"FHist.bincenters","text":"bincenters(h::Hist1D)\n\nGet the bin centers of a histogram.\n\n\n\n\n\n","category":"method"},{"location":"api/#FHist.bincounts-Tuple{Hist1D}","page":"APIs","title":"FHist.bincounts","text":"bincounts(h::Hist1D)\n\nGet the bin counts of a histogram.\n\n\n\n\n\n","category":"method"},{"location":"api/#FHist.binedges-Tuple{Hist1D}","page":"APIs","title":"FHist.binedges","text":"binedges(h::Hist1D)\n\nGet the bin edges of a histogram.\n\n\n\n\n\n","category":"method"},{"location":"api/#FHist.cumulative-Tuple{Hist1D}","page":"APIs","title":"FHist.cumulative","text":"cumulative(h::Hist1D; forward=true)\n\nCreate a cumulative histogram. If forward, start summing from the left.\n\n\n\n\n\n","category":"method"},{"location":"api/#FHist.hists_to_bars-Tuple{Any}","page":"APIs","title":"FHist.hists_to_bars","text":"hists_to_bars(hist1ds)\n\nGiven a vector of Hist1D, return edges (xs), heights (ys), and grps (for grouping) that is useful for plotting stacked histogram.\n\n\n\n\n\n","category":"method"},{"location":"api/#FHist.integral-Tuple{Hist1D}","page":"APIs","title":"FHist.integral","text":"integral(h::Hist1D)\n\nGet the integral a histogram.\n\n\n\n\n\n","category":"method"},{"location":"api/#FHist.lookup-Tuple{Hist1D, Any}","page":"APIs","title":"FHist.lookup","text":"function lookup(h::Hist1D, v)\n\nFor given x-axis valuev, find the corresponding bin and return the bin content. If a value is out of the histogram range, return missing.\n\n\n\n\n\n","category":"method"},{"location":"api/#FHist.nbins-Tuple{Hist1D}","page":"APIs","title":"FHist.nbins","text":"nbins(h::Hist1D)\n\nGet the number of bins of a histogram.\n\n\n\n\n\n","category":"method"},{"location":"api/#FHist.rebin","page":"APIs","title":"FHist.rebin","text":"rebin(h::Hist1D, n::Int=1)\nrebin(n::Int) = h::Hist1D -> rebin(h, n)\n\nMerges n consecutive bins into one. The returned histogram will have nbins(h)/n bins.\n\n\n\n\n\n","category":"function"},{"location":"api/#FHist.sample-Tuple{Hist1D}","page":"APIs","title":"FHist.sample","text":"sample(h::Hist1D, n::Int=1)\n\nSample a histogram's with weights equal to bin count, n times. The sampled values are the bins' lower edges.\n\n\n\n\n\n","category":"method"},{"location":"api/#LinearAlgebra.normalize-Tuple{Hist1D}","page":"APIs","title":"LinearAlgebra.normalize","text":"normalize(h::Hist1D)\n\nCreate a normalized histogram via division by integral(h).\n\n\n\n\n\n","category":"method"},{"location":"api/#Statistics.mean-Tuple{Hist1D}","page":"APIs","title":"Statistics.mean","text":"Statistics.mean(h::Hist1D)\nStatistics.std(h::Hist1D)\nStatistics.median(h::Hist1D)\nStatistics.quantile(h::Hist1D, p)\n\nCompute statistical quantities based on the bin centers weighted by the bin counts.\n\n\n\n\n\n","category":"method"},{"location":"#FHist.jl","page":"Introduction","title":"FHist.jl","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"(Image: Dev) (Image: Build Status) (Image: Codecov)","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Fast, error-aware, and thread-safe 1&2D histograms that are also compatible with StatsBase.Histogram","category":"page"},{"location":"#Quick-Start","page":"Introduction","title":"Quick Start","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"julia> using FHist\n\njulia> a = randn(1000);\n\njulia> h1 = Hist1D(a)\n                ┌                              ┐\n   [-4.0, -3.0) ┤ 1\n   [-3.0, -2.0) ┤▇▇ 25\n   [-2.0, -1.0) ┤▇▇▇▇▇▇▇▇▇▇ 139\n   [-1.0,  0.0) ┤▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇ 341\n   [ 0.0,  1.0) ┤▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇ 347\n   [ 1.0,  2.0) ┤▇▇▇▇▇▇▇▇▇ 123\n   [ 2.0,  3.0) ┤▇▇ 23\n   [ 3.0,  4.0) ┤ 1\n                └                              ┘\nedges: -4.0:1.0:4.0\nbin counts: [1, 25, 139, 341, 347, 123, 23, 1]\ntotal count: 1000\n\njulia> h2 = Hist1D(Int; bins=-3:0.5:3)\n\njulia> Threads.@threads for i in a\n           push!(h2, i)\n       end\n\njulia> h1 == h2\ntrue","category":"page"},{"location":"#Speed","page":"Introduction","title":"Speed","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"Single-threaded filling happens at ~250 MHz","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"julia> a = randn(10^6);\n\njulia> @benchmark Hist1D(a, -3:0.01:3)\n Range (min … max):  4.040 ms …   5.571 ms  ┊ GC (min … max): 0.00% … 0.00%\n Time  (median):     4.393 ms               ┊ GC (median):    0.00%\n Time  (mean ± σ):   4.460 ms ± 198.070 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%","category":"page"},{"location":"#Features","page":"Introduction","title":"Features","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"julia> using FHist, Statistics, Random\n\njulia> Random.seed!(42)\nMersenneTwister(42)\n\njulia> h1 = Hist1D(randn(10^4).+2, -5:0.5:5);\n\njulia> h1 = (h1 + h1*2)\n                ┌                              ┐\n   [-5.0, -4.5) ┤ 0\n   [-4.5, -4.0) ┤ 0\n   [-4.0, -3.5) ┤ 0\n   [-3.5, -3.0) ┤ 0\n   [-3.0, -2.5) ┤ 0\n   [-2.5, -2.0) ┤ 3\n   [-2.0, -1.5) ┤ 3\n   [-1.5, -1.0) ┤ 39\n   [-1.0, -0.5) ┤▇ 153\n   [-0.5,  0.0) ┤▇▇ 540\n   [ 0.0,  0.5) ┤▇▇▇▇▇▇ 1356\n   [ 0.5,  1.0) ┤▇▇▇▇▇▇▇▇▇▇▇ 2640\n   [ 1.0,  1.5) ┤▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇ 4551\n   [ 1.5,  2.0) ┤▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇ 5757\n   [ 2.0,  2.5) ┤▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇ 5727\n   [ 2.5,  3.0) ┤▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇ 4494\n   [ 3.0,  3.5) ┤▇▇▇▇▇▇▇▇▇▇▇ 2652\n   [ 3.5,  4.0) ┤▇▇▇▇▇▇ 1326\n   [ 4.0,  4.5) ┤▇▇ 546\n   [ 4.5,  5.0) ┤▇ 168\n                └                              ┘\nedges: -5.0:0.5:5.0\nbin counts: [0, 0, 0, 0, 0, 3, 3, 39, 153, 540, 1356, 2640, 4551, 5757, 5727, 4494, 2652, 1326, 546, 168]\ntotal count: 29955\n\njulia> bincenters(h1)\n-4.75:0.5:4.75\n\njulia> println(bincounts(h1))\n[0, 0, 0, 0, 0, 3, 3, 39, 153, 540, 1356, 2640, 4551, 5757, 5727, 4494, 2652, 1326, 546, 168]\n\njulia> println(h1.sumw2);\n[0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 5.0, 65.0, 255.0, 900.0, 2260.0, 4400.0, 7585.0, 9595.0, 9545.0, 7490.0, 4420.0, 2210.0, 910.0, 280.0]\n\njulia> nbins(h1), integral(h1)\n(20, 29955)\n\njulia> mean(h1), std(h1)\n(1.993865798698047, 1.0126978885814524)\n\njulia> median(h1), quantile(h1, 0.5)\n(1.7445284002084418, 1.7445284002084418)\n\njulia> lookup.(h1, [-1.5,0,2.5]) # find bin counts for given x-axis values\n3-element Vector{Int64}:\n   39\n 1356\n 4494\n\njulia> Hist1D(sample(h1; n=100))\n                ┌                              ┐\n   [-1.0,  0.0) ┤▇▇▇ 4\n   [ 0.0,  1.0) ┤▇▇▇▇▇▇▇▇▇▇ 15\n   [ 1.0,  2.0) ┤▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇ 41\n   [ 2.0,  3.0) ┤▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇ 28\n   [ 3.0,  4.0) ┤▇▇▇▇▇▇ 10\n   [ 4.0,  5.0) ┤▇ 2\n                └                              ┘\nedges: -1.0:1.0:5.0\nbin counts: [4, 15, 41, 28, 10, 2]\ntotal count: 100\n\njulia> h1 |> normalize |> integral\n1.0","category":"page"}]
}
