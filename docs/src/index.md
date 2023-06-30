# FHist Package

The [*FHist*](https://github.com/Moelf/FHist.jl) package is a fast histogram library for 1D/2D/3D histograms that
integrates with `StatsBase.Histogram`. Among the usesful utilties are:

* Creation of histograms from data
* Creation of empty histograms that are then filled with data with `push!`
* Histogram arithmetic, such as `+`, `-`, `*` and `/` with scalars or other histograms
* Random sampling from histograms
* Projections and profiles of 2D/3D histograms
* Statistics of histograms such as `mean`, `median`, `std` and `quantile`
* Thread-safe operations with `atomic_push!`
* Plot receipes to work with `Makie.jl` and `Plots.jl`

# Hist1D Data Type

The `Hist1D` data type has several different methods for creating a one-dimensional histogram. The 
easiest ways are `Hist1D(array, binedges)` and `Hist1D(type, binedges)`. The first method takes
an array of data and creates the histogram whereas the second method creates an empty histogram and
can then have values `push!`-ed onto it.

You can then perform various operations on these histograms:

* `sample`: return random sample from histogram
* `integral`: return the integral of the histogram
* `lookup`: find the closest bin and return the bin contents
* `normalize`: normalize the histogram so integral value is one
* `cumulative`: return a cumulative histogram
* `rebin`: merge `n` consecutive bins into one for all bins
* `restrict`: return a histogram restricted to certain x-values
* `empty!`: clear bin contents and sumw2

To access various aspects of the histograms:

* `bincenters(h)`: return midpoint of bin edges
* `binedges(h)`: return bin edges
* `bincounts(h)`: return the counts of each bin
* `binerrors(h)`: defaults to sqrt of sum of weights squared, can pass user-defined function also
* `nbins(h)`: return number of bins of histogram

# Hist2D Data Type

The `Hist2D` data type is very similar to `Hist1D` where 2-tuples are passed as arguments instead of 
single vectors for `array` and `binedges`. In addition to the same functions above (expect for `sample`), `Hist2D` has:

* `project`: return a projected `Hist1D` for the given axis
* `profile`: return a `Hist1D` profile of the given axis
* `transpose`: reverse the `x` and `y` axes

# Hist3D Data Type

The `Hist3D` data type provides a third `z` dimension and has the same functions available to it as
the `Hist2D` type.