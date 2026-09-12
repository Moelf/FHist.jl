## Basic attributes

- [`binedges`](@ref): Get the bin edges of the histogram, for 1D histogram, it returns just a vector. For others, it returns a tuple of vectors. If you need a tuple of vectors, use h.binedges at your own risk.
- [`bincounts`](@ref): Get the bin counts (weights) of the histogram.
- [`sumw2`](@ref): Get the sum of weights squared of the histogram, it has the same shape as `bincounts(h)`.
- [`nentries`](@ref): Get the number of times a histogram is filled (`push!`ed)

## Derived attributes
- [`bincenters`](@ref): Get the bin centers of the histogram, for 1D histogram, it returns just a vector. For others, it returns a tuple of vectors.
- [`binerrors`](@ref): Get the error of each bin of the histogram. By default it calls `sqrt()` on
  each entry of `sumw2` as an approximation.
- `mean`, `std`, `median`, `quantile`, weighted by histogram `bincounts()`.
```@docs
integral
```

## Manipulating histogram

This section includes adding data to histogram, and other operations that return a histogram (may
with reduced dimensionality)
```@docs
atomic_push!
cumulative
rebin
restrict
profile
project
normalize
```

All of the above work for `Hist1D`, `Hist2D` and `Hist3D` alike (except `cumulative`, `profile`
and `transpose`); `rebin` and `restrict` take one argument (pair of arguments) per axis. See
also `append!` (bulk `push!`, thread-safe) and `empty!` (reset counts, `sumw2` and `nentries`)
in the [API reference](api.md).

## GPU

Histograms can be filled from GPU arrays (`Hist1D(cu_array; binedges = ...)`), and bin counts can
be computed on the device with [`gpu_bincounts`](@ref) / [`gpu_bincounts!`](@ref), see
[GPU histogramming](@ref).
