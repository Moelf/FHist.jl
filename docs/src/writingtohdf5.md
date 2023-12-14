# Write out a histogram to HDF5

`FHist.jl` provides an `HDF5.jl` extension to dump and read histograms
(`Hist1D`, `Hist2D` and `Hist3D`) to HDF5 files. Once the `HDF5` package is
loaded, the corresponding methods for [`h5dumphist`](@ref) and
[`h5readhist`](@ref) will become available:


```@example 1
using FHist
using HDF5

h = Hist1D(rand(1000), -3:0.3:3)
h5dumphist("foo.h5", "some/path/to/myhist", h)
```

The histogram is now saved in `foo.h5` in the `some/path/to/myhist` group
including all the meta information needed to be able to recreate it.

```@example 1
h5readhist("foo.h5", "some/path/to/myhist", Hist1D)
```

