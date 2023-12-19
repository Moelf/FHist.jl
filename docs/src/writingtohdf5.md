# Write out a histogram to HDF5

`FHist.jl` provides an `HDF5.jl` extension to write and read histograms
(`Hist1D`, `Hist2D` and `Hist3D`) to HDF5 files. Once the `HDF5` package is
loaded, the corresponding methods for `h5writehist` and
`h5readhist` will become available.

The HDF5 format specification used in `FHist.jl` for histograms is
language/library agnostic and easily portable to other libraries. Counts
(weights, including `sumw2`), edges and other histogram parameters are stored as
HDF5 datasets or attributes.

Let's create a `Hist1D`:

```@example 1
using FHist
using HDF5

h = Hist1D(randn(10_000), -3:0.1:3)
```

Now write it to an HDF5 file:

```@example 1
h5writehist("foo.h5", "some/path/to/myhist", h)
```

The histogram is now saved in `foo.h5` in the `some/path/to/myhist` group
including all the meta information needed to be able to recreate it. This is how
it looks like when opening it with `HDF5.jl` (any other HDF5 library will look
similarly):

```
f = HDF5.File: (read-only) foo.h5
🗂️ HDF5.File: (read-only) foo.h5
└─ 📂 some
   └─ 📂 path
      └─ 📂 to
         └─ 📂 myhist
            ├─ 🏷️ _producer
            ├─ 🏷️ closed
            ├─ 🏷️ isdensity
            ├─ 🏷️ nentries
            ├─ 🏷️ overflow
            ├─ 🔢 edges_1
            ├─ 🔢 sumw2
            └─ 🔢 weights
```

Here is how to read it back to an actual `FHist.Hist1D` instance:

```@example 1
h5readhist("foo.h5", "some/path/to/myhist")
```

