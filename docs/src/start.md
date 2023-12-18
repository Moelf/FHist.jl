# Quick Start

This section describes the most important surface level API: creating histogram.
```@contents
Pages = ["start.md"]
```

```@setup fh
using FHist
```

There are two possible ways to make a histogram:
- make an empty histogram to be filled later (via [`atomic_push!`](@ref))
- make a filled histogram given data

We call the first one "constructor like" use case and the second one "fit like" use case (as in "fit
data to a histogram").

Each of the two options imply a set of options that can go with them. The two sets of options
have some overlap, e.g. `binedges`, `overflow`; but there are also options only make sense of one of
them, e.g., `nbins` only works in the "fit like" use case.

## Make an empty histogram

The minimal information needed is `binedges`:
```@example fh
h1 = Hist1D(; binedges = 1:10)
```

You can provide more if you have them available (for example, when you're creating
histogram from another data source):

```@example fh
h2 = Hist1D(; binedges = 1:10, bincounts = collect(1:9))
```

!!! note "Default / infer strategy for empty histogram"

    To see the full list of keyword arguments as well as how they are inferred, see [`Hist1D`](@ref). To
    summarize in words:

    1. `counttype` defaults to `Float64` because `maxintfloat(Float64) ≈ 9e15`
    2. `binedges` must be provided by user
    3. `bincounts` defaults to all zero, with `eltype == counttype` and appropriate shape
    4. `sumw2` defaults to all zero with appropriate shape
    5. `nentries` defaults to 0
    6. `overflow` defaults to `false`



## Make a filled histogram given data

The minimal information needed is `data`, which is a __positional__ argument:
```@example fh
h3 = Hist1D(randn(2000))
```

You can provide more (in keyword arguments) if you want granular control over the binning behavior (e.g. number of bins,
the exact bins to use):

```@example fh
h4 = Hist1D(randn(2000); nbins=4)
```

!!! warning
    `nbins` is not strictly enforced, use `binedges` if you need exact control.

We can do non-uniform binning for example:
```@example fh
h5 = Hist1D(randn(2000); binedges = [0, 0.5, 0.8, 0.9, 1.0])
```

!!! note "Default / infer strategy for filled histogram with data"

    To see the full list of keyword arguments as well as how they are inferred, see [`Hist1D`](@ref). To
    summarize in words:

    1. data `array` must be provided by the user
    2. `counttype` defaults to `Float64` because `maxintfloat(Float64) ≈ 9e15`
    3. `nbins` defaults to sturges approximation
    4. `binedges` defaults to uniform binning given `nbins`
    5. `weights` defaults to all unity (i.e. not weighted)
    6. `overflow` defaults `false`

