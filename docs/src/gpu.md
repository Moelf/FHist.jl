# GPU histogramming

FHist.jl can fill 1D histograms directly from data that lives on a GPU, using
[KernelAbstractions.jl](https://github.com/JuliaGPU/KernelAbstractions.jl) so that the same
kernels run on every supported vendor. The functionality is implemented as a package extension
that is loaded automatically whenever a GPU array package is loaded (CUDA.jl, AMDGPU.jl,
Metal.jl, oneAPI.jl, ...); no extra `using` is needed beyond the GPU package itself.

## Filling a `Hist1D` from GPU data

The regular "fit like" constructor accepts GPU arrays. The binning happens on the device and the
resulting [`Hist1D`](@ref) lives on the CPU like any other FHist histogram:

```julia
using FHist, CUDA   # or Metal, AMDGPU, oneAPI

xs = CUDA.randn(Float32, 10^8)
ws = CUDA.rand(Float32, 10^8)

h  = Hist1D(xs; binedges = -3:0.1:3, counttype = Float32)
hw = Hist1D(xs; binedges = -3:0.1:3, weights = ws, counttype = Float32, overflow = true)
```

Everything a CPU-filled histogram has is filled in: `bincounts`, `sumw2` and `nentries`, and the
result is identical to what `Hist1D(Array(xs); ...)` gives (up to floating point summation
order for weighted data). The semantics for values outside of the edges (`overflow`), `NaN` and
`Inf` are the same as on the CPU, see [`atomic_push!`](@ref).

!!! note "Count type"
    The counts are accumulated on the device with atomics of the histogram's `counttype`.
    Not every device supports `Float64` atomics (Metal does not support `Float64` at all), so
    pass `counttype = Float32` (or `Int32` for unweighted data) on such devices. For weighted
    data `sumw2` is accumulated in `float(counttype)`.

!!! note "Bin edges precision"
    On the device, the bin edges are stored in the floating point type of the data
    (`Float32` for `Float32` data). They are rounded *up* when converting so that a value
    `x` compares against the edge exactly like the `Float64` comparison on the CPU does; the
    bin assignment is therefore identical to the CPU one. Bin edges that are so close that they
    become equal in `Float32` are rejected.

## Lower level: `gpu_bincounts`

If you only need the counts, and want to keep them on the device (e.g. to keep processing them
there, or to accumulate over many chunks of data), use [`gpu_bincounts`](@ref) and
[`gpu_bincounts!`](@ref):

```julia
counts = gpu_bincounts(xs; binedges = -3:0.1:3)               # CuVector{Float32} of length 60
gpu_bincounts!(counts, more_xs; binedges = -3:0.1:3)           # accumulate more data into it
Array(counts)                                                  # bring it back to the CPU
```

Both functions also accept plain `Array`s, in which case the kernels run on the
KernelAbstractions `CPU()` backend (this is how the extension is tested on CI):

```julia
using KernelAbstractions, GPUArraysCore   # what a GPU package would load for you
gpu_bincounts(randn(10^6); binedges = -3:0.1:3) == bincounts(Hist1D(randn(10^6); binedges = -3:0.1:3))
```

### Options

- `binedges`: an `AbstractRange` or a vector, uniform binnings (also when given as a vector)
  use an O(1) lookup, non-uniform ones a binary search.
- `weights`: `nothing` or a device array of the same length as the data.
- `overflow`: as for `Hist1D`, values outside of the edges are discarded, or clamped into the
  first/last bin when `true`.
- `counttype` (only `gpu_bincounts`): element type of the returned counts, `Float32` by default.
- `blocksize`: the work-group size, `256` by default.
- `naive`: by default a work-group first accumulates into a shared-memory histogram and
  merges it into the global one afterwards, which is much faster than one global atomic per
  entry. When the shared-memory buffers would exceed 16 KiB (e.g. more than 4096 `Float32`
  bins), global atomics (`naive = true`) are used automatically. Pass `naive = true/false` to
  force either algorithm.
- `edgetype`: floating point type of the edges on the device, `float(eltype(data))` by default.
- `sync`: whether to `synchronize(backend)` before returning (default `true`).

## Performance

On an Apple M4 (10 GPU cores), histogramming `10^8` `Float32` values into 600 bins with
`gpu_bincounts` takes about 15 ms (the same with `naive = true`, global atomics are fast on
Apple silicon; on discrete GPUs the shared-memory version is typically much faster), and
`Hist1D(xs; ...)` about 30 ms including the transfers back to the host, compared to about 90 ms
for the single-threaded CPU `Hist1D` on the same data. The numbers depend strongly on the
device and on the number of bins.
