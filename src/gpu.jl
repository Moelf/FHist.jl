# GPU support lives in the `FHistGPUExt` package extension, which is loaded automatically when
# both `KernelAbstractions` and `GPUArraysCore` are loaded, i.e. whenever a GPU array package
# such as CUDA.jl, AMDGPU.jl, Metal.jl or oneAPI.jl is loaded. Only the docstrings and the
# generic function stubs are defined here.

export gpu_bincounts, gpu_bincounts!

"""
    gpu_bincounts(data::AbstractArray; binedges, weights=nothing, overflow=false,
                  counttype=Float32, backend=get_backend(data), blocksize=256, sync=true)

Compute the 1D bin counts of `data` on the device it lives on (GPU or the KernelAbstractions
`CPU()` backend) and return them as a new device array of element type `counttype` and length
`length(binedges) - 1`.

Requires a GPU array package (CUDA.jl, AMDGPU.jl, Metal.jl, oneAPI.jl, ...) to be loaded.

- `binedges` can be an `AbstractRange` or a vector, uniform binnings use an O(1) lookup.
- `weights` must be `nothing` or a device array of the same length as `data`.
- `overflow` follows the semantics of `Hist1D`: values outside of the edges are discarded, or
  clamped into the first/last bin when `overflow=true`.
- `counttype` is the element type of the counts (`Float32` by default, as not every device
  supports `Float64` atomics; `Int32` also works for unweighted data).
- `blocksize` is the work-group size, `sync=false` skips the final `synchronize(backend)`.

To get a full [`Hist1D`](@ref) (with `sumw2` and `nentries`) from GPU data, simply call the
regular constructor with a GPU array: `Hist1D(cu_data; binedges = 0:0.1:1)`, the result lives on
the CPU.

See also [`gpu_bincounts!`](@ref).

# Example
```julia
using FHist, CUDA
xs = CUDA.randn(10^7)
counts = gpu_bincounts(xs; binedges = -3:0.1:3)            # CuVector{Float32}
h = Hist1D(xs; binedges = -3:0.1:3, counttype = Float32)   # Hist1D on the CPU
```
"""
function gpu_bincounts end

"""
    gpu_bincounts!(counts::AbstractArray, data::AbstractArray; binedges, weights=nothing,
                   overflow=false, backend=get_backend(data), blocksize=256, sync=true)

In-place version of [`gpu_bincounts`](@ref): accumulate the bin counts of `data` into the device
array `counts` (which must have length `length(binedges) - 1` and is *not* zeroed first).
Returns `counts`.
"""
function gpu_bincounts! end
