module FHistGPUExt

# GPU (and KernelAbstractions `CPU()`) histogramming for `Hist1D`.
#
# This extension is loaded automatically when `KernelAbstractions` and `GPUArraysCore` are both
# loaded, which is the case whenever a GPU array package (CUDA.jl, AMDGPU.jl, Metal.jl,
# oneAPI.jl, ...) is loaded.

using FHist
using FHist: BinEdges, nbins, bincounts, sumw2, _uniform_lookup_ok
using GPUArraysCore: AbstractGPUArray
using KernelAbstractions
using KernelAbstractions: KernelAbstractions as KA
using KernelAbstractions: @atomic

# Shared-memory histograms are used automatically as long as the per-work-group buffers fit in
# `DEFAULT_SHMEM_BYTES`; requesting them explicitly (`naive=false`) is refused above
# `MAX_SHMEM_BYTES` (the limit of e.g. Metal devices) since exceeding it crashes the compiler.
const DEFAULT_SHMEM_BYTES = 16 * 1024
const MAX_SHMEM_BYTES = 32 * 1024

# ---------------------------------------------------------------------------------------------
# Bin lookup, identical to `FHist._binindex` on the CPU (see `polybinedges.jl`)

# `searchsortedlast(edges, x)` with `edges[1] = rfirst`, `edges[end] = rlast`, `L = length - 1`
@inline function _gpu_searchsortedlast(x::F, edges, L::Int32, isuniform::Bool, rfirst::F, rlast::F, inv_step::F) where {F}
    x < rfirst && return Int32(0)
    x < rlast || return L + Int32(1)  # x >= last, or NaN
    if isuniform
        g = unsafe_trunc(Int32, (x - rfirst) * inv_step) + Int32(1)
        g = clamp(g, Int32(1), L)
        @inbounds if x < edges[g]
            g -= Int32(1)
        elseif x >= edges[g+Int32(1)]
            g += Int32(1)
        end
        return g
    else
        # invariant: edges[lo] <= x < edges[hi]
        lo = Int32(1)
        hi = L + Int32(1)
        while hi - lo > Int32(1)
            mid = (lo + hi) >>> Int32(1)
            @inbounds if edges[mid] <= x
                lo = mid
            else
                hi = mid
            end
        end
        return lo
    end
end

@inline function _gpu_binindex(x, edges, L::Int32, isuniform::Bool, rfirst::F, rlast::F, inv_step::F, overflow::Bool) where {F}
    i = _gpu_searchsortedlast(F(x), edges, L, isuniform, rfirst, rlast, inv_step)
    if overflow
        return clamp(i, Int32(1), L)
    else
        return (Int32(1) <= i <= L) ? i : Int32(0)
    end
end

@inline _weight(::Type{T}, ::Nothing, tid) where {T} = one(T)
@inline _weight(::Type{T}, weights, tid) where {T} = @inbounds T(weights[tid])

# ---------------------------------------------------------------------------------------------
# Kernels. `S` (a `Val{Bool}`) says whether `sumw2` is accumulated; `nentries` has one slot per
# work-group and is *set* (not accumulated) by the kernel.

@kernel unsafe_indices = true function _hist1d_naive_kernel!(counts, s2, nentries, ::Val{S},
    @Const(data), @Const(weights), @Const(edges), L::Int32, isuniform::Bool, rfirst, rlast, inv_step, overflow::Bool) where {S}
    gid = @index(Group, Linear)
    lid = @index(Local, Linear)
    gs = Int32(prod(@groupsize()))
    tid = (gid - Int32(1)) * gs + lid
    T = eltype(counts)

    if tid <= length(data)
        bin = _gpu_binindex(@inbounds(data[tid]), edges, L, isuniform, rfirst, rlast, inv_step, overflow)
        if bin != Int32(0)
            w = _weight(T, weights, tid)
            @atomic counts[bin] += w
            if S
                @atomic s2[bin] += eltype(s2)(w) * eltype(s2)(w)
            end
            @atomic nentries[gid] += Int32(1)
        end
    end
end

@kernel unsafe_indices = true function _hist1d_shmem_kernel!(counts, s2, nentries, ::Val{S}, ::Val{Nb}, ::Val{Nb2},
    @Const(data), @Const(weights), @Const(edges), isuniform::Bool, rfirst, rlast, inv_step, overflow::Bool) where {S,Nb,Nb2}
    @uniform L = Int32(Nb)
    @uniform gs = Int32(prod(@groupsize()))

    # (`@localmem` declarations are hoisted to the top of the kernel, so no local aliases here)
    sh_counts = @localmem eltype(counts) (Nb)
    sh_s2 = @localmem eltype(s2) (Nb2)  # Nb2 == Nb when S, 1 otherwise
    sh_n = @localmem Int32 (1)

    # zero the shared buffers
    lid = @index(Local, Linear)
    i = lid
    while i <= L
        @inbounds sh_counts[i] = zero(eltype(counts))
        if S
            @inbounds sh_s2[i] = zero(eltype(s2))
        end
        i += gs
    end
    if lid == Int32(1)
        @inbounds sh_n[1] = Int32(0)
    end
    @synchronize()

    # bin this work-item's value into shared memory
    gid = @index(Group, Linear)
    lid = @index(Local, Linear)
    tid = (gid - Int32(1)) * gs + lid
    if tid <= length(data)
        bin = _gpu_binindex(@inbounds(data[tid]), edges, L, isuniform, rfirst, rlast, inv_step, overflow)
        if bin != Int32(0)
            w = _weight(eltype(counts), weights, tid)
            @atomic sh_counts[bin] += w
            if S
                @atomic sh_s2[bin] += eltype(s2)(w) * eltype(s2)(w)
            end
            @atomic sh_n[1] += Int32(1)
        end
    end
    @synchronize()

    # flush shared memory into the global histogram (skipping empty bins)
    gid = @index(Group, Linear)
    lid = @index(Local, Linear)
    i = lid
    while i <= L
        c = @inbounds sh_counts[i]
        if c != zero(eltype(counts))
            @atomic counts[i] += c
        end
        if S
            v = @inbounds sh_s2[i]
            if v != zero(eltype(s2))
                @atomic s2[i] += v
            end
        end
        i += gs
    end
    if lid == Int32(1)
        @inbounds nentries[gid] = sh_n[1]
    end
end

# ---------------------------------------------------------------------------------------------
# Host side

# Device-friendly copy of the bin edges in precision `F`. Edges are rounded *up* when converting to
# a narrower float type so that `x >= edge` gives the same answer for every `x::F` as the exact
# comparison against the Float64 edge does.
function _device_edges(b::BinEdges, ::Type{F}, backend) where {F}
    v = F === Float64 ? copy(b.edges) : F[F(e, RoundUp) for e in b.edges]
    for i in 1:length(v)-1
        v[i] < v[i+1] || throw(ArgumentError("bin edges are not distinct anymore after conversion to $F: $(v[i]) and $(v[i+1])"))
    end
    inv_step = F(length(v) - 1) / (last(v) - first(v))
    isuniform = b.isuniform && _uniform_lookup_ok(v, inv_step)
    d_edges = KA.allocate(backend, F, length(v))
    copyto!(d_edges, v)
    return d_edges, isuniform, first(v), last(v), inv_step
end

# Fill `counts` (and `s2` unless `nothing`) with the histogram of `data`; returns the number of
# accepted entries. `counts`/`s2` are accumulated into (not zeroed).
function _gpu_fill!(counts, s2, data, weights, b::BinEdges, overflow::Bool;
        backend=get_backend(data), blocksize::Int=256, sync::Bool=true, naive::Union{Nothing,Bool}=nothing,
        edgetype::Type=float(eltype(data)))
    L = length(counts)
    L == length(b) - 1 || throw(DimensionMismatch("`counts` must have `length(binedges) - 1 == $(length(b) - 1)` elements, got $L"))
    L <= typemax(Int32) || throw(ArgumentError("too many bins"))
    if !isnothing(weights)
        length(weights) == length(data) || throw(DimensionMismatch("data and weights must have the same length"))
        get_backend(weights) == backend || throw(ArgumentError("weights must live on the same backend as the data"))
    end
    get_backend(counts) == backend || throw(ArgumentError("`counts` must live on the same backend as the data"))
    S = !isnothing(s2)
    S && (length(s2) == L || throw(DimensionMismatch("`sumw2` must have the same length as `counts`")))
    s2_arg = S ? s2 : counts

    d_edges, isuniform, rfirst, rlast, inv_step = _device_edges(b, edgetype, backend)
    N = length(data)
    ngroups = max(1, cld(N, blocksize))
    nentries = KA.zeros(backend, Int32, ngroups)

    shmem_bytes = L * (sizeof(eltype(counts)) + (S ? sizeof(eltype(s2_arg)) : 0))
    use_naive = isnothing(naive) ? shmem_bytes > DEFAULT_SHMEM_BYTES : naive
    if !use_naive && shmem_bytes > MAX_SHMEM_BYTES
        throw(ArgumentError("the shared-memory histogram needs $shmem_bytes bytes per work-group which exceeds the $MAX_SHMEM_BYTES bytes limit, use `naive=true` (global atomics) instead"))
    end
    if N > 0
        if use_naive
            kernel! = _hist1d_naive_kernel!(backend, blocksize)
            kernel!(counts, s2_arg, nentries, Val(S), data, weights, d_edges, Int32(L), isuniform, rfirst, rlast, inv_step, overflow;
                ndrange=ngroups * blocksize)
        else
            kernel! = _hist1d_shmem_kernel!(backend, blocksize)
            kernel!(counts, s2_arg, nentries, Val(S), Val(L), Val(S ? L : 1), data, weights, d_edges, isuniform, rfirst, rlast, inv_step, overflow;
                ndrange=ngroups * blocksize)
        end
    end
    n = Int(sum(Int64, Array(nentries)))  # also acts as a synchronization point for the counts
    sync && synchronize(backend)
    return n
end

_binedges(binedges) = binedges isa BinEdges ? binedges : BinEdges(binedges)

function FHist.gpu_bincounts!(counts::AbstractArray, data::AbstractArray; binedges, weights=nothing, overflow::Bool=false,
        backend=get_backend(data), blocksize::Int=256, sync::Bool=true, naive::Union{Nothing,Bool}=nothing,
        edgetype::Type=float(eltype(data)))
    _gpu_fill!(counts, nothing, data, weights, _binedges(binedges), overflow; backend, blocksize, sync, naive, edgetype)
    return counts
end

function FHist.gpu_bincounts(data::AbstractArray; binedges, weights=nothing, overflow::Bool=false, counttype::Type=Float32,
        backend=get_backend(data), blocksize::Int=256, sync::Bool=true, naive::Union{Nothing,Bool}=nothing,
        edgetype::Type=float(eltype(data)))
    b = _binedges(binedges)
    counts = KA.zeros(backend, counttype, length(b) - 1)
    return FHist.gpu_bincounts!(counts, data; binedges=b, weights, overflow, backend, blocksize, sync, naive, edgetype)
end

# Fill an (empty) `Hist1D` living on the host from device `data`/`weights`.
function _gpu_hist1d!(h::Hist1D, data::AbstractArray, weights; backend=get_backend(data), kws...)
    T = eltype(bincounts(h))
    b = h.binedges[1]
    d_counts = KA.zeros(backend, T, nbins(h))
    if isnothing(weights)
        n = _gpu_fill!(d_counts, nothing, data, nothing, b, h.overflow; backend, kws...)
        bincounts(h) .+= Array(d_counts)
        sumw2(h) .= bincounts(h)  # valid because `h` starts empty
    else
        T2 = float(T)
        d_s2 = KA.zeros(backend, T2, nbins(h))
        n = _gpu_fill!(d_counts, d_s2, data, weights, b, h.overflow; backend, kws...)
        bincounts(h) .+= Array(d_counts)
        sumw2(h) .+= Array(d_s2)
    end
    h.nentries[] += n
    return h
end

# `Hist1D(gpu_array; binedges = ...)` goes through here
function FHist._fast_bincounts!(h::Hist1D, A::Tuple{<:AbstractGPUArray}, weights)
    isnothing(weights) || weights isa AbstractGPUArray || throw(ArgumentError("`weights` must be a GPU array when the data is a GPU array"))
    return _gpu_hist1d!(h, A[1], weights)
end

end
