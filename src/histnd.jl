auto_bins(ary, n::Val{N}; nbins = nothing) where {N} = _auto_bin_n(ary, n, nbins)

function _auto_bin_n(ary, n::Val{N}, ::Nothing) where {N}
    return ntuple(n) do i
        nbins = _sturges(ary[i])
        _auto_range(ary[i], nbins)
    end
end
function _auto_bin_n(ary, n::Val{N}, nbins::Integer) where {N}
    return ntuple(n) do i
        _auto_range(ary[i], nbins)
    end
end
function _auto_bin_n(ary, n::Val{N}, nbins::NTuple) where {N}
    length(nbins) == N || throw(ArgumentError("`nbins` must be an integer or a tuple of $N integers, got $nbins"))
    return ntuple(n) do i
        _auto_range(ary[i], nbins[i])
    end
end

"""
    sample(h::Hist; n::Int=1)

Sample a histogram's with weights equal to bin count, `n` times.
The sampled values are the bins' lower edges.
"""
function StatsBase.sample(h::HistND{T, N}; n::Int = 1) where {T, N}
    edges = binedges(h)
    counts = bincounts(h)
    cis = CartesianIndices(counts)
    sampled = StatsBase.sample(cis, Weights(vec(counts)), n)
    return ntuple(Val{N}()) do i
        [edges[i][I[i]] for I in sampled]
    end
end

"""
    nbins(h::Hist)

Get a N-tuple of the number of bins of each axes of a histogram.
"""
function nbins(h::HistND{T, N}) where {T, N}
    return size(bincounts(h))
end

function integral(h::HistND{T, N}; width = false) where {T, N}
    if width
        widths = map(diff, h.binedges)
        counts = bincounts(h)
        return sum(eachindex(IndexCartesian(), counts)) do ci
            volume = prod(zip(widths, Tuple(ci))) do (width, i)
                width[i]
            end
            counts[ci] * volume
        end
    else
        return sum(bincounts(h))
    end
end

"""
    push!(h::Hist{T, N}, vals::NTuple{N, Real}, w::Real = 1) where {T, N}

Adding one value at a time into histogram.
`sumw2` (sum of weights^2) accumulates `wgt^2` with a default weight of 1.

Note that unlike `push!` for `Hist1D`/`Hist2D`/`Hist3D`, this method requires added `vals` to be
a single `NTuple`.
"""
@inline function Base.push!(h::HistND{T, N}, vals::NTuple{N, Real}, w::Real = one(T)) where {T, N}
    Ls = nbins(h)
    is = ntuple(Val(N)) do i
        @inline
        @inbounds _binindex(h.binedges[i], Ls[i], h.overflow, vals[i])
    end
    any(==(0), is) && return nothing
    h.nentries[] += 1
    @inbounds bincounts(h)[is...] += w
    @inbounds sumw2(h)[is...] += w^2
    return nothing
end

"""
    atomic_push!(h::Hist{T, N}, vals::NTuple{N, Real}, w::Real = 1) where {T, N}

Slower but thread-safe version of [`Base.push!(::HistND, ::NTuple, ::Real)`](@ref)
"""
@inline function atomic_push!(h::HistND{T, N}, vals::NTuple{N, Real}, w::Real = one(T)) where {T, N}
    lock(h)
    push!(h, vals, w)
    unlock(h)
    return nothing
end

function Base.append!(h::HistND{T, N}, vals::NTuple{N, <:AbstractVector}, wgts::AbstractVector) where {T, N}
    allequal(length, vals) && length(wgts) == length(first(vals)) || throw(DimensionMismatch("append! to histogram expect same length values and weights"))
    lock(h)
    try
        for i in eachindex(first(vals))
            @inbounds push!(
                h, ntuple(Val(N)) do j
                    @inline
                    vals[j][i]
                end, wgts[i]
            )
        end
    finally
        unlock(h)
    end
    return h
end
function Base.append!(h::HistND{T, N}, vals::NTuple{N, <:AbstractVector}) where {T, N}
    allequal(length, vals) || throw(DimensionMismatch("append! to histogram expect same length values and weights"))
    lock(h)
    try
        for val in zip(vals...)
            @inbounds push!(h, val)
        end
    finally
        unlock(h)
    end
    return h
end

function _project1d(h::HistND{T, N}, dim::Int) where {T, N}
    dims_remove = ntuple(Val(N - 1)) do i
        i < dim ? i : i + 1
    end
    return project(h, dims_remove)
end

for op in (:mean, :std, :median)
    @eval function Statistics.$op(h::HistND{T, N}) where {T, N}
        return ntuple(Val(N)) do i
            $op(convert(Hist1D, _project1d(h, i)))
        end
    end
end

"""
    lookup(h::Hist{T, N}, vals::NTuple{N, Real}) where {T, N}

For given values `vals`, find the corresponding bin and return the bin content.
If `vals` is out of the histogram range, return `missing`.
"""
function lookup(h::HistND{T, N}, vals::NTuple{N, Any}) where {T, N}
    edges = binedges(h)
    any_out_of_range = any(Base.OneTo(N)) do i
        !(first(edges[i]) <= vals[i] < last(edges[i]))
    end
    any_out_of_range && return missing
    ids = ntuple(Val(N)) do i
        searchsortedlast(edges[i], vals[i])
    end
    return bincounts(h)[ids...]
end

"""
    normalize(h::Hist3D; width=false)

Create a normalized histogram via division by `integral(h)`. When `width==true`, each bin is
additionally divided by its volume such that `integral(normalize(h; width=true); width=true) == 1`.

!!! note
    Unlike for `Hist1D`, `width` defaults to `false` for backward compatibility.
"""
function normalize(h::HistND{T, N}; width = false) where {T, N}
    hn = h * (1 / integral(h; width = false))
    if width
        widths = map(diff, h.binedges)
        for ci in eachindex(IndexCartesian(), hn.bincounts)
            volume = prod(zip(widths, Tuple(ci))) do (width, i)
                width[i]
            end
            hn.bincounts[ci] /= volume
            hn.sumw2[ci] /= volume^2
        end
    end
    return hn
end

"""
    rebin(h::Hist{T, N}, ns::NTuple{N, Int}) where {T, N}
    rebin(h::Hist{T, N}, n::Int) where {T, N}
    rebin(ns::NTuple{N, Int}) where {T, N} = h::Hist{T, N} -> rebin(h, ns)

Merges `ns` consecutive bins into one along the axis by summing.
Alternatively, provide the new bin edges along each axis; they must be a subset of the existing
edges (see the `Hist1D` method of [`rebin`](@ref)).
"""
function rebin(h::HistND{T, N}, ns::NTuple{N, Int}) where {T, N}
    ss = nbins(h)
    bin_num_dividable = all(Base.OneTo(N)) do i
        ss[i] % ns[i] == 0
    end
    !bin_num_dividable && _rebin_error(h, ns)
    bs = h.binedges
    blocks = ntuple(Val(N)) do i
        _rebin_blocks(ss[i], ns[i])
    end
    counts = _block_sum(bincounts(h), blocks)
    s2 = _block_sum(sumw2(h), blocks)
    es = ntuple(Val(N)) do i
        _subedges(bs[i], 1:ns[i]:length(bs[i]))
    end
    return HistND{T, N}(; binedges = es, bincounts = counts, sumw2 = s2, nentries = nentries(h), overflow = h.overflow)
end
rebin(h::HistND{T, N}, n::Int) where {T, N} = rebin(h, ntuple(_ -> n, Val(N)))
rebin(ns::NTuple) = Base.Fix{2}(rebin, ns)
# Fix method ambiguity: this should not exist though
rebin(h::HistND{T, 0}, ::Tuple{}) where {T} = h

function rebin(h::HistND{T, N}, edges::NTuple{N, AbstractVector{<:Real}}) where {T, N}
    b_and_spans = ntuple(Val(N)) do i
        _edge_blocks(h.binedges[i], edges[i])
    end
    bs = map(first, b_and_spans)
    spans = map(last, b_and_spans)
    counts = _block_sum(bincounts(h), bs)
    s2 = _block_sum(sumw2(h), bs)
    return HistND{T, N}(;
        binedges = edges, bincounts = counts, sumw2 = s2,
        nentries = nentries(h), overflow = h.overflow && all(spans)
    )
end

"""
    project(h::Hist{T, N}, axis::NTuple{M, Integer}) where {T, N, M}
    project(axis::NTuple{M, Integer}) = h::Hist{T, N} -> project(h, axis)

Compute projetion of a `N`-D histogram by summing over `M`-tuple `axis`, returning `N - M` D histogram.

!!! note
    `axis` convention follows `reduce` for array, which means they are numbers of the dimension/axis
    to be reduced over. And they must be all unique unlike `reduce`.
"""
function project(h::HistND{T, N}, axis::NTuple{M, Integer}) where {T, N, M}
    allunique(axis) || throw(ArgumentError("axis must be unique, got $axis"))
    all(≤(N), axis) || throw(ArgumentError("axis must be smaller than $N, got $axis"))
    dims_all = Base.OneTo(N)
    dims_keep = filter(!∈(axis), dims_all)
    counts = dropdims(sum(bincounts(h); dims = axis); dims = axis)
    s2 = dropdims(sum(sumw2(h); dims = axis); dims = axis)
    edges = ntuple(Val(N - M)) do i
        h.binedges[dims_keep[i]]
    end
    return HistND{T, N - M}(; binedges = edges, bincounts = counts, sumw2 = s2, nentries = nentries(h), overflow = h.overflow)
end
project(axis::NTuple{M, Integer}) where {M} = Base.Fix{2}(project, axis)

"""
    restrict(h::Hist{T, N}, lows = (-Inf, -Inf, ...), highs = (Inf, Inf, ...))

Returns a new histogram with restricted axes: the slice of `h` where the bin centers are within
the given (inclusive) intervals along each axis.
"""
function restrict(h::HistND{T, N}, lows = ntuple(_ -> -Inf, Val(N)), highs = ntuple(_ -> Inf, Val(N))) where {T, N}
    lows::NTuple{N, T} = lows
    highs::NTuple{N, T} = highs
    bs = h.binedges
    sels = ntuple(Val(N)) do i
        _restrict_bins(bs[i], lows[i], highs[i])
    end
    edges = ntuple(Val(N)) do i
        _subedges(bs[i], first(sels[i]):(last(sels[i]) + 1))
    end
    c = bincounts(h)[sels...]
    s2 = sumw2(h)[sels...]
    return HistND{T, N}(; binedges = edges, bincounts = c, sumw2 = s2, nentries = nentries(h), overflow = h.overflow)
end
# Duplicated methods with `Hist1D`
#restrict(lows, highs) = Base.Fix{2}(Base.Fix{3}(restrict, highs), lows)

function Base.convert(::Type{Hist1D}, h::HistND{T, 1}) where {T}
    return Hist1D(; counttype = T, binedges = only(binedges(h)), bincounts = bincounts(h), sumw2 = sumw2(h), nentries = nentries(h), overflow = h.overflow)
end
function Base.convert(::Type{Hist2D}, h::HistND{T, 2}) where {T}
    return Hist2D(; counttype = T, binedges = binedges(h), bincounts = bincounts(h), sumw2 = sumw2(h), nentries = nentries(h), overflow = h.overflow)
end
function Base.convert(::Type{Hist3D}, h::HistND{T, 3}) where {T}
    return Hist3D(; counttype = T, binedges = binedges(h), bincounts = bincounts(h), sumw2 = sumw2(h), nentries = nentries(h), overflow = h.overflow)
end

"""
    slices(h::HistND{T, N}, axis::Integer)

Get N-1 dim sliced histogram along `axis`. Returns an iterator of tuples of sliced axis edges (`(low, high)`) and sliced histogram.

```jldoctest; setup = :(using FHist)
julia> h = HistND{Float64, 2}(([0.25, 0.75], [0.25, 0.75]); binedges=(0:0.5:1, 0:0.5:1))
HistND{Float64, 2}, edges=(0.0:0.5:1.0, 0.0:0.5:1.0), integral=2.0

julia> collect(slices(h, 2))
2-element Vector{Tuple{Tuple{Float64, Float64}, HistND{Float64, 1}}}:
 ((0.0, 0.5), HistND{Float64, 1}, edges=(0.0:0.5:1.0,), integral=1.0)
 ((0.5, 1.0), HistND{Float64, 1}, edges=(0.0:0.5:1.0,), integral=1.0)
```
"""
function slices(h::HistND{T, N}, axis::Integer) where {T, N}
    binedges_sliced = h.binedges[axis]
    return Iterators.map(1:(length(binedges_sliced) - 1)) do i
        lows = ntuple(Val(N)) do j
            @inline
            j == axis ? binedges_sliced[i] : -Inf
        end
        highs = ntuple(Val(N)) do j
            @inline
            j == axis ? binedges_sliced[i + 1] : Inf
        end
        h_restricted = restrict(h, lows, highs)
        (binedges_sliced[i], binedges_sliced[i + 1]), project(h_restricted, (axis,))
    end
end
