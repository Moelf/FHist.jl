"""
    hists_to_bars(hist1ds)

Given a vector of Hist1D, return `edges` (xs), `heights` (ys),
and `grps` (for grouping) that is useful for plotting stacked
histogram.
"""
function hists_to_bars(hist1ds)
    L = length(hist1ds)
    oneedge = binedges(hist1ds[1])[1:end-1]
    edges = repeat(oneedge, L)
    heights = mapreduce(bincounts, vcat, hist1ds)
    grps = repeat(1:L; inner=length(oneedge))

    edges, heights, grps
end

@inline function pearson_err(n::Real)
    s = sqrt(n+0.25)
    s+0.5, s-0.5
end
@inline function sqrt_err(n::Real)
    s = sqrt(n)
    s,s
end

_sturges(x) = StatsBase.sturges(length(x))

# `nbins` keyword: a single integer applies to every axis, otherwise one entry per axis
_nbins_per_axis(nbins::Integer, ::Val{N}) where {N} = ntuple(_ -> Int(nbins), Val(N))
function _nbins_per_axis(nbins, ::Val{N}) where {N}
    length(nbins) == N || throw(ArgumentError("`nbins` must be an integer or a tuple of $N integers, got $nbins"))
    return ntuple(i -> Int(nbins[i]), Val(N))
end

_float_eltype(xs) = (E = eltype(xs); E <: Number ? float(E) : Float64)

function _auto_range(xs, nbins)
    F = _float_eltype(xs)
    isempty(xs) && throw(ArgumentError("cannot infer bin edges from empty data, provide `binedges`"))
    lo, hi = extrema(xs)
    if !(isfinite(lo) && isfinite(hi))
        # ignore NaN/Inf for the purpose of choosing the edges (they are discarded / clamped when filling)
        finite = filter(isfinite, xs)
        isempty(finite) && throw(ArgumentError("cannot infer bin edges from data without finite values, provide `binedges`"))
        lo, hi = extrema(finite)
    end
    return StatsBase.histrange(F(lo), F(hi), nbins)
end

@inline function _is_uniform_bins(A::AbstractVector{T}) where T<:Real
    diffs = diff(A)
    diff1 = first(diffs)
    all(isapprox.(diff1, diffs; atol = 1e-9)) #hack
end
function _is_uniform_bins(A::AbstractRange{T}) where T<:Real
    true
end

# A sub-selection (by index) of bin edges that keeps ranges as ranges (for display)
_subedges(b::BinEdges, idx::AbstractRange{<:Integer}) = b.isrange ? b.range[idx] : b.edges[idx]

# Index ranges of the old bins that make up each new bin when merging `n` consecutive bins
_rebin_blocks(L::Int, n::Int) = [((k-1)*n+1):(k*n) for k in 1:(L ÷ n)]

# Index ranges of the old bins that make up each new bin, for new edges that are a subset of
# the old ones. Also returns whether the new edges span the full old range.
function _edge_blocks(old::BinEdges, new_edges::AbstractVector{<:Real})
    length(new_edges) >= 2 || throw(ArgumentError("`edges` must contain at least two elements"))
    allunique(new_edges) || throw(ArgumentError("`edges` must be all unique"))
    issorted(new_edges) || throw(ArgumentError("`edges` must be sorted"))
    idx = map(new_edges) do ne
        i = findfirst(==(ne), old.edges)
        isnothing(i) && throw(ArgumentError("`edges` must be composed of existing histogram bin edges, $ne not found"))
        i
    end
    blocks = [idx[i]:(idx[i+1]-1) for i in 1:length(idx)-1]
    spans_all = first(idx) == 1 && last(idx) == length(old)
    return blocks, spans_all
end

# Sum `A` over the blocks (index ranges) along each axis
function _block_sum(A::AbstractArray{T,N}, blocks::NTuple{N,Vector{UnitRange{Int}}}) where {T,N}
    out = similar(A, length.(blocks))
    for I in CartesianIndices(out)
        @inbounds out[I] = sum(view(A, ntuple(d -> blocks[d][I[d]], Val(N))...))
    end
    return out
end

# Index range of the bins whose centers are within [low, high]
function _restrict_bins(b::BinEdges, low, high)
    sel = low .<= StatsBase.midpoints(b.edges) .<= high
    @assert count(sel) > 0 "No bin centers contained in [$(low), $(high)]"
    return findfirst(sel):findlast(sel)
end

_range_or_vector(b::BinEdges) = b.isrange ? b.range : b.edges
Base.convert(::Type{StatsBase.Histogram}, h::Union{Hist1D, Hist2D, Hist3D}) = StatsBase.Histogram(map(_range_or_vector, h.binedges), bincounts(h))

"""
    valid_rebin_values(h::Union{Hist1D, Hist2D, Hist3D})

Calculates the legal values for rebinning, essentially the prime factors of
the number of bins. For a 1D histogram, a `Set` of numbers is return, for higher
dimensional histograms a `Vector{Set}` for each dimension.
"""
valid_rebin_values(h::Hist1D) = _factor(nbins(h))
valid_rebin_values(h::Union{Hist2D, Hist3D}) = [_factor(x) for x in nbins(h)]

"""
    _factor(n::Integer)

Helper function to calculate the prime factors of a given integer.
"""
function _factor(n::Integer)
    factors = Set{Int}()
    n = Int(n)
    factor = 2
    while n > 1 && factor * factor <= n
        while n % factor == 0
            n ÷= factor
            push!(factors, factor)
        end
        factor += 1
    end
    n > 1 && push!(factors, n)
    return factors
end

function _rebin_error(h, ns)
    vals = valid_rebin_values(h)
    valid = map(x -> join(sort(collect(x)), ", ", " or "), vals isa Set ? (vals,) : Tuple(vals))
    error("Invalid rebin value(s) $(ns) for a histogram with $(nbins(h)) bins. Each has to be a divisor of the number of bins along that axis, i.e. a product of the prime factors: $(valid)")
end
