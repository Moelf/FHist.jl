
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
    heights = mapreduce(h->bincounts(h), vcat, hist1ds)
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

@inline function _is_uniform_bins(A::AbstractVector{T}) where T<:Real
    diffs = diff(A)
    diff1 = first(diffs)
    all(isapprox.(diff1, diffs; atol = 1e-9)) #hack
end
function _is_uniform_bins(A::AbstractRange{T}) where T<:Real
    true
end
