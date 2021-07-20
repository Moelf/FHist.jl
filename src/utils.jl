
"""
    hists_to_bars(hist1ds)

Given a vector of Hist1D, return `edges` (xs), `heights` (ys),
and `grps` (for grouping) that is useful for plotting stacked
histogram.
"""
function hists_to_bars(hist1ds)
    L = length(hist1ds)
    oneedge = hist1ds[1].hist.edges[1][1:end-1]
    edges = repeat(oneedge, L)
    heights = mapreduce(h->h.hist.weights, vcat, hist1ds)
    grps = repeat(1:L; inner=length(oneedge))
    
    edges, heights, grps
end
