function auto_bins(ary, ::Val{1}; nbins=nothing)
    xs = only(ary)
    nb = isnothing(nbins) ? _sturges(xs) : only(_nbins_per_axis(nbins, Val(1)))
    return (_auto_range(xs, nb),)
end

"""
    sample(h::Hist1D; n::Int=1)

Sample a histogram's with weights equal to bin count, `n` times.
The sampled values are the bins' lower edges.
"""
StatsBase.sample(h::Hist1D; n::Int=1) = StatsBase.sample(binedges(h)[1:end-1], Weights(bincounts(h)), n)

"""
    nbins(h::Hist1D)

Get the number of bins of a histogram.
"""
function nbins(h::Hist1D)
    length(bincounts(h))
end

"""
    integral(h; width=false)

Get the integral a histogram; `width` means multiply each bincount
by their bin width (bin area for `Hist2D`, bin volume for `Hist3D`) when calculating the integral.

!!! warning
    Be aware of the approximation you make
    when using `width=true` with histogram with overflow bins, the overflow
    bins (i.e. the left/right most bins) width will be taken "as is".
"""
function integral(h::Hist1D; width=false)
    if width
        mapreduce(*, +, bincounts(h), diff(binedges(h)))
    else
        sum(bincounts(h))
    end
end

"""
    push!(h::Hist1D, val::Real, wgt::Real=1)
    atomic_push!(h::Hist1D, val::Real, wgt::Real=1)

Adding one value at a time into histogram.
`sumw2` (sum of weights^2) accumulates `wgt^2` with a default weight of 1.
`atomic_push!` is a slower version of `push!` that is thread-safe.

Values outside of the bin edges are discarded (and not counted in `nentries`), unless the
histogram was created with `overflow=true`, in which case they are clamped into the first/last
bin. `NaN` is treated like `+Inf`.

N.B. To append multiple values at once, use broadcasting via
`push!.(h, [-3.0, -2.9, -2.8])` or `push!.(h, [-3.0, -2.9, -2.8], 2.0)`, or [`append!`](@ref).
"""
@inline function atomic_push!(h::Hist1D, val::Real, wgt::Real=1)
    lock(h)
    push!(h, val, wgt)
    unlock(h)
    return nothing
end

@inline function Base.push!(h::Hist1D, val::Real, wgt::Real=1)
    i = _binindex(h.binedges[1], nbins(h), h.overflow, val)
    i == 0 && return nothing
    h.nentries[] += 1
    @inbounds bincounts(h)[i] += wgt
    @inbounds sumw2(h)[i] += wgt^2
    return nothing
end

"""
    append!(h::Hist1D, vals[, wgts])
    append!(h::Hist2D, xs, ys[, wgts])
    append!(h::Hist3D, xs, ys, zs[, wgts])

`push!` many values (optionally with weights) into the histogram at once. The histogram lock
is held for the duration of the call, so this is thread-safe like [`atomic_push!`](@ref).
Returns `h`.
"""
function Base.append!(h::Hist1D, val::AbstractVector, wgt::AbstractVector)
    length(val) == length(wgt) || throw(DimensionMismatch("append! to histogram expect same length values and weights"))
    lock(h)
    try
        for (v, w) in zip(val, wgt)
            push!(h, v, w)
        end
    finally
        unlock(h)
    end
    return h
end

function Base.append!(h::Hist1D, val::AbstractVector)
    lock(h)
    try
        for v in val
            push!(h, v)
        end
    finally
        unlock(h)
    end
    return h
end

"""
    Hist1D(data::AbstractVector; kws...)

Convenience method: a non-tuple `data` is wrapped into a 1-tuple, see [`Hist1D`](@ref) for the
keyword arguments.
"""
function Hist1D(ary; kws...)
    Hist1D((ary, ); kws...)
end

"""
    Statistics.mean(h)
    Statistics.std(h)
    Statistics.median(h)
    Statistics.quantile(h::Hist1D, p)

Compute statistical quantities based on the bin centers weighted
by the bin counts.

When the histogram is `Hist2D` (`Hist3D`), return a 2-tuple (3-tuple) instead, e.g
`(mean(project(h, :x)), mean(project(h, :y)))` etc.
"""
Statistics.mean(h::Hist1D) = Statistics.mean(bincenters(h), Weights(bincounts(h)))
Statistics.std(h::Hist1D) = begin
    var = Statistics.var(bincenters(h), Weights(bincounts(h)); corrected=false)
    var < 0 ? NaN : sqrt(var)
end
Statistics.median(h::Hist1D) = Statistics.median(bincenters(h), Weights(bincounts(h)))
Statistics.quantile(h::Hist1D, p) = Statistics.quantile(bincenters(h), Weights(bincounts(h)), p)

"""
    lookup(h::Hist1D, x)

For given x-axis value `x`, find the corresponding bin and return the bin content.
If a value is out of the histogram range, return `missing`.
"""
function lookup(h::Hist1D, x)
    r = binedges(h)
    !(first(r) <= x < last(r)) && return missing
    return bincounts(h)[searchsortedlast(r, x)]
end

"""
    normalize(h::Hist1D; width=true)

Create a normalized histogram via division by `integral(h)`, when `width==true`, the
resultant histogram has area under the curve equals 1.

!!! warning
    Implicit approximation is made when using `width=true` with histograms
    that have overflow bins: the overflow data lives inthe left/right most bins
    and the bin width is taken "as is".
"""
function normalize(h::Hist1D; width=true)
    h_prob_normalized = h*(1/integral(h; width=false))
    if width
        h_prob_normalized.bincounts ./= diff(binedges(h_prob_normalized))
    end
    return h_prob_normalized
end

"""
    cumulative(h::Hist1D; forward=true)

Create a cumulative histogram. If `forward`, start
summing from the left.
"""
function cumulative(h::Hist1D; forward=true)
    # https://root.cern.ch/doc/master/TH1_8cxx_source.html#l02608
    f = forward ? identity : reverse
    h = deepcopy(h)
    bc = bincounts(h)
    bc .= f(cumsum(f(bc)))

    s2 = sumw2(h)
    s2 .= f(cumsum(f(s2)))
    return h
end


"""
    rebin(h::Hist1D, n::Int=1)
    rebin(h::Hist1D, edges::AbstractVector{<:Real})
    rebin(n::Int)
    rebin(edges::AbstractVector{<:Real})

Rebin a histogram by merging existing bins. When provided an integer `n`, the
function merges `n` consecutive bins and returns `nbins(h) / n` bins. When
provided a collection of bin edges `edges`, the function returns a new
histogram whose bin edges match `edges`; every element of `edges` must align
with the original bin edges.

If the `edges` is an array and doesn't include original histogram's leftmost
and rightmost edges, those bins will be ignored (and `overflow` is set to `false`
for the result).

The curried forms `rebin(n)` / `rebin(edges)` return a function `h -> rebin(h, ...)`,
they also work for `Hist2D` and `Hist3D` (using `n` along every axis).
"""
function rebin(h::Hist1D, n::Int=1)
    nbins(h) % n == 0 || _rebin_error(h, n)
    b = h.binedges[1]
    blocks = _rebin_blocks(nbins(h), n)
    counts = _block_sum(bincounts(h), (blocks,))
    s2 = _block_sum(sumw2(h), (blocks,))
    edges = _subedges(b, 1:n:length(b))
    return Hist1D(; binedges = edges, bincounts = counts, sumw2 = s2, nentries = nentries(h), overflow = h.overflow)
end
rebin(n::Int) = h -> rebin(h, n)

function rebin(h::Hist1D, new_edges::AbstractVector{<:Real})
    blocks, spans_all = _edge_blocks(h.binedges[1], new_edges)
    counts = _block_sum(bincounts(h), (blocks,))
    s2 = _block_sum(sumw2(h), (blocks,))
    return Hist1D(; binedges = new_edges, bincounts = counts, sumw2 = s2,
        nentries = nentries(h), overflow = h.overflow && spans_all)
end
rebin(edges::AbstractVector{<:Real}) = h -> rebin(h, edges)

"""
    bayes_rebin_edges(h::Hist1D; prior=BayesHistogram.Geometric(0.995))

Find optimal bin edges for a histogram using Bayesian rebinning algorithm.
This function only find edges, it doesn't return a new histogram.

For possible priors, see [`BayesHistogram.jl`](https://github.com/francescoalemanno/BayesHistogram.jl/blob/main/src/BayesHistogram.jl).

"""
function bayes_rebin_edges(h::Hist1D; prior=BayesHistogram.Geometric(0.995))
    old_edges = binedges(h)
    length(old_edges) < 4 && error("too little bins to rebin")
    fake_xs = [first(old_edges); bincenters(h); last(old_edges)]
    weights = [0; bincounts(h); 0]
    res = BayesHistogram.bayesian_blocks(fake_xs; weights=weights, prior=prior)
    return res.edges
end

"""
    restrict(h::Hist1D, low=-Inf, high=Inf)
    restrict(low=-Inf, high=Inf) = h::Hist1D -> restrict(h, low, high)

Returns a new histogram with a restricted x-axis.
`restrict(h, 0, 3)` (or `h |> restrict(0, 3)`)
will return a slice of `h` where the bin centers are in `[0, 3]` (inclusive).
"""
function restrict(h::Hist1D, low=-Inf, high=Inf)
    sel = _restrict_bins(h.binedges[1], low, high)
    edges = _subedges(h.binedges[1], first(sel):last(sel)+1)
    c = bincounts(h)[sel]
    s2 = sumw2(h)[sel]
    Hist1D(; binedges = edges, bincounts = c, sumw2 = s2, nentries = nentries(h), overflow=h.overflow)
end
restrict(low, high) = h::Hist1D->restrict(h, low, high)
