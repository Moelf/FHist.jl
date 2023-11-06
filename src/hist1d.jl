Base.lock(h::Hist1D) = lock(h.hlock)
Base.unlock(h::Hist1D) = unlock(h.hlock)

"""
    sample(h::Hist1D, n::Int=1)

Sample a histogram's with weights equal to bin count, `n` times.
The sampled values are the bins' lower edges.
"""
function sample(h::Hist1D; n::Int=1)
    StatsBase.sample(binedges(h)[1:end-1], Weights(bincounts(h)), n)
end

@inline function _edge_binindex(r::AbstractRange, x::Real)
    return floor(Int, (x - first(r)) * inv(step(r))) + 1
    # # 20% faster and assigns -Inf, Inf, NaN to typemin(Int64)
    # return Base.unsafe_trunc(Int, round((x - first(r)) * inv(step(r)), RoundDown)) + 1
end

@inline function _edge_binindex(r::AbstractRange{<:Integer}, x::Integer)
    return (x - first(r))÷step(r) + 1
end

@inline function _edge_binindex(v::AbstractVector, x::Real)
    return searchsortedlast(v, x)
end

"""
    bincounts(h::Hist1D)

Get the bin counts of a histogram.
"""
@inline bincounts(h::Hist1D) = h.hist.weights

"""
    binedges(h::Hist1D)

Get the bin edges of a histogram.
"""
@inline binedges(h::Hist1D) = h.hist.edges[1]

"""
    bincenters(h::Hist1D)

Get the bin centers of a histogram.
"""
function bincenters(h::Hist1D)
    StatsBase.midpoints(binedges(h))
end

"""
    binerrors(f::T, h::Hist1D) where T<:Function = f.(h.sumw2)
    binerrors(h::Hist1D) = binerrors(sqrt, h)

Calculate the bin errors from `sumw2` with a Gaussian default.
"""
binerrors(f::T, h::Hist1D) where T<:Function = f.(h.sumw2)
binerrors(h::Hist1D) = binerrors(sqrt, h)


"""
    nbins(h::Hist1D)

Get the number of bins of a histogram.
"""
function nbins(h::Hist1D)
    length(bincounts(h))
end

"""
    integral(h::Hist1D; width=false)

Get the integral a histogram; `width` means multiply each bincount
by their bin width when calculating the integral.

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
    empty!(h::Hist1D)

Resets a histogram's bin counts and `sumw2`.
"""
function Base.empty!(h::Hist1D{T,E}) where {T,E}
    h.hist.weights .= zero(T)
    h.sumw2 .= 0.0
    return h
end

"""
    push!(h::Hist1D, val::Real, wgt::Real=1)
    atomic_push!(h::Hist1D, val::Real, wgt::Real=1)

Adding one value at a time into histogram.
`sumw2` (sum of weights^2) accumulates `wgt^2` with a default weight of 1.
`atomic_push!` is a slower version of `push!` that is thread-safe.

N.B. To append multiple values at once, use broadcasting via
`push!.(h, [-3.0, -2.9, -2.8])` or `push!.(h, [-3.0, -2.9, -2.8], 2.0)`
"""
@inline function atomic_push!(h::Hist1D{T,E}, val::Real, wgt::Real=1) where {T,E}
    lock(h)
    push!(h, val, wgt)
    unlock(h)
    return nothing
end

@inline function Base.push!(h::Hist1D{T,E}, val::Real, wgt::Real=1) where {T,E}
    r = binedges(h)
    L = nbins(h)
    binidx = _edge_binindex(r, val)
    if h.overflow
        binidx = clamp(binidx, 1, L)
        h.nentries[] += 1
        @inbounds h.hist.weights[binidx] += wgt
        @inbounds h.sumw2[binidx] += wgt^2
    else
        if unsigned(binidx - 1) < L
            h.nentries[] += 1
            @inbounds h.hist.weights[binidx] += wgt
            @inbounds h.sumw2[binidx] += wgt^2
        end
    end
    return nothing
end

Base.broadcastable(h::Hist1D) = Ref(h)

"""
    Hist1D(elT::Type{T}=Float64; bins, overflow) where {T}

Initialize an empty histogram with bin content typed as `T` and bin edges.
To be used with [`push!`](@ref). Default overflow behavior (`false`)
will exclude values that are outside of `binedges`.
"""
function Hist1D(elT::Type{T}=Float64; bins, overflow=_default_overflow) where {T}
    counts = zeros(elT, length(bins) - 1)
    return Hist1D(Histogram(bins, counts); overflow=overflow)
end

"""
    Hist1D(array, edges::AbstractRange; overflow)
    Hist1D(array, edges::AbstractVector; overflow)

Create a `Hist1D` with given bin `edges` and vlaues from
array. Weight for each value is assumed to be 1.
"""
function Hist1D(A, r::AbstractRange; overflow=_default_overflow)
    h = Hist1D(Int; bins=r, overflow=overflow)
    firstr = first(r)
    invstep = inv(step(r))
    L = nbins(h)
    nentries = 0
    for val in A
        cursor = floor(Int, (val - firstr) * invstep)
        binidx = cursor + 1
        if overflow
            binidx = clamp(binidx, 1, L)
            nentries += 1
            @inbounds h.hist.weights[binidx] += 1
        else
            if unsigned(cursor) < L
                nentries += 1
                @inbounds h.hist.weights[binidx] += 1
            end
        end
    end
    h.nentries[] = nentries
    h.sumw2 .= h.hist.weights
    return h
end
function Hist1D(A, edges::AbstractVector; overflow=_default_overflow)
    if _is_uniform_bins(edges)
        r = range(first(edges), last(edges), length=length(edges))
        return Hist1D(A, r; overflow=overflow)
    else
        h = Hist1D(Int; bins=edges, overflow=overflow)
        for x in A
            push!(h, x)
        end
        return h
    end
end

"""
    Hist1D(array, wgts::AbstractWeights, edges::AbstractRange; overflow)
    Hist1D(array, wgts::AbstractWeights, edges::AbstractVector; overflow)

Create a `Hist1D` with given bin `edges` and vlaues from
array. `wgts` should have the same `size` as `array`.
"""
function Hist1D(A, wgts::AbstractWeights, r::AbstractRange; overflow=_default_overflow)
    @boundscheck @assert size(A) == size(wgts)
    h = Hist1D(eltype(wgts); bins=r, overflow=overflow)
    push!.(h, A, wgts)
    return h
end
function Hist1D(A, wgts::AbstractWeights, edges::AbstractVector, overflow=_default_overflow)
    @inbounds if _is_uniform_bins(edges)
        r = range(first(edges), last(edges), length=length(edges))
        return Hist1D(A, wgts, r; overflow=overflow)
    else
        h = Hist1D(eltype(wgts); bins=edges, overflow=overflow)
        push!.(h, A, wgts)
        return h
    end
end

"""
    Hist1D(A::AbstractVector{T}; nbins::Integer=_sturges(A), overflow) where T
    Hist1D(A::AbstractVector{T}, wgts::AbstractWeights; nbins::Integer=_sturges(A), overflow) where T

Automatically determine number of bins based on `Sturges` algo.
"""
function Hist1D(A::AbstractVector{T}; nbins::Integer=_sturges(A), overflow=_default_overflow) where {T}
    F = float(T)
    lo, hi = minimum(A), maximum(A)
    r = StatsBase.histrange(F(lo), F(hi), nbins)
    return Hist1D(A, r; overflow=overflow)
end

function Hist1D(
    A::AbstractVector{T},
    wgts::AbstractWeights;
    nbins::Integer=_sturges(A),
    overflow=_default_overflow,
) where {T}
    F = float(T)
    lo, hi = minimum(A), maximum(A)
    r = StatsBase.histrange(F(lo), F(hi), nbins)
    return Hist1D(A, wgts, r; overflow=overflow)
end


"""
    Statistics.mean(h)
    Statistics.std(h)
    Statistics.median(h)
    Statistics.quantile(h::Hist1D, p)

Compute statistical quantities based on the bin centers weighted
by the bin counts.

When the histogram is `Hist2D`, return tuple instead, e.g `(mean(project(h, :x)), mean(project(h, :y)))` etc.
"""
Statistics.mean(h::Hist1D) = Statistics.mean(bincenters(h), Weights(bincounts(h)))
Statistics.std(h::Hist1D) = Statistics.var(bincenters(h), Weights(bincounts(h)); corrected=false) |> sqrt
Statistics.median(h::Hist1D) = Statistics.median(bincenters(h), Weights(bincounts(h)))
Statistics.quantile(h::Hist1D, p) = Statistics.quantile(bincenters(h), Weights(bincounts(h)), p)

"""
    lookup(h::Hist1D, x)

For given x-axis value `x`, find the corresponding bin and return the bin content.
If a value is out of the histogram range, return `missing`.
"""
function lookup(h::Hist1D, x)
    r = binedges(h)
    !(first(r) <= x <= last(r)) && return missing
    return bincounts(h)[_edge_binindex(r, x)]
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
    return h*(1/integral(h; width=width))
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
    h.hist.weights .= f(cumsum(f(h.hist.weights)))
    h.sumw2 .= f(cumsum(f(h.sumw2)))
    return h
end


"""
    rebin(h::Hist1D, n::Int=1)
    rebin(n::Int) = h::Hist1D -> rebin(h, n)

Merges `n` consecutive bins into one.
The returned histogram will have `nbins(h)/n` bins.
"""
function rebin(h::Hist1D, n::Int=1)
    @assert nbins(h) % n == 0
    p = x->Iterators.partition(x, n)
    counts = sum.(p(bincounts(h)))
    sumw2 = sum.(p(h.sumw2))
    edges = first.(p(binedges(h)))
    if _is_uniform_bins(edges)
        edges = range(first(edges), last(edges), length=length(edges))
    end
    return Hist1D(Histogram(edges, counts), sumw2, nentries(h); overflow=h.overflow)
end
rebin(n::Int) = h::Hist1D -> rebin(h, n)

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
    sel = low .<= bincenters(h) .<= high
    @assert sum(sel) > 0 "No bin centers contained in [$(low), $(high)]"
    edgesel = push!(copy(sel), false)

    # include the right edge of the rightmost selected bin
    lastidx = findlast(edgesel)
    if lastidx != nothing
        edgesel[lastidx+1] = 1
    end

    c = bincounts(h)[sel]
    edges = binedges(h)[edgesel]
    sumw2 = h.sumw2[sel]
    if _is_uniform_bins(edges)
        edges = range(first(edges), last(edges), length=length(edges))
    end
    Hist1D(Histogram(edges, c), sumw2, nentries(h); overflow=h.overflow)
end
restrict(low, high) = h::Hist1D->restrict(h, low, high)
