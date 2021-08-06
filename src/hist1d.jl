struct Hist1D{T<:Real,E} <: AbstractHistogram{T,1,E}
    hist::Histogram{T,1,E}
    sumw2::Vector{Float64}
    hlock::SpinLock
    # most concrete inner constructor
    function Hist1D(h::Histogram{T,1,E}, sw2 = copy(h.weights)) where {T,E}
        return new{T,E}(h, sw2, SpinLock())
    end
end
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
    return floor(Int, (x - first(r)) / step(r)) + 1
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
    nbins(h::Hist1D)

Get the number of bins of a histogram.
"""
function nbins(h::Hist1D)
    length(bincounts(h))
end

"""
    integral(h::Hist1D)

Get the integral a histogram.
"""
function integral(h::Hist1D)
    sum(bincounts(h))
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
    unsafe_push!(h::Hist1D, val::Real, wgt::Real=1)
    push!(h::Hist1D, val::Real, wgt::Real=1)

Adding one value at a time into histogram. 
`sumw2` (sum of weights^2) accumulates `wgt^2` with a default weight of 1.
`unsafe_push!` is a faster version of `push!` that is not thread-safe.

N.B. To append multiple values at once, use broadcasting via
`push!.(h, [-3.0, -2.9, -2.8])` or `push!.(h, [-3.0, -2.9, -2.8], 2.0)`
"""
@inline function Base.push!(h::Hist1D{T,E}, val::Real, wgt::Real=1) where {T,E}
    lock(h)
    unsafe_push!(h, val, wgt)
    unlock(h)
    return nothing
end

@inline function unsafe_push!(h::Hist1D{T,E}, val::Real, wgt::Real=1) where {T,E}
    r = @inbounds h.hist.edges[1]
    L = length(r) - 1
    start = first(r)
    stop = last(r)
    c = ifelse(val > stop, 0, 1)
    c = ifelse(val < start, 0, c)
    binidx = _edge_binindex(r, val)
    binidx = ifelse(binidx > L, L, binidx)
    binidx = ifelse(binidx < 1, 1, binidx)
    @inbounds h.hist.weights[binidx] += c*wgt
    @inbounds h.sumw2[binidx] += c*wgt^2
    return nothing
end

Base.broadcastable(h::Hist1D) = Ref(h)

"""
    Hist1D(elT::Type{T}=Float64; binedges) where {T}

Initialize an empty histogram with bin content typed as `T` and bin edges.
To be used with [`push!`](@ref)
"""
function Hist1D(elT::Type{T}=Float64; bins) where {T}
    counts = zeros(elT, length(bins) - 1)
    return Hist1D(Histogram(bins, counts))
end

"""
    Hist1D(array, edges::AbstractRange)
    Hist1D(array, edges::AbstractVector)

Create a `Hist1D` with given bin `edges` and vlaues from
array. Weight for each value is assumed to be 1.
"""
function Hist1D(A::AbstractVector, r::AbstractRange)
    h = Hist1D(Int; bins=r)
    unsafe_push!.(h, A)
    return h
end
function Hist1D(A::AbstractVector, edges::AbstractVector)
    if _is_uniform_bins(edges)
        s = edges[2] - first(edges)
        r = first(edges):s:last(edges)
        return Hist1D(A, r)
    else
        h = Hist1D(Int; bins=edges)
        unsafe_push!.(h, A)
        return h
    end
end

"""
    Hist1D(array, wgts::AbstractWeights, edges::AbstractRange)
    Hist1D(array, wgts::AbstractWeights, edges::AbstractVector)

Create a `Hist1D` with given bin `edges` and vlaues from
array. `wgts` should have the same `size` as `array`.
"""
function Hist1D(A, wgts::AbstractWeights, r::AbstractRange)
    @boundscheck @assert size(A) == size(wgts)
    h = Hist1D(eltype(wgts); bins=r)
    unsafe_push!.(h, A, wgts)
    return h
end
function Hist1D(A, wgts::AbstractWeights, edges::AbstractVector)
    @inbounds if _is_uniform_bins(edges)
        s = edges[2] - first(edges)
        r = first(edges):s:last(edges)
        return Hist1D(A, wgts, r)
    else
        h = Hist1D(eltype(wgts); bins=edges)
        unsafe_push!.(h, A, wgts)
        return h
    end
end

"""
    Hist1D(A::AbstractVector{T}; nbins::Integer=StatsBase.sturges(length(A))) where T
    Hist1D(A::AbstractVector{T}, wgts::AbstractWeights; nbins::Integer=StatsBase.sturges(length(A))) where T

Automatically determine number of bins based on `Sturges` algo.
"""
function Hist1D(A::AbstractVector{T}; nbins::Integer=StatsBase.sturges(length(A))) where {T}
    F = float(T)
    lo, hi = extrema(A)
    r = StatsBase.histrange(F(lo), F(hi), nbins)
    return Hist1D(A, r)
end

function Hist1D(
    A::AbstractVector{T},
    wgts::AbstractWeights;
    nbins::Integer=StatsBase.sturges(length(A)),
) where {T}
    F = float(T)
    lo, hi = extrema(A)
    r = StatsBase.histrange(F(lo), F(hi), nbins)
    return Hist1D(A, wgts, r)
end


"""
    Statistics.mean(h::Hist1D)
    Statistics.std(h::Hist1D)
    Statistics.median(h::Hist1D)
    Statistics.quantile(h::Hist1D, p)

Compute statistical quantities based on the bin centers weighted
by the bin counts.
"""
Statistics.mean(h::Hist1D) = Statistics.mean(bincenters(h), Weights(bincounts(h)))
Statistics.std(h::Hist1D) = sqrt(Statistics.var(bincenters(h), Weights(bincounts(h))))
Statistics.median(h::Hist1D) = Statistics.median(bincenters(h), Weights(bincounts(h)))
Statistics.quantile(h::Hist1D, p) = Statistics.quantile(bincenters(h), Weights(bincounts(h)), p)

"""
    function lookup(h::Hist1D, v) 

For given x-axis value`v`, find the corresponding bin and return the bin content.
If a value is out of the histogram range, return `missing`.
"""
function lookup(h::Hist1D, v)
    r = binedges(h)
    !(first(r) <= v <= last(r)) && return missing
    binidx = searchsortedlast(r, v) # TODO replace with `_edges_binindex`
    return bincounts(h)[binidx]
end

"""
    normalize(h::Hist1D)

Create a normalized histogram via division by `integral(h)`.
"""
function normalize(h::Hist1D)
    return h*(1/integral(h))
end

"""
    cumulative(h::Hist1D; forward=true)::Hist1D

Create a cumulative histogram. If `forward`, start
summing from the left.
"""
function cumulative(h::Hist1D; forward=true)::Hist1D
    # https://root.cern.ch/doc/master/TH1_8cxx_source.html#l02608
    f = forward ? identity : reverse
    h = deepcopy(h)
    h.hist.weights .= f(cumsum(h.hist.weights))
    h.sumw2 .= f(cumsum(h.sumw2))
    return h
end

function Base.show(io::IO, h::Hist1D)
    _e = binedges(h)
    if _e isa AbstractRange && length(_e) < 50
        _h = Histogram(float(_e), bincounts(h))
        show(io, UnicodePlots.histogram(_h; width=30, xlabel=""))
    end
    println()
    println(io, "edges: ", binedges(h))
    println(io, "bin counts: ", bincounts(h))
    print(io, "total count: ", sum(bincounts(h)))
end
