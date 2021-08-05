struct Hist1D{T<:Real,E} <: AbstractHistogram{T,1,E}
    hist::Histogram{T,1,E}
    sumw2::Vector{Float64}
    hlock::SpinLock
    # most concrete inner constructor
    function Hist1D(h::Histogram{T,1,E}, sw2 = h.weights) where {T,E}
        return new{T,E}(h, sw2, SpinLock())
    end
end
Base.lock(h::Hist1D) = lock(h.hlock)
Base.unlock(h::Hist1D) = unlock(h.hlock)

"""
    sample(h::Hist1D)
    sample(h::Hist1D, n::Int)

Sample a histogram's with weights equal to bin count, one or `n` times.
The returned sample value will be one of the bin's left edge.
"""
function sample(h::Hist1D)
    @inbounds sample(@view(only(h.hist.edges)[1:end-1]), Weights(h.hist.weights))
end
function sample(h::Hist1D, n::Int)
    @inbounds sample(@view(only(h.hist.edges)[1:end-1]), Weights(h.hist.weights), n)
end


@inline function _edge_binindex(r::AbstractRange{T}, x::Real) where T <:  AbstractFloat
    s = step(r)
    start = first(r) + 0.5s
    return round(Int, (x - start) / s) + 1
end

@inline function _edge_binindex(r::AbstractRange{T}, x::Real) where T <:  Integer
    s = step(r)
    start = first(r)
    return Int(fld(x-start, s)) + 1
end

@inline function _edge_binindex(v::AbstractVector, x::Real)
    return searchsortedlast(v, x)
end

"""
    bincenters(h::Hist1D)

Get the bin centers of a histogram.
"""
function bincenters(h::Hist1D)
    edges = h.hist.edges[1]
    edges[2:end] - diff(edges) ./ 2
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
function Hist1D(A::AbstractVector, r::AbstractRange{T}) where T <: Union{AbstractFloat,Integer}
    h = Hist1D(Int; bins=r)
    unsafe_push!.(Ref(h), A)
    return h
end
function Hist1D(A::AbstractVector, edges::AbstractVector)
    if _is_uniform_bins(edges)
        s = edges[2] - first(edges)
        r = first(edges):s:last(edges)
        Hist1D(A, r)
    else
        Hist1D(fit(Histogram, A, edges))
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
    s = step(r)
    start = first(r)
    start2 = start + s / 2
    stop = last(r)
    wgt_zero = zero(eltype(wgts))
    L = length(r) - 1
    counts = zeros(L)
    sumw2 = zeros(L)
    @inbounds for i in eachindex(A)
        # skip overflow
        c = ifelse(A[i] < start || A[i] > stop, wgt_zero, wgts[i])
        id = round(Int, (A[i] - start2) / s) + 1
        idx = clamp(id, 1, L)
        counts[idx] += c
        sumw2[idx] += c^2
    end
    return Hist1D(Histogram(r, counts), sumw2)
end
function Hist1D(A, wgts::AbstractWeights, edges::AbstractVector)
    @inbounds if _is_uniform_bins(edges)
        s = edges[2] - first(edges)
        r = first(edges):s:last(edges)
        Hist1D(A, wgts, r)
    else
        hist = Hist1D(fit(Histogram, A, wgts, edges)).hist
        sw2 = Hist1D(fit(Histogram, A, Weights(wgts.^2), edges)).hist.weights
        Hist1D(hist, sw2)
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

function Base.show(io::IO, h::Hist1D)
    _e = h.hist.edges[1]
    if _e isa AbstractRange && length(_e) < 50
        _h = Histogram(float(_e), h.hist.weights)
        show(io, UnicodePlots.histogram(_h; width=30, xlabel=""))
    end
    println()
    println(io, "edges: ", h.hist.edges[1])
    println(io, "bin counts: ", h.hist.weights)
    print(io, "total count: ", sum(h.hist.weights))
end
