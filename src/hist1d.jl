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

@inline function _edge_binindex(r::AbstractRange{<:Integer}, x::Integer)
    return x - first(r) + 1
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
    binidx = _edge_binindex(r, val)
    if unsigned(binidx - 1) < L
        @inbounds h.hist.weights[binidx] += wgt
        @inbounds h.sumw2[binidx] += wgt^2
    end
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
        r = range(first(edges), last(edges), length=length(edges))
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
        r = range(first(edges), last(edges), length=length(edges))
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
    lo, hi = minimum(A), maximum(A)
    r = StatsBase.histrange(F(lo), F(hi), nbins)
    return Hist1D(A, r)
end

function Hist1D(
    A::AbstractVector{T},
    wgts::AbstractWeights;
    nbins::Integer=StatsBase.sturges(length(A)),
) where {T}
    F = float(T)
    lo, hi = minimum(A), maximum(A)
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
    function lookup(h::Hist1D, x)

For given x-axis value `x`, find the corresponding bin and return the bin content.
If a value is out of the histogram range, return `missing`.
"""
function lookup(h::Hist1D, x)
    r = binedges(h)
    !(first(r) <= x <= last(r)) && return missing
    return bincounts(h)[_edge_binindex(r, x)]
end

"""
    normalize(h::Hist1D)

Create a normalized histogram via division by `integral(h)`.
"""
function normalize(h::Hist1D)
    return h*(1/integral(h))
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
    return Hist1D(Histogram(edges, counts), sumw2)
end
rebin(n::Int) = h::Hist1D -> rebin(h, n)

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
    Hist1D(Histogram(edges, c), sumw2)
end
restrict(low=-Inf, high=Inf) = h::Hist1D->restrict(h, low, high)


function _svg(h::Hist1D)
    paddingx, paddingy = 0.05, 0.10
    framewidth, frameheight = 250, 200
    fill = "#ffffff00" # 0 alpha
    strokecolor, strokewidth = "black", 1
    strokewidth = 1
    _c, _e = bincounts(h), binedges(h)
    ys = frameheight * ((2*paddingy-1)/maximum(_c) .* _c .+ (1-paddingy))
    xs = framewidth * (
        (1 - 2 * paddingx) / (maximum(_e) - minimum(_e))
        .* (_e .- minimum(_e))
        .+ paddingx
    )
    points = [(paddingx*framewidth, (1-paddingy)*frameheight)] # bottom left
    for i in 1:(length(xs)-1)
        push!(points, (xs[i],ys[i])) # left bin edge
        push!(points, (xs[i+1],ys[i])) # right bin edge
    end
    push!(points, ((1-paddingx)*framewidth,(1-paddingy)*frameheight))
    push!(points, points[1]) # close path
    pathstr = join(["$(x),$(y)" for (x,y) in points],",")
    return """
    <svg width="$(framewidth)" height="$(frameheight)" version="1.1" xmlns="http://www.w3.org/2000/svg">
        <polyline points="$(pathstr)" stroke="$(strokecolor)" fill="$(fill)" stroke-width="$(strokewidth)"/>
        <polyline points="$(framewidth*paddingx),$(frameheight*(1-paddingy)),$(framewidth*(1-paddingx)),$(frameheight*(1-paddingy))" stroke="black" stroke-width="1"/>
        <text x="$(framewidth*paddingx)" y="$(frameheight*(1-0.5*paddingy))" dominant-baseline="middle" text-anchor="start" fill="black">$(minimum(_e))</text>
        <text x="$(framewidth*(1-paddingx))" y="$(frameheight*(1-0.5*paddingy))" dominant-baseline="middle" text-anchor="end" fill="black">$(maximum(_e))</text>
    </svg>
    """
end

function Base.show(io::IO, h::Hist1D)
    if (nbins(h) < 50) && all(bincounts(h) .>= 0)
        _e = binedges(h)
        _h = Histogram(float.(_e), bincounts(h))
        show(io, UnicodePlots.histogram(_h; width=30, xlabel=""))
    end
    println(io)
    println(io, "edges: ", binedges(h))
    println(io, "bin counts: ", bincounts(h))
    print(io, "total count: ", integral(h))
end

function Base.show(io::IO, m::MIME"text/html", h::Hist1D)
    println(io, """
    <div style="display: flex;">
        <div style="float:left; margin:5px">$(_svg(h))</div>
        <div style="float:left; margin:5px; max-width: 50%; display:flex; justify-content:center; align-items:center;">
            <ul>
                <li>edges: $(repr(binedges(h), context=:limit => true))</li>
                <li>bin counts: $(repr(bincounts(h), context=:limit => true))</li>
                <li>maximum count: $(maximum(bincounts(h)))</li>
                <li>total count: $(integral(h))</li>
            </ul>
        </div>
    </div>
    """)
end
