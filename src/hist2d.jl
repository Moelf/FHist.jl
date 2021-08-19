Base.lock(h::Hist2D) = lock(h.hlock)
Base.unlock(h::Hist2D) = unlock(h.hlock)

"""
    bincounts(h::Hist2D)

Get the bin counts of a histogram.
"""
@inline bincounts(h::Hist2D) = h.hist.weights

"""
    binedges(h::Hist2D)

Get the bin edges of a histogram.
"""
@inline binedges(h::Hist2D) = h.hist.edges

"""
    bincenters(h::Hist2D)

Get the bin centers of a histogram.
"""
function bincenters(h::Hist2D)
    StatsBase.midpoints.(binedges(h))
end


"""
    nbins(h::Hist2D)

Get a tuple of the number of x and y bins of a histogram.
"""
function nbins(h::Hist2D)
    size(bincounts(h))
end

"""
    integral(h::Hist2D)

Get the integral a histogram.
"""
function integral(h::Hist2D)
    sum(bincounts(h))
end

"""
    empty!(h::Hist2D)

Resets a histogram's bin counts and `sumw2`.
"""
function Base.empty!(h::Hist2D{T,E}) where {T,E}
    h.hist.weights .= zero(T)
    h.sumw2 .= 0.0
    return h
end

"""
    unsafe_push!(h::Hist2D, valx::Real, valy::Real, wgt::Real=1)
    push!(h::Hist2D, valx::Real, valy::Real, wgt::Real=1)

Adding one value at a time into histogram.
`sumw2` (sum of weights^2) accumulates `wgt^2` with a default weight of 1.
`unsafe_push!` is a faster version of `push!` that is not thread-safe.

"""
@inline function Base.push!(h::Hist2D{T,E}, valx::Real, valy::Real, wgt::Real=1) where {T,E}
    lock(h)
    unsafe_push!(h, valx, valy, wgt)
    unlock(h)
    return nothing
end

@inline function unsafe_push!(h::Hist2D{T,E}, valx::Real, valy::Real, wgt::Real=1) where {T,E}
    rx = @inbounds h.hist.edges[1]
    ry = @inbounds h.hist.edges[2]
    Lx = length(rx) - 1
    Ly = length(ry) - 1
    binidxx = _edge_binindex(rx, valx)
    binidxy = _edge_binindex(ry, valy)
    if (1 <= binidxx <= Lx) && (1 <= binidxy <= Ly)
        @inbounds h.hist.weights[binidxx,binidxy] += wgt
        @inbounds h.sumw2[binidxx,binidxy] += wgt^2
    end
    return nothing
end

Base.broadcastable(h::Hist2D) = Ref(h)

"""
    Hist2D(elT::Type{T}=Float64; binedges) where {T}

Initialize an empty histogram with bin content typed as `T` and bin edges.
To be used with [`push!`](@ref)
"""
function Hist2D(elT::Type{T}=Float64; bins) where {T}
    counts = zeros(elT, (length(bins[1]) - 1, length(bins[2])-1))
    return Hist2D(Histogram(bins, counts))
end

"""
    Hist2D(tuple, edges::Tuple{AbstractRange,AbstractRange})
    Hist2D(tuple, edges::Tuple{AbstractVector,AbstractVector})

Create a `Hist2D` with given bin `edges` and values from
a 2-tuple of arrays of x, y values. Weight for each value is assumed to be 1.
"""
function Hist2D(A::Tuple{AbstractVector,AbstractVector}, r::Tuple{AbstractRange,AbstractRange})
    h = Hist2D(Int; bins=r)
    unsafe_push!.(h, A[1], A[2])
    return h
end
function Hist2D(A::Tuple{AbstractVector,AbstractVector}, edges::Tuple{AbstractVector,AbstractVector})
    if all(_is_uniform_bins.(edges))
        r = (range(first(edges[1]), last(edges[1]), length=length(edges[1])),
             range(first(edges[2]), last(edges[2]), length=length(edges[2])))
        return Hist2D(A, r)
    else
        h = Hist2D(Int; bins=edges)
        unsafe_push!.(h, A[1], A[2])
        return h
    end
end

"""
    Hist2D(tuple, wgts::AbstractWeights, edges::Tuple{AbstractRange,AbstractRange})
    Hist2D(tuple, wgts::AbstractWeights, edges::Tuple{AbstractVector,AbstractVector})

Create a `Hist2D` with given bin `edges` and values from
a 2-tuple of arrays of x, y values.
`wgts` should have the same `size` as elements of `tuple`.
"""
function Hist2D(A::Tuple{AbstractVector,AbstractVector}, wgts::AbstractWeights, r::Tuple{AbstractRange,AbstractRange})
    @boundscheck @assert size(A[1]) == size(A[2]) == size(wgts)
    h = Hist2D(eltype(wgts); bins=r)
    unsafe_push!.(h, A[1], A[2], wgts)
    return h
end
function Hist2D(A::Tuple{AbstractVector,AbstractVector}, wgts::AbstractWeights, r::Tuple{AbstractVector,AbstractVector})
    if all(_is_uniform_bins.(edges))
        r = (range(first(edges[1]), last(edges[1]), length=length(edges[1])),
             range(first(edges[2]), last(edges[2]), length=length(edges[2])))
        return Hist2D(A, wgts, r)
    else
        h = Hist2D(Int; bins=edges)
        unsafe_push!.(h, A[1], A[2], wgts)
        return h
    end
end

"""
    Hist2D(A::AbstractVector{T}; nbins::Tuple{Integer,Integer}) where T
    Hist2D(A::AbstractVector{T}, wgts::AbstractWeights; nbins::Tuple{Integer,Integer}) where T

Automatically determine number of bins based on `Sturges` algo.
"""
function Hist2D(A::Tuple{AbstractVector{T},AbstractVector{T}};
        nbins::Tuple{Integer,Integer}=(StatsBase.sturges(length(A[1])),
                                       StatsBase.sturges(length(A[2])))) where {T}
    F = float(T)
    nbinsx, nbinsy = nbins
    lox, hix = extrema(A[1])
    loy, hiy = extrema(A[2])
    rx = StatsBase.histrange(F(lox), F(hix), nbinsx)
    ry = StatsBase.histrange(F(loy), F(hiy), nbinsy)
    r = (rx, ry)
    return Hist2D(A, r)
end

function Hist2D(A::Tuple{AbstractVector{T},AbstractVector{T}}, wgts::AbstractWeights;
        nbins::Tuple{Integer,Integer}=(StatsBase.sturges(length(A[1])),
                                       StatsBase.sturges(length(A[2])))) where {T}
    F = float(T)
    nbinsx, nbinsy = nbins
    lox, hix = extrema(A[1])
    loy, hiy = extrema(A[2])
    rx = StatsBase.histrange(F(lox), F(hix), nbinsx)
    ry = StatsBase.histrange(F(loy), F(hiy), nbinsy)
    r = (rx, ry)
    return Hist2D(A, wgts, r)
end

"""
    function lookup(h::Hist1D, x)

For given x-axis and y-axis value `x`, `y`, find the corresponding bin and return the bin content.
If a value is out of the histogram range, return `missing`.
"""
function lookup(h::Hist2D, x, y)
    rx, ry = binedges(h)
    !(first(rx) <= x <= last(rx)) && return missing
    !(first(ry) <= y <= last(ry)) && return missing
    return bincounts(h)[_edge_binindex(rx, x), _edge_binindex(ry, y)]
end


"""
    normalize(h::Hist2D)

Create a normalized histogram via division by `integral(h)`.
"""
function normalize(h::Hist2D)
    return h*(1/integral(h))
end

# TODO profile

"""
    project(h::Hist2D, axis::Symbol=:x)
    project(axis::Symbol=:x) = h::Hist2D -> project(h, axis)

Computes the `:x` (`:y`) axis projection of the 2D histogram by
summing over the y (x) axis. Returns a `Hist1D`.
"""
function project(h::Hist2D, axis::Symbol=:x)
    @assert axis ∈ (:x, :y)
    dim = axis == :x ? 2 : 1
    ex, ey = binedges(h)
    counts = [sum(bincounts(h), dims=dim)...]
    sumw2 = [sum(h.sumw2, dims=dim)...]
    edges = axis == :x ? ex : ey
    return Hist1D(Histogram(edges, counts), sumw2)
end
project(axis::Symbol=:x) = h::Hist2D -> project(h, axis)

"""
    rebin(h::Hist2D, nx::Int=1, ny::Int=nx)
    rebin(nx::Int, ny::Int) = h::Hist2D -> rebin(h, nx, ny)

Merges `nx` (`ny`) consecutive bins into one along the x (y) axis by summing.
"""
function rebin(h::Hist2D, nx::Int=1, ny::Int=nx)
    sx, sy = nbins(h)
    @assert sx % nx == sy % ny == 0
    p1d = (x,n)->Iterators.partition(x, n)
    p2d = x->(x[i:i+(nx-1),j:j+(ny-1)] for i=1:nx:sx, j=1:ny:sy)
    counts = sum.(p2d(bincounts(h)))
    sumw2 = sum.(p2d(h.sumw2))
    ex = first.(p1d(binedges(h)[1], nx))
    ey = first.(p1d(binedges(h)[2], ny))
    _is_uniform_bins(ex) && (ex = range(first(ex), last(ex), length=length(ex)))
    _is_uniform_bins(ey) && (ey = range(first(ey), last(ey), length=length(ey)))
    return Hist2D(Histogram((ex,ey), counts), sumw2)
end
rebin(nx::Int, ny::Int) = h::Hist2D -> rebin(h, nx, ny)

function Base.show(io::IO, h::Hist2D)
    println(io)
    println(io, "edges: ", h.hist.edges)
    println(io, "bin counts: ", bincounts(h))
    print(io, "total count: ", integral(h))
end

