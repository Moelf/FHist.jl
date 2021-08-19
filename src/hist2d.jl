struct Hist2D{T<:Real,E} <: AbstractHistogram{T,2,E}
    hist::Histogram{T,2,E}
    sumw2::Array{Float64, 2}
    hlock::SpinLock
    # most concrete inner constructor
    function Hist2D(h::Histogram{T,2,E}, sw2 = copy(h.weights)) where {T,E}
        return new{T,E}(h, sw2, SpinLock())
    end
end

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
    Hist1D(elT::Type{T}=Float64; binedges) where {T}

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
    normalize(h::Hist2D)

Create a normalized histogram via division by `integral(h)`.
"""
function normalize(h::Hist2D)
    return h*(1/integral(h))
end


function Base.show(io::IO, h::Hist2D)
    println(io)
    println(io, "edges: ", h.hist.edges)
    println(io, "bin counts: ", bincounts(h))
    print(io, "total count: ", integral(h))
end

