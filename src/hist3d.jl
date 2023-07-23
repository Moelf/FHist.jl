"""
    Hist3D(elT::Type{T}, binedges; overflow) where {T}

Initialize an empty histogram with bin content typed as `T` and bin edges.
To be used with [`push!`](@ref). Default overflow behavior (`false`)
will exclude values that are outside of `binedges`.
"""
@depkws function Hist3D(elT::Type{T}=Float64; binedges=(-1:1, -1:1, -1:1), overflow=_default_overflow,
    @deprecate(bins, binedges)) where {T}
    counts = zeros(elT, length.(binedges) .- 1)
    return Hist3D(Histogram(binedges, counts); overflow=overflow)
end

"""
    Hist3D(tuple, binedges::NTuple{3,AbstractRange}; overflow)
    Hist3D(tuple, binedges::NTuple{3,AbstractVector}; overflow)

Create a `Hist3D` with given bin `binedges` and values from
a 2-tuple of arrays of x, y values. Weight for each value is assumed to be 1.
"""
function Hist3D(A::NTuple{3,AbstractVector}, binedges::NTuple{3,AbstractRange}; overflow=_default_overflow)
    h = Hist3D(Int; binedges=binedges, overflow=overflow)
    push!.(h, A[1], A[2], A[3])
    return h
end

function Hist3D(A::NTuple{3,AbstractVector}, binedges::NTuple{3,AbstractVector}; overflow=_default_overflow)
    if all(_is_uniform_bins.(binedges))
        r = (range(first(binedges[1]), last(binedges[1]), length=length(binedges[1])),
            range(first(binedges[2]), last(binedges[2]), length=length(binedges[2])),
            range(first(binedges[3]), last(binedges[3]), length=length(binedges[3])))
        return Hist3D(A, r; overflow=overflow)
    else
        h = Hist3D(Int; binedges=binedges, overflow=overflow)
        push!.(h, A[1], A[2], A[3])
        return h
    end
end

"""
    Hist3D(tuple, wgts::AbstractWeights, binedges::NTuple{3,AbstractRange}; overflow)
    Hist3D(tuple, wgts::AbstractWeights, binedges::NTuple{3,AbstractVector}; overflow)

Create a `Hist3D` with given bin `binedges` and values from
a 2-tuple of arrays of x, y values.
`wgts` should have the same `size` as elements of `tuple`.
"""
function Hist3D(A::NTuple{3,AbstractVector}, wgts::AbstractWeights, binedges::NTuple{3,AbstractRange}; overflow=_default_overflow)
    @boundscheck @assert size(A[1]) == size(A[2]) == size(A[3]) == size(wgts)
    h = Hist3D(eltype(wgts); binedges=binedges, overflow=overflow)
    push!.(h, A[1], A[2], A[3], wgts)
    return h
end

function Hist3D(A::NTuple{3,AbstractVector}, wgts::AbstractWeights, binedges::NTuple{3,AbstractVector}; overflow=_default_overflow)
    if all(_is_uniform_bins.(binedges))
        r = (range(first(binedges[1]), last(binedges[1]), length=length(binedges[1])),
            range(first(binedges[2]), last(binedges[2]), length=length(binedges[2])),
            range(first(binedges[3]), last(binedges[3]), length=length(binedges[3])))
        return Hist3D(A, wgts, r; overflow=overflow)
    else
        h = Hist3D(Int; binedges=binedges, overflow=overflow)
        push!.(h, A[1], A[2], A[3], wgts)
        return h
    end
end

"""
    Hist3D(A::AbstractVector{T}; nbins::NTuple{3,Integer}, overflow) where T
    Hist3D(A::AbstractVector{T}, wgts::AbstractWeights; nbins::NTuple{3,Integer}, overflow) where T

Automatically determine number of bins based on `Sturges` algo.
"""
function Hist3D(A::NTuple{3,AbstractVector{T}};
    nbins::NTuple{3,Integer}=_sturges.(A),
    overflow=_default_overflow
) where {T}
    F = float(T)
    nbinsx, nbinsy, nbinsz = nbins
    lox, hix = minimum(A[1]), maximum(A[1])
    loy, hiy = minimum(A[2]), maximum(A[2])
    loz, hiz = minimum(A[3]), maximum(A[3])
    rx = StatsBase.histrange(F(lox), F(hix), nbinsx)
    ry = StatsBase.histrange(F(loy), F(hiy), nbinsy)
    rz = StatsBase.histrange(F(loz), F(hiz), nbinsz)
    r = (rx, ry, rz)
    return Hist3D(A, r; overflow=overflow)
end

function Hist3D(A::NTuple{3,AbstractVector{T}}, wgts::AbstractWeights;
    nbins::NTuple{3,Integer}=_sturges.(A),
    overflow=_default_overflow
) where {T}
    F = float(T)
    nbinsx, nbinsy, nbinsz = nbins
    lox, hix = minimum(A[1]), maximum(A[1])
    loy, hiy = minimum(A[2]), maximum(A[2])
    loz, hiz = minimum(A[3]), maximum(A[3])
    rx = StatsBase.histrange(F(lox), F(hix), nbinsx)
    ry = StatsBase.histrange(F(loy), F(hiy), nbinsy)
    rz = StatsBase.histrange(F(loz), F(hiz), nbinsz)
    r = (rx, ry, rz)
    return Hist3D(A, wgts, r; overflow=overflow)
end

Base.lock(h::Hist3D) = lock(h.hlock)
Base.unlock(h::Hist3D) = unlock(h.hlock)

"""
    bincounts(h::Hist3D)

Get the bin counts of a histogram.
"""
@inline bincounts(h::Hist3D) = h.hist.weights

"""
    binedges(h::Hist3D)

Get a 3-tuple of the bin edges of a histogram.
"""
@inline binedges(h::Hist3D) = h.hist.edges

"""
    bincenters(h::Hist3D)

Get a 2-tuple of the bin centers of a histogram.
"""
function bincenters(h::Hist3D)
    StatsBase.midpoints.(binedges(h))
end

"""
    binerrors(f::T, h::Hist3D) where T<:Function = f.(h.sumw2)
    binerrors(h::Hist3D) = binerrors(sqrt, h)

Calculate the bin errors from `sumw2` with a Gaussian default.
"""
binerrors(f::T, h::Hist3D) where {T<:Function} = f.(h.sumw2)
binerrors(h::Hist3D) = binerrors(sqrt, h)

"""
    nbins(h::Hist3D)

Get a 3-tuple of the number of x and y bins of a histogram.
"""
function nbins(h::Hist3D)
    size(bincounts(h))
end

"""
    integral(h::Hist3D)

Get the integral a histogram.
"""
function integral(h::Hist3D)
    sum(bincounts(h))
end

"""
    empty!(h::Hist3D)

Resets a histogram's bin counts and `sumw2`.
"""
function Base.empty!(h::Hist3D{T,E}) where {T,E}
    h.hist.weights .= zero(T)
    h.sumw2 .= 0.0
    return h
end

"""
    push!(h::Hist3D, valx::Real, valy::Real, wgt::Real=1)
    atomic_push!(h::Hist3D, valx::Real, valy::Real, wgt::Real=1)

Adding one value at a time into histogram.
`sumw2` (sum of weights^2) accumulates `wgt^2` with a default weight of 1.
`atomic_push!` is a slower version of `push!` that is thread-safe.
"""
@inline function Base.push!(h::Hist3D{T,E}, valx::Real, valy::Real, valz::Real, wgt::Real=1) where {T,E}
    rx, ry, rz = binedges(h)
    Lx, Ly, Lz = nbins(h)
    binidxx = _edge_binindex(rx, valx)
    binidxy = _edge_binindex(ry, valy)
    binidxz = _edge_binindex(rz, valz)
    h.nentries[] += 1
    if h.overflow
        binidxx = clamp(binidxx, 1, Lx)
        binidxy = clamp(binidxy, 1, Ly)
        binidxz = clamp(binidxz, 1, Lz)
        @inbounds h.hist.weights[binidxx, binidxy, binidxz] += wgt
        @inbounds h.sumw2[binidxx, binidxy, binidxz] += wgt^2
    else
        if (unsigned(binidxx - 1) < Lx) && (unsigned(binidxy - 1) < Ly) && (unsigned(binidxz - 1) < Lz)
            @inbounds h.hist.weights[binidxx, binidxy, binidxz] += wgt
            @inbounds h.sumw2[binidxx, binidxy, binidxz] += wgt^2
        end
    end
    return nothing
end

@inline function atomic_push!(h::Hist3D{T,E}, valx::Real, valy::Real, valz::Real, wgt::Real=1) where {T,E}
    lock(h)
    push!(h, valx, valy, valz, wgt)
    unlock(h)
    return nothing
end

Base.broadcastable(h::Hist3D) = Ref(h)

"""
    function lookup(h::Hist3D, x, y)

For given x/y/z-axis value `x`, `y`, `z`, find the corresponding bin and return the bin content.
If a value is out of the histogram range, return `missing`.
"""
function lookup(h::Hist3D, x, y, z)
    rx, ry, rz = binedges(h)
    !(first(rx) <= x <= last(rx)) && return missing
    !(first(ry) <= y <= last(ry)) && return missing
    !(first(rz) <= z <= last(rz)) && return missing
    return bincounts(h)[_edge_binindex(rx, x), _edge_binindex(ry, y), _edge_binindex(rz, z)]
end

"""
    normalize(h::Hist3D)

Create a normalized histogram via division by `integral(h)`.
"""
function normalize(h::Hist3D)
    return h * (1 / integral(h))
end

"""
    project(h::Hist3D, axis::Symbol=:x)
    project(axis::Symbol=:x) = h::Hist3D -> project(h, axis)

Computes the `:x`/`:y`/`:z` axis projection of the 3D histogram by
summing over the specified axis. Returns a `Hist2D`.
"""
function project(h::Hist3D, axis::Symbol=:x)
    @assert axis ∈ (:x, :y, :z)
    dimremove, dimskeep = if axis == :z
        3, [1, 2]
    elseif axis == :y
        2, [1, 3]
    else
        1, [2, 3]
    end
    counts = dropdims(sum(bincounts(h), dims=dimremove), dims=dimremove)
    sumw2 = dropdims(sum(h.sumw2, dims=dimremove), dims=dimremove)
    edges = binedges(h)[dimskeep]
    return Hist2D(Histogram(edges, counts), sumw2; overflow=h.overflow)
end

project(axis::Symbol=:x) = h::Union{Hist2D,Hist3D} -> project(h, axis)