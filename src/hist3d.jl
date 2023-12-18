function auto_bins(ary, ::Val{3}; nbins=nothing)
    xs, ys, zs = ary
    E = eltype(xs)
    xnbins, ynbins, znbins = isnothing(nbins) ? _sturges.(xs) : nbins
    F = E <: Number ? float(E) : Float64
    lo, hi = minimum(xs), maximum(xs)
    loy, hiy = minimum(ys), maximum(ys)
    loz, hiz = minimum(zs), maximum(zs)
    (StatsBase.histrange(F(lo), F(hi), xnbins),
        StatsBase.histrange(F(loy), F(hiy), ynbins),
        StatsBase.histrange(F(loz), F(hiz), znbins),)
end

"""
    binerrors(f::T, h::Hist3D) where T<:Function = f.(h.sumw2)
    binerrors(h::Hist3D) = binerrors(sqrt, h)

Calculate the bin errors from `sumw2` with a Gaussian default.
"""
binerrors(f::T, h::Hist3D) where T<:Function = f.(h.sumw2)
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
    push!(h::Hist3D, valx::Real, valy::Real, wgt::Real=1)
    atomic_push!(h::Hist3D, valx::Real, valy::Real, wgt::Real=1)

Adding one value at a time into histogram.
`sumw2` (sum of weights^2) accumulates `wgt^2` with a default weight of 1.
`atomic_push!` is a slower version of `push!` that is thread-safe.

"""
@inline function atomic_push!(h::Hist3D{T,E}, valx::Real, valy::Real, valz::Real, wgt::Real=1) where {T,E}
    lock(h)
    push!(h, valx, valy, valz, wgt)
    unlock(h)
    return nothing
end

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
        @inbounds bincounts(h)[binidxx,binidxy,binidxz] += wgt
        @inbounds sumw2(h)[binidxx,binidxy,binidxz] += wgt^2
    else
        if (unsigned(binidxx - 1) < Lx) && (unsigned(binidxy - 1) < Ly) && (unsigned(binidxz - 1) < Lz)
            @inbounds bincounts(h)[binidxx,binidxy,binidxz] += wgt
            @inbounds sumw2(h)[binidxx,binidxy,binidxz] += wgt^2
        end
    end
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
    return h*(1/integral(h))
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
        3, [1,2]
    elseif axis == :y
        2, [1,3]
    else
        1, [2,3]
    end
    counts = dropdims(sum(bincounts(h), dims=dimremove), dims=dimremove)
    sumw2 = dropdims(sum(h.sumw2, dims=dimremove), dims=dimremove)
    edges = binedges(h)[dimskeep]
    return Hist2D(; binedges = edges, bincounts = counts, sumw2, nentries=nentries(h), overflow=h.overflow)
end
project(axis::Symbol=:x) = h::Union{Hist2D,Hist3D} -> project(h, axis)
