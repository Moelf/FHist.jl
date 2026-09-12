function auto_bins(ary, ::Val{3}; nbins=nothing)
    xs, ys, zs = ary
    xnbins, ynbins, znbins = isnothing(nbins) ? _sturges.((xs, ys, zs)) : _nbins_per_axis(nbins, Val(3))
    return (_auto_range(xs, xnbins), _auto_range(ys, ynbins), _auto_range(zs, znbins))
end

"""
    sample(h::Hist3D; n::Int=1)

Sample a histogram's with weights equal to bin count, `n` times.
The sampled values are the bins' lower edges.
"""
function StatsBase.sample(h::Hist3D; n::Int=1)
    xedges, yedges, zedges = binedges(h)
    counts = bincounts(h)
    cis = CartesianIndices(counts)
    sampled = StatsBase.sample(cis, Weights(vec(counts)), n)
    xs = [xedges[I[1]] for I in sampled]
    ys = [yedges[I[2]] for I in sampled]
    zs = [zedges[I[3]] for I in sampled]
    return (xs, ys, zs)
end

"""
    nbins(h::Hist3D)

Get a 3-tuple of the number of x, y and z bins of a histogram.
"""
function nbins(h::Hist3D)
    size(bincounts(h))
end

function integral(h::Hist3D; width=false)
    if width
        wx, wy, wz = map(diff, h.binedges)
        return sum(bincounts(h) .* wx .* wy' .* reshape(wz, 1, 1, :))
    else
        return sum(bincounts(h))
    end
end

"""
    push!(h::Hist3D, valx::Real, valy::Real, valz::Real, wgt::Real=1)
    atomic_push!(h::Hist3D, valx::Real, valy::Real, valz::Real, wgt::Real=1)

Adding one value at a time into histogram.
`sumw2` (sum of weights^2) accumulates `wgt^2` with a default weight of 1.
`atomic_push!` is a slower version of `push!` that is thread-safe.

Entries where any coordinate is outside of the bin edges are discarded (and not counted in
`nentries`), unless the histogram was created with `overflow=true`, in which case the
coordinates are clamped into the first/last bin along each axis. `NaN` is treated like `+Inf`.
"""
@inline function atomic_push!(h::Hist3D, valx::Real, valy::Real, valz::Real, wgt::Real=1)
    lock(h)
    push!(h, valx, valy, valz, wgt)
    unlock(h)
    return nothing
end

@inline function Base.push!(h::Hist3D, valx::Real, valy::Real, valz::Real, wgt::Real=1)
    Lx, Ly, Lz = nbins(h)
    ix = _binindex(h.binedges[1], Lx, h.overflow, valx)
    iy = _binindex(h.binedges[2], Ly, h.overflow, valy)
    iz = _binindex(h.binedges[3], Lz, h.overflow, valz)
    (ix == 0 || iy == 0 || iz == 0) && return nothing
    h.nentries[] += 1
    @inbounds bincounts(h)[ix, iy, iz] += wgt
    @inbounds sumw2(h)[ix, iy, iz] += wgt^2
    return nothing
end

function Base.append!(h::Hist3D, xs::AbstractVector, ys::AbstractVector, zs::AbstractVector, wgts::AbstractVector)
    length(xs) == length(ys) == length(zs) == length(wgts) || throw(DimensionMismatch("append! to histogram expect same length values and weights"))
    lock(h)
    try
        for (x, y, z, w) in zip(xs, ys, zs, wgts)
            push!(h, x, y, z, w)
        end
    finally
        unlock(h)
    end
    return h
end

function Base.append!(h::Hist3D, xs::AbstractVector, ys::AbstractVector, zs::AbstractVector)
    length(xs) == length(ys) == length(zs) || throw(DimensionMismatch("append! to histogram expect same length values along each axis"))
    lock(h)
    try
        for (x, y, z) in zip(xs, ys, zs)
            push!(h, x, y, z)
        end
    finally
        unlock(h)
    end
    return h
end

# 1D projection onto a single axis (summing over the other two)
function _project1d(h::Hist3D, dim::Int)
    others = filter(!=(dim), (1, 2, 3))
    counts = vec(sum(bincounts(h), dims=others))
    s2 = vec(sum(sumw2(h), dims=others))
    return Hist1D(; binedges = h.binedges[dim], bincounts = counts, sumw2 = s2, nentries = nentries(h), overflow = h.overflow)
end

for op in (:mean, :std, :median)
    @eval function Statistics.$op(h::Hist3D)
        return $op(_project1d(h, 1)), $op(_project1d(h, 2)), $op(_project1d(h, 3))
    end
end

"""
    lookup(h::Hist3D, x, y, z)

For given x/y/z-axis value `x`, `y`, `z`, find the corresponding bin and return the bin content.
If a value is out of the histogram range, return `missing`.
"""
function lookup(h::Hist3D, x, y, z)
    rx, ry, rz = binedges(h)
    !(first(rx) <= x < last(rx)) && return missing
    !(first(ry) <= y < last(ry)) && return missing
    !(first(rz) <= z < last(rz)) && return missing
    return bincounts(h)[searchsortedlast(rx, x), searchsortedlast(ry, y), searchsortedlast(rz, z)]
end


"""
    normalize(h::Hist3D; width=false)

Create a normalized histogram via division by `integral(h)`. When `width==true`, each bin is
additionally divided by its volume such that `integral(normalize(h; width=true); width=true) == 1`.

!!! note
    Unlike for `Hist1D`, `width` defaults to `false` for backward compatibility.
"""
function normalize(h::Hist3D; width=false)
    hn = h*(1/integral(h; width=false))
    if width
        wx, wy, wz = map(diff, h.binedges)
        vol = wx .* wy' .* reshape(wz, 1, 1, :)
        hn.bincounts ./= vol
        hn.sumw2 ./= vol .^ 2
    end
    return hn
end

"""
    rebin(h::Hist3D, nx::Int=1, ny::Int=nx, nz::Int=nx)
    rebin(h::Hist3D, xedges::AbstractVector{<:Real}, yedges::AbstractVector{<:Real}, zedges::AbstractVector{<:Real})
    rebin(nx::Int, ny::Int, nz::Int) = h::Hist3D -> rebin(h, nx, ny, nz)

Merges `nx` (`ny`, `nz`) consecutive bins into one along the x (y, z) axis by summing.
Alternatively, provide the new bin edges along each axis; they must be a subset of the existing
edges (see the `Hist1D` method of [`rebin`](@ref)).
"""
function rebin(h::Hist3D, nx::Int=1, ny::Int=nx, nz::Int=nx)
    sx, sy, sz = nbins(h)
    (sx % nx == 0 && sy % ny == 0 && sz % nz == 0) || _rebin_error(h, (nx, ny, nz))
    bx, by, bz = h.binedges
    blocks = (_rebin_blocks(sx, nx), _rebin_blocks(sy, ny), _rebin_blocks(sz, nz))
    counts = _block_sum(bincounts(h), blocks)
    s2 = _block_sum(sumw2(h), blocks)
    ex = _subedges(bx, 1:nx:length(bx))
    ey = _subedges(by, 1:ny:length(by))
    ez = _subedges(bz, 1:nz:length(bz))
    return Hist3D(; binedges = (ex, ey, ez), bincounts = counts, sumw2 = s2, nentries = nentries(h), overflow=h.overflow)
end
rebin(nx::Int, ny::Int, nz::Int) = h::Hist3D -> rebin(h, nx, ny, nz)

function rebin(h::Hist3D, xedges::AbstractVector{<:Real}, yedges::AbstractVector{<:Real}, zedges::AbstractVector{<:Real})
    bx, spans_x = _edge_blocks(h.binedges[1], xedges)
    by, spans_y = _edge_blocks(h.binedges[2], yedges)
    bz, spans_z = _edge_blocks(h.binedges[3], zedges)
    counts = _block_sum(bincounts(h), (bx, by, bz))
    s2 = _block_sum(sumw2(h), (bx, by, bz))
    return Hist3D(; binedges = (xedges, yedges, zedges), bincounts = counts, sumw2 = s2,
        nentries = nentries(h), overflow = h.overflow && spans_x && spans_y && spans_z)
end

"""
    project(h::Hist3D, axis::Symbol=:x)
    project(axis::Symbol=:x) = h::Hist3D -> project(h, axis)

Computes the `:x`/`:y`/`:z` axis projection of the 3D histogram by
summing over the specified axis. Returns a `Hist2D`.

!!! note
    Beware that this is the opposite convention of the `Hist2D` method of `project`, where the
    given axis is the one that is *kept*. E.g. `project(h3, :z)` gives the `(x, y)` `Hist2D`, and
    `project(project(h3, :z), :x)` gives the x-axis `Hist1D`.
"""
function project(h::Hist3D, axis::Symbol=:x)
    axis ∈ (:x, :y, :z) || throw(ArgumentError("axis must be ∈ `(:x, :y, :z)`, got $axis"))
    dimremove, dimskeep = if axis == :z
        3, (1,2)
    elseif axis == :y
        2, (1,3)
    else
        1, (2,3)
    end
    counts = dropdims(sum(bincounts(h), dims=dimremove), dims=dimremove)
    s2 = dropdims(sum(sumw2(h), dims=dimremove), dims=dimremove)
    edges = (h.binedges[dimskeep[1]], h.binedges[dimskeep[2]])
    return Hist2D(; binedges = edges, bincounts = counts, sumw2 = s2, nentries=nentries(h), overflow=h.overflow)
end
project(axis::Symbol=:x) = h::Union{Hist2D,Hist3D} -> project(h, axis)

"""
    restrict(h::Hist3D, xlow=-Inf, xhigh=Inf, ylow=-Inf, yhigh=Inf, zlow=-Inf, zhigh=Inf)
    restrict(xlow, xhigh, ylow, yhigh, zlow, zhigh) = h::Hist3D -> restrict(h, xlow, xhigh, ylow, yhigh, zlow, zhigh)

Returns a new histogram with restricted axes: the slice of `h` where the bin centers are within
the given (inclusive) intervals along each axis.
"""
function restrict(h::Hist3D, xlow=-Inf, xhigh=Inf, ylow=-Inf, yhigh=Inf, zlow=-Inf, zhigh=Inf)
    bx, by, bz = h.binedges
    xsel = _restrict_bins(bx, xlow, xhigh)
    ysel = _restrict_bins(by, ylow, yhigh)
    zsel = _restrict_bins(bz, zlow, zhigh)
    xedges = _subedges(bx, first(xsel):last(xsel)+1)
    yedges = _subedges(by, first(ysel):last(ysel)+1)
    zedges = _subedges(bz, first(zsel):last(zsel)+1)
    c = bincounts(h)[xsel, ysel, zsel]
    s2 = sumw2(h)[xsel, ysel, zsel]
    Hist3D(; binedges = (xedges, yedges, zedges), bincounts = c, sumw2 = s2, nentries = nentries(h), overflow=h.overflow)
end
restrict(xlow, xhigh, ylow, yhigh, zlow, zhigh) = h::Hist3D -> restrict(h, xlow, xhigh, ylow, yhigh, zlow, zhigh)
