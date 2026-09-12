function auto_bins(ary, ::Val{2}; nbins=nothing)
    xs, ys = ary
    xnbins, ynbins = isnothing(nbins) ? _sturges.((xs, ys)) : _nbins_per_axis(nbins, Val(2))
    return (_auto_range(xs, xnbins), _auto_range(ys, ynbins))
end

"""
    sample(h::Hist2D; n::Int=1)

Sample a histogram's with weights equal to bin count, `n` times.
The sampled values are the bins' lower edges.
"""
function StatsBase.sample(h::Hist2D; n::Int=1)
    xedges, yedges = binedges(h)
    counts = bincounts(h)
    cis = CartesianIndices(counts)
    sampled = StatsBase.sample(cis, Weights(vec(counts)), n)
    xs = [xedges[I[1]] for I in sampled]
    ys = [yedges[I[2]] for I in sampled]
    return (xs, ys)
end

"""
    nbins(h::Hist2D)

Get a 2-tuple of the number of x and y bins of a histogram.
"""
function nbins(h::Hist2D)
    size(bincounts(h))
end

function integral(h::Hist2D; width=false)
    if width
        wx, wy = map(diff, h.binedges)
        return sum(bincounts(h) .* wx .* wy')
    else
        return sum(bincounts(h))
    end
end

"""
    push!(h::Hist2D, valx::Real, valy::Real, wgt::Real=1)
    atomic_push!(h::Hist2D, valx::Real, valy::Real, wgt::Real=1)

Adding one value at a time into histogram.
`sumw2` (sum of weights^2) accumulates `wgt^2` with a default weight of 1.
`atomic_push!` is a slower version of `push!` that is thread-safe.

Entries where any coordinate is outside of the bin edges are discarded (and not counted in
`nentries`), unless the histogram was created with `overflow=true`, in which case the
coordinates are clamped into the first/last bin along each axis. `NaN` is treated like `+Inf`.
"""
@inline function atomic_push!(h::Hist2D, valx::Real, valy::Real, wgt::Real=1)
    lock(h)
    push!(h, valx, valy, wgt)
    unlock(h)
    return nothing
end

@inline function Base.push!(h::Hist2D, valx::Real, valy::Real, wgt::Real=1)
    Lx, Ly = nbins(h)
    ix = _binindex(h.binedges[1], Lx, h.overflow, valx)
    iy = _binindex(h.binedges[2], Ly, h.overflow, valy)
    (ix == 0 || iy == 0) && return nothing
    h.nentries[] += 1
    @inbounds bincounts(h)[ix, iy] += wgt
    @inbounds sumw2(h)[ix, iy] += wgt^2
    return nothing
end

function Base.append!(h::Hist2D, xs::AbstractVector, ys::AbstractVector, wgts::AbstractVector)
    length(xs) == length(ys) == length(wgts) || throw(DimensionMismatch("append! to histogram expect same length values and weights"))
    lock(h)
    try
        for (x, y, w) in zip(xs, ys, wgts)
            push!(h, x, y, w)
        end
    finally
        unlock(h)
    end
    return h
end

function Base.append!(h::Hist2D, xs::AbstractVector, ys::AbstractVector)
    length(xs) == length(ys) || throw(DimensionMismatch("append! to histogram expect same length values along each axis"))
    lock(h)
    try
        for (x, y) in zip(xs, ys)
            push!(h, x, y)
        end
    finally
        unlock(h)
    end
    return h
end

for op in (:mean, :std, :median)
    @eval function Statistics.$op(h::Hist2D)
        px = project(h, :x)
        py = project(h, :y)
        return $op(px), $op(py)
    end
end

"""
    function lookup(h::Hist2D, x, y)

For given x-axis and y-axis value `x`, `y`, find the corresponding bin and return the bin content.
If a value is out of the histogram range, return `missing`.
"""
function lookup(h::Hist2D, x, y)
    rx, ry = binedges(h)
    !(first(rx) <= x < last(rx)) && return missing
    !(first(ry) <= y < last(ry)) && return missing
    return bincounts(h)[searchsortedlast(rx, x), searchsortedlast(ry, y)]
end


"""
    normalize(h::Hist2D; width=false)

Create a normalized histogram via division by `integral(h)`. When `width==true`, each bin is
additionally divided by its area such that `integral(normalize(h; width=true); width=true) == 1`.

!!! note
    Unlike for `Hist1D`, `width` defaults to `false` for backward compatibility.
"""
function normalize(h::Hist2D; width=false)
    hn = h*(1/integral(h; width=false))
    if width
        wx, wy = map(diff, h.binedges)
        hn.bincounts ./= wx .* wy'
        hn.sumw2 ./= (wx .* wy') .^ 2
    end
    return hn
end

"""
    rebin(h::Hist2D, nx::Int=1, ny::Int=nx)
    rebin(h::Hist2D, xedges::AbstractVector{<:Real}, yedges::AbstractVector{<:Real})
    rebin(nx::Int, ny::Int) = h::Hist2D -> rebin(h, nx, ny)

Merges `nx` (`ny`) consecutive bins into one along the x (y) axis by summing. Alternatively,
provide the new bin edges along each axis; they must be a subset of the existing edges (see
the `Hist1D` method of [`rebin`](@ref)).
"""
function rebin(h::Hist2D, nx::Int=1, ny::Int=nx)
    sx, sy = nbins(h)
    (sx % nx == 0 && sy % ny == 0) || _rebin_error(h, (nx, ny))
    bx, by = h.binedges
    blocks = (_rebin_blocks(sx, nx), _rebin_blocks(sy, ny))
    counts = _block_sum(bincounts(h), blocks)
    s2 = _block_sum(sumw2(h), blocks)
    ex = _subedges(bx, 1:nx:length(bx))
    ey = _subedges(by, 1:ny:length(by))
    return Hist2D(; binedges = (ex, ey), bincounts = counts, sumw2 = s2, nentries = nentries(h), overflow=h.overflow)
end
rebin(nx::Int, ny::Int) = h -> rebin(h, nx, ny)

function rebin(h::Hist2D, xedges::AbstractVector{<:Real}, yedges::AbstractVector{<:Real})
    bx, spans_x = _edge_blocks(h.binedges[1], xedges)
    by, spans_y = _edge_blocks(h.binedges[2], yedges)
    counts = _block_sum(bincounts(h), (bx, by))
    s2 = _block_sum(sumw2(h), (bx, by))
    return Hist2D(; binedges = (xedges, yedges), bincounts = counts, sumw2 = s2,
        nentries = nentries(h), overflow = h.overflow && spans_x && spans_y)
end


"""
    project(h::Hist2D, axis::Symbol=:x)
    project(axis::Symbol=:x) = h::Hist2D -> project(h, axis)

Computes the `:x` (`:y`) axis projection of the 2D histogram by
summing over the y (x) axis. Returns a `Hist1D`.

!!! note
    Beware that the `Hist3D` method of `project` has the opposite convention: there the given
    axis is the one being summed over (removed).
"""
function project(h::Hist2D, axis::Symbol=:x)
    axis ∈ (:x, :y) || throw(ArgumentError("axis must be ∈ `(:x, :y)`, got $axis"))
    dim = axis == :x ? 2 : 1
    counts = vec(sum(bincounts(h), dims=dim))
    s2 = vec(sum(sumw2(h), dims=dim))
    edges = axis == :x ? h.binedges[1] : h.binedges[2]
    return Hist1D(; binedges = edges, bincounts = counts, sumw2 = s2, nentries = nentries(h), overflow=h.overflow)
end

"""
    transpose(h::Hist2D)

Reverses the x and y axes.
"""
function transpose(h::Hist2D)
    edges = reverse(h.binedges)
    counts = permutedims(bincounts(h))
    s2 = permutedims(sumw2(h))
    return Hist2D(; binedges = edges, bincounts = counts, sumw2 = s2, nentries = nentries(h), overflow=h.overflow)
end

"""
    profile(h::Hist2D, axis::Symbol=:x)
    profile(axis::Symbol=:x) = h::Hist2D -> profile(h, axis)

Returns the `axis`-profile of the 2D histogram by
calculating the weighted mean over the other axis.
`profile(h, :x)` will return a `Hist1D` with the y-axis edges of `h`.
"""
function profile(h::Hist2D, axis::Symbol=:x)
    axis ∈ (:x, :y) || throw(ArgumentError("axis must be ∈ `(:x, :y)`, got $axis"))
    if axis == :y
        h = transpose(h)
    end

    edges = h.binedges[1]
    centers = bincenters(h)[2]
    counts = bincounts(h)
    _sumw2 = sumw2(h)

    num = counts*centers
    den = sum(counts, dims=2)
    numerr2 = _sumw2 * centers.^2
    denerr2 = sum(_sumw2, dims=2)
    val = vec(num ./ den)
    sw2 = vec(@. numerr2/den^2 - denerr2*(num/den^2)^2)

    # ROOT sets the NaN entries and their error to 0
    val[isnan.(val)] .= zero(eltype(val))
    sw2[isnan.(sw2)] .= zero(eltype(sw2))

    return Hist1D(; binedges = edges, bincounts = val, sumw2 = sw2, nentries = nentries(h), overflow=h.overflow)
end
profile(axis::Symbol=:x) = h::Hist2D -> profile(h, axis)

"""
    restrict(h::Hist2D, xlow=-Inf, xhigh=Inf, ylow=-Inf, yhigh=Inf)
    restrict(xlow=-Inf, xhigh=Inf, ylow=-Inf, yhigh=Inf) = h::Hist2D -> restrict(h, xlow, xhigh, ylow, yhigh)

Returns a new histogram with a restricted x-axis.
`restrict(h, 0, 3)` (or `h |> restrict(0, 3)`)
will return a slice of `h` where the bin centers are in `[0, 3]` (inclusive).
"""
function restrict(h::Hist2D, xlow=-Inf, xhigh=Inf, ylow=-Inf, yhigh=Inf)
    bx, by = h.binedges
    xsel = _restrict_bins(bx, xlow, xhigh)
    ysel = _restrict_bins(by, ylow, yhigh)
    xedges = _subedges(bx, first(xsel):last(xsel)+1)
    yedges = _subedges(by, first(ysel):last(ysel)+1)
    c = bincounts(h)[xsel, ysel]
    s2 = sumw2(h)[xsel, ysel]
    Hist2D(; binedges = (xedges, yedges), bincounts = c, sumw2 = s2, nentries = nentries(h), overflow=h.overflow)
end
restrict(xlow, xhigh, ylow, yhigh) = h::Hist2D -> restrict(h, xlow, xhigh, ylow, yhigh)
