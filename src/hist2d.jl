Base.lock(h::Hist2D) = lock(h.hlock)
Base.unlock(h::Hist2D) = unlock(h.hlock)

"""
    bincounts(h::Hist2D)

Get the bin counts of a histogram.
"""
@inline bincounts(h::Hist2D) = h.hist.weights

"""
    binedges(h::Hist2D)

Get a 2-tuple of the bin edges of a histogram.
"""
@inline binedges(h::Hist2D) = h.hist.edges

"""
    bincenters(h::Hist2D)

Get a 2-tuple of the bin centers of a histogram.
"""
function bincenters(h::Hist2D)
    StatsBase.midpoints.(binedges(h))
end

"""
    binerrors(f::T, h::Hist2D) where T<:Function = f.(h.sumw2)
    binerrors(h::Hist2D) = binerrors(sqrt, h)

Calculate the bin errors from `sumw2` with a Gaussian default.
"""
binerrors(f::T, h::Hist2D) where T<:Function = f.(h.sumw2)
binerrors(h::Hist2D) = binerrors(sqrt, h)


"""
    nbins(h::Hist2D)

Get a 2-tuple of the number of x and y bins of a histogram.
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
    push!(h::Hist2D, valx::Real, valy::Real, wgt::Real=1)
    atomic_push!(h::Hist2D, valx::Real, valy::Real, wgt::Real=1)

Adding one value at a time into histogram.
`sumw2` (sum of weights^2) accumulates `wgt^2` with a default weight of 1.
`atomic_push!` is a slower version of `push!` that is thread-safe.

"""
@inline function atomic_push!(h::Hist2D{T,E}, valx::Real, valy::Real, wgt::Real=1) where {T,E}
    lock(h)
    push!(h, valx, valy, wgt)
    unlock(h)
    return nothing
end

@inline function Base.push!(h::Hist2D{T,E}, valx::Real, valy::Real, wgt::Real=1) where {T,E}
    rx, ry = binedges(h)
    Lx, Ly = nbins(h)
    binidxx = _edge_binindex(rx, valx)
    binidxy = _edge_binindex(ry, valy)
    h.nentries[] += 1
    if h.overflow
        binidxx = clamp(binidxx, 1, Lx)
        binidxy = clamp(binidxy, 1, Ly)
        @inbounds h.hist.weights[binidxx,binidxy] += wgt
        @inbounds h.sumw2[binidxx,binidxy] += wgt^2
    else
        if (unsigned(binidxx - 1) < Lx) && (unsigned(binidxy - 1) < Ly)
            @inbounds h.hist.weights[binidxx,binidxy] += wgt
            @inbounds h.sumw2[binidxx,binidxy] += wgt^2
        end
    end
    return nothing
end

Base.broadcastable(h::Hist2D) = Ref(h)

"""
    Hist2D(elT::Type{T}=Float64; binedges, overflow) where {T}

Initialize an empty histogram with bin content typed as `T` and bin edges.
To be used with `push!` or [`atomic_push!`](@ref). Default overflow behavior (`false`)
will exclude values that are outside of `binedges`.
"""
function Hist2D(elT::Type{T}=Float64; bins, overflow=_default_overflow) where {T}
    counts = zeros(elT, length.(bins) .- 1)
    return Hist2D(Histogram(bins, counts); overflow=overflow)
end

"""
    Hist2D(Tuple{<:iterable, <:iterable}, edges::NTuple{2,AbstractRange}; overflow)
    Hist2D(Tuple{<:iterable, <:iterable}, edges::NTuple{2,AbstractVector}; overflow)

Create a `Hist2D` with given bin `edges` and values from
a 2-tuple of arrays of x, y values. Weight for each value is assumed to be 1.
"""
function Hist2D(A, r::NTuple{2,AbstractRange}; overflow=_default_overflow)
    A1, A2 = A
    h = Hist2D(Int; bins=r, overflow=overflow)
    push!.(h, A1, A2)
    return h
end
function Hist2D(A, edges::NTuple{2,AbstractVector}; overflow=_default_overflow)
    A1, A2 = A
    if all(_is_uniform_bins.(edges))
        r = (range(first(edges[1]), last(edges[1]), length=length(edges[1])),
             range(first(edges[2]), last(edges[2]), length=length(edges[2])))
        return Hist2D(A, r; overflow=overflow)
    else
        h = Hist2D(Int; bins=edges, overflow=overflow)
        push!.(h, A1, A2)
        return h
    end
end

"""
    Hist2D(tuple, wgts::AbstractWeights, edges::NTuple{2,AbstractRange}; overflow)
    Hist2D(tuple, wgts::AbstractWeights, edges::NTuple{2,AbstractVector}; overflow)

Create a `Hist2D` with given bin `edges` and values from
a 2-tuple of arrays of x, y values.
`wgts` should have the same `size` as elements of `tuple`.
"""
function Hist2D(A::NTuple{2,AbstractVector}, wgts::AbstractWeights, r::NTuple{2,AbstractRange}; overflow=_default_overflow)
    @boundscheck @assert size(A[1]) == size(A[2]) == size(wgts)
    h = Hist2D(eltype(wgts); bins=r, overflow=overflow)
    push!.(h, A[1], A[2], wgts)
    return h
end
function Hist2D(A::NTuple{2,AbstractVector}, wgts::AbstractWeights, edges::NTuple{2,AbstractVector}; overflow=_default_overflow)
    if all(_is_uniform_bins.(edges))
        r = (range(first(edges[1]), last(edges[1]), length=length(edges[1])),
             range(first(edges[2]), last(edges[2]), length=length(edges[2])))
        return Hist2D(A, wgts, r; overflow=overflow)
    else
        h = Hist2D(Int; bins=edges, overflow=overflow)
        push!.(h, A[1], A[2], wgts)
        return h
    end
end

"""
    Hist2D(XsYs::Tuple; nbins::NTuple{2,Integer}, overflow) where T
    Hist2D(XsYs::Tuple, wgts::AbstractWeights; nbins::NTuple{2,Integer}, overflow) where T

Automatically determine number of bins based on `Sturges` algo.
"""
function Hist2D(A::NTuple{2,AbstractVector{T}};
        nbins::NTuple{2,Integer}=_sturges.(A),
        overflow=_default_overflow,
    ) where {T}
    F = float(T)
    nbinsx, nbinsy = nbins
    lox, hix = minimum(A[1]), maximum(A[1])
    loy, hiy = minimum(A[2]), maximum(A[2])
    rx = StatsBase.histrange(F(lox), F(hix), nbinsx)
    ry = StatsBase.histrange(F(loy), F(hiy), nbinsy)
    r = (rx, ry)
    return Hist2D(A, r; overflow=overflow)
end

function Hist2D(A::NTuple{2,AbstractVector{T}}, wgts::AbstractWeights;
        nbins::NTuple{2,Integer}=_sturges.(A),
        overflow=_default_overflow,
    ) where {T}
    F = float(T)
    nbinsx, nbinsy = nbins
    lox, hix = minimum(A[1]), maximum(A[1])
    loy, hiy = minimum(A[2]), maximum(A[2])
    rx = StatsBase.histrange(F(lox), F(hix), nbinsx)
    ry = StatsBase.histrange(F(loy), F(hiy), nbinsy)
    r = (rx, ry)
    return Hist2D(A, wgts, r; overflow=overflow)
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

"""
    rebin(h::Hist2D, nx::Int=1, ny::Int=nx)
    rebin(nx::Int, ny::Int) = h::Hist2D -> rebin(h, nx, ny)

Merges `nx` (`ny`) consecutive bins into one along the x (y) axis by summing.
"""
function rebin(h::Hist2D, nx::Int=1, ny::Int=nx)
    sx, sy = nbins(h)
    if !(sx % nx == sy % ny == 0)
        rebin_values_x, rebin_values_y = map(x->join(sort(collect(x)), ", ", " or "), valid_rebin_values(h))
        error("Invalid rebin values (nx: $nx, ny: $ny) for a 2D histogram with $(nbins(h)) bins. They have to be a multple of nx: $(rebin_values_x) and ny: $(rebin_values_y)")
    end
    p1d = (x,n)->Iterators.partition(x, n)
    p2d = x->(x[i:i+(nx-1),j:j+(ny-1)] for i=1:nx:sx, j=1:ny:sy)
    counts = sum.(p2d(bincounts(h)))
    sumw2 = sum.(p2d(h.sumw2))
    ex = first.(p1d(binedges(h)[1], nx))
    ey = first.(p1d(binedges(h)[2], ny))
    _is_uniform_bins(ex) && (ex = range(first(ex), last(ex), length=length(ex)))
    _is_uniform_bins(ey) && (ey = range(first(ey), last(ey), length=length(ey)))
    return Hist2D(Histogram((ex,ey), counts), sumw2, nentries(h); overflow=h.overflow)
end
rebin(nx::Int, ny::Int) = h::Hist2D -> rebin(h, nx, ny)


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
    return Hist1D(Histogram(edges, counts), sumw2, nentries(h); overflow=h.overflow)
end

"""
    transpose(h::Hist2D)

Reverses the x and y axes.
"""
function transpose(h::Hist2D)
    edges = reverse(binedges(h))
    counts = collect(bincounts(h)')
    sumw2 = collect(h.sumw2')
    return Hist2D(Histogram(edges, counts), sumw2, nentries(h); overflow=h.overflow)
end

"""
    profile(h::Hist2D, axis::Symbol=:x)
    profile(axis::Symbol=:x) = h::Hist2D -> profile(h, axis)

Returns the `axis`-profile of the 2D histogram by
calculating the weighted mean over the other axis.
`profile(h, :x)` will return a `Hist1D` with the y-axis edges of `h`.
"""
function profile(h::Hist2D, axis::Symbol=:x)
    axis ∈ (:x, :y) || throw("axis must be ∈ `(:x, :y)`, got $axis")
    if axis == :y
        h = transpose(h)
    end

    edges = binedges(h)[1]
    centers = bincenters(h)[2]
    counts = bincounts(h)
    sumw2 = h.sumw2

    num = counts*centers
    den = sum(counts, dims=2)
    numerr2 = sumw2 * centers.^2
    denerr2 = sum(sumw2, dims=2)
    val = vec(num ./ den)
    sw2 = vec(@. numerr2/den^2 - denerr2*(num/den^2)^2)

    # ROOT sets the NaN entries and their error to 0
    val[isnan.(val)] .= zero(eltype(val))
    sw2[isnan.(sw2)] .= zero(eltype(sw2))

    return Hist1D(Histogram(edges, val), sw2, nentries(h); overflow=h.overflow)
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
    xsel = xlow .<= bincenters(h)[1] .<= xhigh
    ysel = ylow .<= bincenters(h)[2] .<= yhigh
    @assert count(xsel) > 0 "No bin centers contained in [$(xlow), $(xhigh)]"
    @assert count(ysel) > 0 "No bin centers contained in [$(ylow), $(yhigh)]"
    xedgesel = push!(copy(xsel), false)
    yedgesel = push!(copy(ysel), false)

    xlastidx = findlast(xedgesel)
    if xlastidx !== nothing
        xedgesel[xlastidx+1] = 1
    end

    ylastidx = findlast(yedgesel)
    if ylastidx !== nothing
        yedgesel[ylastidx+1] = 1
    end

    xedges = binedges(h)[1][xedgesel]
    if _is_uniform_bins(xedges)
        xedges = range(first(xedges), last(xedges), length=length(xedges))
    end

    yedges = binedges(h)[2][yedgesel]
    if _is_uniform_bins(yedges)
        yedges = range(first(yedges), last(yedges), length=length(yedges))
    end

    c = bincounts(h)[xsel,ysel]
    sumw2 = h.sumw2[xsel,ysel]

    Hist2D(Histogram((xedges,yedges), c), sumw2, nentries(h); overflow=h.overflow)
end
restrict(xlow, xhigh, ylow, yhigh) = h::Hist2D -> restrict(h, xlow, xhigh, ylow, yhigh)
