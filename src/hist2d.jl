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
To be used with [`push!`](@ref). Default overflow behavior (`false`)
will exclude values that are outside of `binedges`.
"""
function Hist2D(elT::Type{T}=Float64; bins, overflow=_default_overflow) where {T}
    counts = zeros(elT, (length(bins[1]) - 1, length(bins[2])-1))
    return Hist2D(Histogram(bins, counts); overflow=overflow)
end

"""
    Hist2D(tuple, edges::Tuple{AbstractRange,AbstractRange}; overflow)
    Hist2D(tuple, edges::Tuple{AbstractVector,AbstractVector}; overflow)

Create a `Hist2D` with given bin `edges` and values from
a 2-tuple of arrays of x, y values. Weight for each value is assumed to be 1.
"""
function Hist2D(A::Tuple{AbstractVector,AbstractVector}, r::Tuple{AbstractRange,AbstractRange}; overflow=_default_overflow)
    h = Hist2D(Int; bins=r, overflow=overflow)
    unsafe_push!.(h, A[1], A[2])
    return h
end
function Hist2D(A::Tuple{AbstractVector,AbstractVector}, edges::Tuple{AbstractVector,AbstractVector}; overflow=_default_overflow)
    if all(_is_uniform_bins.(edges))
        r = (range(first(edges[1]), last(edges[1]), length=length(edges[1])),
             range(first(edges[2]), last(edges[2]), length=length(edges[2])))
        return Hist2D(A, r; overflow=overflow)
    else
        h = Hist2D(Int; bins=edges, overflow=overflow)
        unsafe_push!.(h, A[1], A[2])
        return h
    end
end

"""
    Hist2D(tuple, wgts::AbstractWeights, edges::Tuple{AbstractRange,AbstractRange}; overflow)
    Hist2D(tuple, wgts::AbstractWeights, edges::Tuple{AbstractVector,AbstractVector}; overflow)

Create a `Hist2D` with given bin `edges` and values from
a 2-tuple of arrays of x, y values.
`wgts` should have the same `size` as elements of `tuple`.
"""
function Hist2D(A::Tuple{AbstractVector,AbstractVector}, wgts::AbstractWeights, r::Tuple{AbstractRange,AbstractRange}; overflow=_default_overflow)
    @boundscheck @assert size(A[1]) == size(A[2]) == size(wgts)
    h = Hist2D(eltype(wgts); bins=r, overflow=overflow)
    unsafe_push!.(h, A[1], A[2], wgts)
    return h
end
function Hist2D(A::Tuple{AbstractVector,AbstractVector}, wgts::AbstractWeights, edges::Tuple{AbstractVector,AbstractVector}; overflow=_default_overflow)
    if all(_is_uniform_bins.(edges))
        r = (range(first(edges[1]), last(edges[1]), length=length(edges[1])),
             range(first(edges[2]), last(edges[2]), length=length(edges[2])))
        return Hist2D(A, wgts, r; overflow=overflow)
    else
        h = Hist2D(Int; bins=edges, overflow=overflow)
        unsafe_push!.(h, A[1], A[2], wgts)
        return h
    end
end

"""
    Hist2D(A::AbstractVector{T}; nbins::Tuple{Integer,Integer}, overflow) where T
    Hist2D(A::AbstractVector{T}, wgts::AbstractWeights; nbins::Tuple{Integer,Integer}, overflow) where T

Automatically determine number of bins based on `Sturges` algo.
"""
function Hist2D(A::Tuple{AbstractVector{T},AbstractVector{T}};
        nbins::Tuple{Integer,Integer}=_sturges.(A),
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

function Hist2D(A::Tuple{AbstractVector{T},AbstractVector{T}}, wgts::AbstractWeights;
        nbins::Tuple{Integer,Integer}=_sturges.(A),
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
    @assert sx % nx == sy % ny == 0
    p1d = (x,n)->Iterators.partition(x, n)
    p2d = x->(x[i:i+(nx-1),j:j+(ny-1)] for i=1:nx:sx, j=1:ny:sy)
    counts = sum.(p2d(bincounts(h)))
    sumw2 = sum.(p2d(h.sumw2))
    ex = first.(p1d(binedges(h)[1], nx))
    ey = first.(p1d(binedges(h)[2], ny))
    _is_uniform_bins(ex) && (ex = range(first(ex), last(ex), length=length(ex)))
    _is_uniform_bins(ey) && (ey = range(first(ey), last(ey), length=length(ey)))
    return Hist2D(Histogram((ex,ey), counts), sumw2; overflow=h.overflow)
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
    return Hist1D(Histogram(edges, counts), sumw2; overflow=h.overflow)
end
project(axis::Symbol=:x) = h::Hist2D -> project(h, axis)

"""
    transpose(h::Hist2D)

Reverses the x and y axes.
"""
function transpose(h::Hist2D)
    edges = reverse(binedges(h))
    counts = collect(bincounts(h)')
    sumw2 = collect(h.sumw2')
    return Hist2D(Histogram(edges, counts), sumw2; overflow=h.overflow)
end

"""
    profile(h::Hist2D, axis::Symbol=:x)
    profile(axis::Symbol=:x) = h::Hist2D -> profile(h, axis)

Returns the `axis`-profile of the 2D histogram by
calculating the weighted mean over the other axis.
`profile(h, :x)` will return a `Hist1D` with the y-axis edges of `h`.
"""
function profile(h::Hist2D, axis::Symbol=:x)
    @assert axis ∈ (:x, :y)
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

    return Hist1D(Histogram(edges, val), sw2; overflow=h.overflow)
end
profile(axis::Symbol=:x) = h::Hist2D -> profile(h, axis)

function _svg(h::Hist2D)
    paddingx1, paddingy1 = 0.15, 0.15 # left, bottom
    paddingx2, paddingy2 = 0.05, 0.05 # right, top
    framewidth, frameheight = 250, 200

    (xlow, xhigh), (ylow, yhigh) = extrema.(binedges(h))
    function transform(x::Real, y::Real)
        xfrac = (x - xlow)/(xhigh - xlow)
        xdraw = xfrac * framewidth*(1-paddingx1-paddingx2) + framewidth*paddingx1
        yfrac = (y - ylow)/(yhigh - ylow)
        ydraw = frameheight - (yfrac * frameheight*(1-paddingy1-paddingy2) + frameheight*paddingy1)
        return xdraw, ydraw
    end

    function colorscale(x::Real)
        x = isnan(x) ? 0. : x
        # x is fraction between [0,1]
        # 6th deg polynomial fit to viridis color palette
        # See https://gist.github.com/aminnj/dba84777718613d3a37291a43659feff
        # to generate other color palettes. Black and white is trivial:
        #   return 255 .* (1-x,1-x,1-x)
        r = @evalpoly(x, 0.2731, 0.1270, -0.3617, -4.7456, 6.7092, 4.2491, -5.2678)
        g = @evalpoly(x, 0.0039, 1.3810, 0.3969, -6.4246, 15.3241, -14.7345, 4.9587)
        b = @evalpoly(x, 0.3305, 1.3726, 0.3948, -20.7431, 59.7287, -68.3565, 27.4051)
        return round.(Int, 255 .* (r,g,b))
    end

    counts = bincounts(h)
    mincount, maxcount = extrema(counts)
    counts = clamp.(counts, 0, maxcount)
    ex, ey = binedges(h)
    sx, sy = size(counts)
    rectlines = String[]
    textstyle = """ dominant-baseline="middle" text-anchor="middle" font-size="85%" font-family="sans-serif" """
    for i in 1:sx, j in 1:sy
        tx1, ty1 = ceil.(Int, transform(ex[i], ey[j+1]))
        tx2, ty2 = ceil.(Int, transform(ex[i+1], ey[j]))
        tw, th = tx2-tx1, ty2-ty1
        c = counts[i,j]
        (c == 0) && continue
        r,g,b = colorscale((c-mincount)/(maxcount-mincount))
        rcolor = "rgb($(r),$(g),$(b))"
        line = """<rect x="$(tx1)" y="$(ty1)" width="$(tw)" height="$(th)" fill="$(rcolor)" stroke="none" />"""
        if (sx <= 15) && (sy <= 15)
            tcolor = ((0.299*r + 0.587*g + 0.114*b) < 0.60*255) ? "#fff" : "#000"
            line = """<g>$(line)<text class="svgplotlabels" x="$(tx1+tw/2)" y="$(ty1+th/2)" fill="$(tcolor)" $(textstyle)>$(c)</text></g>"""
        end
        push!(rectlines, line * "\n")
    end

    return """
    <svg width="$(framewidth)" height="$(frameheight)" version="1.1" xmlns="http://www.w3.org/2000/svg" class="svgplot">
        <style>
          .svgplotlabels { display: none; }
          .svgplot g:hover text { display: block; cursor: default; }
        </style>
        $(join(rectlines))
        <rect x="$(round(Int,framewidth*paddingx1))" y="$(round(Int,frameheight*paddingy2))" width="$(round(Int,framewidth*(1-paddingx1-paddingx2)))" height="$(frameheight*(1-paddingy1-paddingy2))" fill="none" stroke="#000" stroke-width="1" />
        <text x="$(framewidth*paddingx1)" y="$(frameheight*(1-0.5*paddingy1))" fill="black" $(textstyle)>$(minimum(ex))</text>
        <text x="$(framewidth*(1-paddingx2))" y="$(frameheight*(1-0.5*paddingy1))" fill="black" $(textstyle)>$(maximum(ex))</text>
        <text x="$(framewidth*0.5*paddingx1)" y="$(frameheight*(1-paddingy1))" fill="black" $(textstyle)>$(minimum(ey))</text>
        <text x="$(framewidth*0.5*paddingx1)" y="$(frameheight*paddingy2)" fill="black" $(textstyle)>$(maximum(ey))</text>
    </svg>
    """
end

function Base.show(io::IO, h::Hist2D)
    ex, ey = binedges(h)
    nx, ny = nbins(h)
    xscale = nx > 1 ? (maximum(ex)-minimum(ex))/(nx-1) : 0.0
    yscale = ny > 1 ? (maximum(ey)-minimum(ey))/(ny-1) : 0.0
    show(io, UnicodePlots.heatmap(bincounts(h)'; xscale=xscale, xoffset=minimum(ex), yscale=yscale, yoffset=minimum(ey)))
    println(io)
    println(io, "edges: ", binedges(h))
    println(io, "bin counts: ", bincounts(h))
    print(io, "total count: ", integral(h))
end

function Base.show(io::IO, m::MIME"text/html", h::Hist2D)
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
