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
    if h.overflow
        binidxx = clamp(binidxx, 1, Lx)
        binidxy = clamp(binidxy, 1, Ly)
        binidxz = clamp(binidxz, 1, Lz)
        @inbounds h.hist.weights[binidxx,binidxy,binidxz] += wgt
        @inbounds h.sumw2[binidxx,binidxy,binidxz] += wgt^2
    else
        if (unsigned(binidxx - 1) < Lx) && (unsigned(binidxy - 1) < Ly) && (unsigned(binidxz - 1) < Lz)
            @inbounds h.hist.weights[binidxx,binidxy,binidxz] += wgt
            @inbounds h.sumw2[binidxx,binidxy,binidxz] += wgt^2
        end
    end
    return nothing
end

Base.broadcastable(h::Hist3D) = Ref(h)

"""
    Hist3D(elT::Type{T}=Float64; binedges, overflow) where {T}

Initialize an empty histogram with bin content typed as `T` and bin edges.
To be used with [`push!`](@ref). Default overflow behavior (`false`)
will exclude values that are outside of `binedges`.
"""
function Hist3D(elT::Type{T}=Float64; bins, overflow=_default_overflow) where {T}
    counts = zeros(elT, length.(bins) .- 1)
    return Hist3D(Histogram(bins, counts); overflow=overflow)
end

"""
    Hist3D(tuple, edges::Tuple{AbstractRange,AbstractRange}; overflow)
    Hist3D(tuple, edges::Tuple{AbstractVector,AbstractVector}; overflow)

Create a `Hist3D` with given bin `edges` and values from
a 2-tuple of arrays of x, y values. Weight for each value is assumed to be 1.
"""
function Hist3D(A::NTuple{3,AbstractVector}, r::NTuple{3,AbstractRange}; overflow=_default_overflow)
    h = Hist3D(Int; bins=r, overflow=overflow)
    push!.(h, A[1], A[2], A[3])
    return h
end
function Hist3D(A::NTuple{3,AbstractVector}, edges::NTuple{3,AbstractVector}; overflow=_default_overflow)
    if all(_is_uniform_bins.(edges))
        r = (range(first(edges[1]), last(edges[1]), length=length(edges[1])),
             range(first(edges[2]), last(edges[2]), length=length(edges[2])),
             range(first(edges[3]), last(edges[3]), length=length(edges[3])))
        return Hist3D(A, r; overflow=overflow)
    else
        h = Hist3D(Int; bins=edges, overflow=overflow)
        push!.(h, A[1], A[2], A[3])
        return h
    end
end

"""
    Hist3D(tuple, wgts::AbstractWeights, edges::Tuple{AbstractRange,AbstractRange}; overflow)
    Hist3D(tuple, wgts::AbstractWeights, edges::Tuple{AbstractVector,AbstractVector}; overflow)

Create a `Hist3D` with given bin `edges` and values from
a 2-tuple of arrays of x, y values.
`wgts` should have the same `size` as elements of `tuple`.
"""
function Hist3D(A::NTuple{3,AbstractVector}, wgts::AbstractWeights, r::NTuple{3,AbstractRange}; overflow=_default_overflow)
    @boundscheck @assert size(A[1]) == size(A[2]) == size(wgts)
    h = Hist3D(eltype(wgts); bins=r, overflow=overflow)
    push!.(h, A[1], A[2], wgts)
    return h
end
function Hist3D(A::NTuple{3,AbstractVector}, wgts::AbstractWeights, edges::NTuple{3,AbstractVector}; overflow=_default_overflow)
    if all(_is_uniform_bins.(edges))
        r = (range(first(edges[1]), last(edges[1]), length=length(edges[1])),
             range(first(edges[2]), last(edges[2]), length=length(edges[2])),
             range(first(edges[3]), last(edges[3]), length=length(edges[3])))
        return Hist3D(A, wgts, r; overflow=overflow)
    else
        h = Hist3D(Int; bins=edges, overflow=overflow)
        push!.(h, A[1], A[2], A[3], wgts)
        return h
    end
end

"""
    Hist3D(A::AbstractVector{T}; nbins::Tuple{Integer,Integer}, overflow) where T
    Hist3D(A::AbstractVector{T}, wgts::AbstractWeights; nbins::Tuple{Integer,Integer}, overflow) where T

Automatically determine number of bins based on `Sturges` algo.
"""
function Hist3D(A::NTuple{3,AbstractVector{T}};
        nbins::NTuple{3,Integer}=_sturges.(A),
        overflow=_default_overflow,
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
        overflow=_default_overflow,
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

function _svg(h::Hist3D)
    framewidth = 250
    frameheight = 200

    x1, x2 = 60, 180
    y1, y2 = 30, 130
    z1, z2 = 10, 120

    # cabinet projection
    α = pi/4
    P = [1 0 1/2*cos(α) ; 0 1 1/2*sin(α); 0 0 0]
    projection(points) = map(p->begin
                                 x,y,z = (P*p)
                                 [round(x),round(frameheight-y)]
                             end, points)
    flatstring(points) = join([points...;], " ")

    front = [[x1,y1,z1],[x2,y1,z1],[x2,y2,z1],[x1,y2,z1],[x1,y1,z1]]
    right = [[x2,y1,z2],[x2,y2,z2],[x2,y2,z1],[x2,y1,z1],[x2,y1,z2]]
    top = [[x1,y2,z2],[x1,y2,z1],[x2,y2,z1],[x2,y2,z2],[x1,y2,z2]]

    frontface = projection(front)
    rightface = projection(right)
    topface = projection(top)

    xaxisline = projection([[x1,y1,z1],[x2,y1,z1]])
    yaxisline = projection([[x1,y1,z1],[x1,y2,z1]])
    zaxisline = projection([[x1,y2,z1],[x1,y2,z2]])

    (x1text, x2text), (y1text, y2text), (z1text, z2text) = map(x->(first(x),last(x)), binedges(h))

    c1, c2, c3 = "#1f77b4", "#d62728", "#2ca02c"
    textstyle = """ font-size="85%" font-family="sans-serif" """
    return """
    <svg width="$(framewidth)" height="$(frameheight)" version="1.1" xmlns="http://www.w3.org/2000/svg">
        <polyline points="$(flatstring(frontface))" stroke="black" stroke-width="1" fill="none"/>
        <polyline points="$(flatstring(rightface))" stroke="black" stroke-width="1" fill="none"/>
        <polyline points="$(flatstring(topface))" stroke="black" stroke-width="1" fill="none"/>
        <polyline points="$(flatstring(xaxisline))" stroke="$(c1)" stroke-width="2" fill="none"/>
        <polyline points="$(flatstring(yaxisline))" stroke="$(c2)" stroke-width="2" fill="none"/>
        <polyline points="$(flatstring(zaxisline))" stroke="$(c3)" stroke-width="2" fill="none"/>
        <text x="$(xaxisline[1][1])" y="$(xaxisline[1][2]+5)" fill="$(c1)" dominant-baseline="hanging" text-anchor="left" $(textstyle)>$(x1text)</text>
        <text x="$(xaxisline[2][1])" y="$(xaxisline[2][2]+5)" fill="$(c1)" dominant-baseline="hanging" text-anchor="end" $(textstyle)>$(x2text)</text>
        <text x="$(yaxisline[1][1]-5)" y="$(yaxisline[1][2])" fill="$(c2)" dominant-baseline="text-bottom" text-anchor="end" $(textstyle)>$(y1text)</text>
        <text x="$(yaxisline[2][1]-5)" y="$(yaxisline[2][2])" fill="$(c2)" dominant-baseline="hanging" text-anchor="end" $(textstyle)>$(y2text)</text>
        <text x="$(zaxisline[1][1]+15)" y="$(zaxisline[1][2]-4)" fill="$(c3)" dominant-baseline="text-bottom" text-anchor="start" $(textstyle)>$(z1text)</text>
        <text x="$(zaxisline[2][1]+5)" y="$(zaxisline[2][2]+3)" fill="$(c3)" dominant-baseline="hanging" text-anchor="start" $(textstyle)>$(z2text)</text>
    </svg>
         """
end

function Base.show(io::IO, h::Hist3D)
    print(io, "Hist3D{$(eltype(bincounts(h)))}, ")
    print(io, "edges=$(repr(binedges(h), context=:limit => true)), ")
    print(io, "integral=$(integral(h))")
end

function Base.show(io::IO, m::MIME"text/html", h::Hist3D)
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
