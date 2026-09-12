const _UniformRange = StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64},Int}
const _EMPTY_RANGE = convert(_UniformRange, 0.0:-1.0)

"""
    BinEdges <: AbstractVector{Float64}

This type implements a vector-like data structure to be used for histogram bin edges, it can
handle both uniform and non-uniform binnings in a single type to reduce the amount of parametric
types.

Bin lookup (`searchsortedlast(edges, x)`) is O(1) whenever the edges are (numerically) uniform,
regardless of whether they were given as an `AbstractRange` or as a plain `Vector`; otherwise a
binary search is used. In both cases the result is exact, i.e. identical to
`searchsortedlast(collect(edges), x)`, including for values that sit exactly on a bin edge.

The edges are always copied on construction, mutating the input afterwards does not affect
the histogram.

!!! note
    Due to the usage of Float64, bin edges shouldn't contain element with absolute value larger than
    9007199254740992, which is the `maxintfloat(Float64)`.
"""
struct BinEdges <: AbstractVector{Float64}
    edges::Vector{Float64}   # materialized edges, always populated
    isuniform::Bool          # O(1) lookup is possible (and exact)
    isrange::Bool            # constructed from an AbstractRange (only affects display)
    range::_UniformRange     # the original range when `isrange`, empty otherwise
    inv_step::Float64
    rfirst::Float64
    rlast::Float64
end

BinEdges(b::BinEdges) = b

function BinEdges(edges::AbstractRange)
    length(edges) >= 2 || throw(ArgumentError("BinEdges must have at least two edges"))
    r = convert(_UniformRange, edges)
    step(r) > 0 || throw(ArgumentError("BinEdges must be strictly increasing"))
    v = collect(r)
    _check_edges(v)
    inv_step = inv(step(r))
    return BinEdges(v, _uniform_lookup_ok(v, inv_step), true, r, inv_step, first(v), last(v))
end

function BinEdges(edges::AbstractVector)
    length(edges) >= 2 || throw(ArgumentError("BinEdges must have at least two edges"))
    v = Vector{Float64}(edges) # always a copy
    for i in 1:length(v)-1
        @inbounds v[i] < v[i+1] || throw(ArgumentError("BinEdges must be strictly increasing (sorted and unique), got $(v[i]) followed by $(v[i+1])"))
    end
    _check_edges(v)
    inv_step = (length(v) - 1) / (last(v) - first(v))
    return BinEdges(v, _uniform_lookup_ok(v, inv_step), false, _EMPTY_RANGE, inv_step, first(v), last(v))
end

function _check_edges(v::Vector{Float64})
    (isfinite(first(v)) && isfinite(last(v))) || throw(ArgumentError("BinEdges must be finite"))
    if abs(first(v)) > maxintfloat(Float64) || abs(last(v)) > maxintfloat(Float64)
        throw(ArgumentError("BinEdges cannot contain element with absolute value larger than $(maxintfloat(Float64))"))
    end
    return nothing
end

# The O(1) lookup computes a guess `g = trunc((x - first) * inv_step) + 1` and then corrects
# it by at most one bin using the real edge values. That is exact as long as, for every edge
# `e_i`, the guess lands on `i` or `i-1`; the guess is monotonic in `x` so every value inside a
# bin is then also within ±1 of the truth.
function _uniform_lookup_ok(v::AbstractVector{F}, inv_step::F) where {F<:AbstractFloat}
    isfinite(inv_step) || return false
    x1 = first(v)
    for i in eachindex(v)
        f = (v[i] - x1) * inv_step
        (0 <= f < length(v) + 1) || return false
        g = unsafe_trunc(Int, f) + 1
        (i - 1 <= g <= i) || return false
    end
    return true
end

isuniform(b::BinEdges) = b.isuniform

Base.size(b::BinEdges) = size(b.edges)
Base.IndexStyle(::Type{BinEdges}) = IndexLinear()
Base.@propagate_inbounds Base.getindex(b::BinEdges, i::Int) = b.edges[i]
Base.first(b::BinEdges) = b.rfirst
Base.last(b::BinEdges) = b.rlast
Base.diff(b::BinEdges) = diff(b.edges)
Base.copy(b::BinEdges) = BinEdges(copy(b.edges), b.isuniform, b.isrange, b.range, b.inv_step, b.rfirst, b.rlast)
Base.:(==)(a::BinEdges, b::BinEdges) = a.edges == b.edges
Base.hash(b::BinEdges, h::UInt) = hash(b.edges, h)

Base.convert(::Type{BinEdges}, edges::AbstractRange) = BinEdges(edges)
Base.convert(::Type{BinEdges}, edges::AbstractVector) = BinEdges(edges)

@inline function _searchsortedlast_uniform(b::BinEdges, x::Float64)
    x < b.rfirst && return 0
    # `x >= last` as well as NaN: same as `searchsortedlast(::Vector, x)`
    x < b.rlast || return length(b.edges)
    i = unsafe_trunc(Int, (x - b.rfirst) * b.inv_step) + 1
    e = b.edges
    @inbounds if x < e[i]
        i -= 1
    elseif x >= e[i+1]
        i += 1
    end
    return i
end

# Branchless binary search (the number of iterations only depends on the length), ~30% faster
# than `searchsortedlast(::Vector, x)` for a few hundred edges and identical results.
@inline function _searchsortedlast_nonuniform(b::BinEdges, x::Float64)
    x < b.rfirst && return 0
    x < b.rlast || return length(b.edges)  # x >= last, or NaN
    v = b.edges
    lo = 0
    len = length(v)
    @inbounds while len > 0
        half = len >>> 1
        mid = lo + half
        c = v[mid+1] <= x
        lo = ifelse(c, mid + 1, lo)
        len = ifelse(c, len - half - 1, half)
    end
    return lo
end

"""
    searchsortedlast(b::BinEdges, x::Real)

Index `i` such that `b[i] <= x < b[i+1]`; `0` if `x < first(b)` and `length(b)` if
`x >= last(b)` (also for `NaN`). O(1) for uniform edges, binary search otherwise.
"""
@inline function Base.searchsortedlast(b::BinEdges, x::Real)
    if isuniform(b)
        return _searchsortedlast_uniform(b, Float64(x))
    else
        return _searchsortedlast_nonuniform(b, Float64(x))
    end
end

Base.show(io::IO, b::BinEdges) = show(io, b.isrange ? b.range : b.edges)
function Base.show(io::IO, mime::MIME"text/plain", b::BinEdges)
    print(io, (b.isrange || _is_uniform_bins(b.edges)) ? "Uniform FHist.BinEdges: " : "Non-uniform FHist.BinEdges: ")
    show(io, b)
end
