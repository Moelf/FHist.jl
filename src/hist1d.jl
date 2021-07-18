struct Hist1D{T<:Real,E} <: AbstractHistogram{T,1,E}
    hist::Histogram{T,1,E}
    errors_up::Vector{Float64}
    errors_down::Vector{Float64}
    error_mode::Symbol #default is :sqrt
    hlock::SpinLock
    # most concrete inner constructor
    function Hist1D(
        h::Histogram{T,1,E}, errors_up, errors_down, error_mode=:sqrt
    ) where {T,E}
        h.closed == :left || error("Hist1D must be only closed on :left")
        size(errors_up) == size(errors_down) == size(h.weights) ||
            error("Bin shapes don't match")
        return new{T,E}(h, errors_up, errors_down, error_mode, SpinLock())
    end
end
Base.lock(h::Hist1D) = lock(h.hlock)
Base.unlock(h::Hist1D) = unlock(h.hlock)

"""
    update_error!(h::Hist1D, error_fun = sqrt_err)

Update the error (up and down) of a histogram according to a specific error_mode.
Remember to call this function after updaing a histogram with `push!()` in a loop.
"""
function update_error!(h::Hist1D, error_fun=sqrt_err)
    lock(h)
    _make_error!(error_fun, h.hist.weights, h.errors_up, h.errors_down)
    unlock(h)
    return h
end

"""
    push!(h::Hist1D, val::Real, wgt::Real=one{T})

Adding one value at a time into histogram. Remember to call [`update_error!`](@ref) after
if you need errors.
"""
function Base.push!(h::Hist1D{T,E}, val::Real, wgt::Real=one(T)) where {T,E}
    @inbounds binidx = searchsortedlast(h.hist.edges[1], val)
    lock(h)
    @inbounds h.hist.weights[binidx] += wgt
    unlock(h)
    return h
end

"""
    Hist1D(elT::Type{T}=Float64; binedges) where {T}

Initialize an empty histogram with bin content typed as `T` and bin edges.
To be used with [`push!`](@ref)
"""
function Hist1D(elT::Type{T}=Float64; bins) where {T}
    counts = zeros(elT, length(bins) - 1)
    e_up = similar(counts, Float64)
    e_down = similar(counts, Float64)
    return Hist1D(Histogram(bins, counts), e_up, e_down)
end

"""
    Hist1D(h::Histogram{T, 1, E}; error_mode=:sqrt) where {T,E}

Convert an existing 1D `StatsBase.Histogram` to a `Hist1D`. Adds
error according to `error_mode`.
"""
function Hist1D(h::Histogram{T,1,E}; error_mode=:sqrt) where {T,E}
    e_up, e_down = _make_error(h.weights, error_mode)
    return Hist1D(h, e_up, e_down, error_mode)
end

"""
    Hist1D(array, edges::AbstractRange; kwgs...)
    Hist1D(array, edges::AbstractVector; error_mode=:sqrt)

Create a `Hist1D` with given bin `edges` and vlaues from
array. Weight for each value is assumed to be 1.
"""
function Hist1D(A, r::AbstractRange; kwgs...)
    s = step(r)
    start = first(r)
    start2 = start + 0.5s
    stop = last(r)
    L = length(r) - 1
    counts = zeros(Int, L)
    @inbounds for idx in eachindex(A)
        # skip overflow
        i = A[idx]
        c = ifelse(i > stop, 0, 1)
        c = ifelse(i < start, 0, c)
        id = round(Int, (i - start2) / s) + 1
        counts[clamp(id, 1, L)] += c
    end
    return Hist1D(Histogram(r, counts); kwgs...)
end
function Hist1D(A, edges::AbstractVector; error_mode=:sqrt)
    if _is_uniform_bins(edges)
        s = edges[2] - first(edges)
        r = first(edges):s:last(edges)
        Hist1D(A, r; error_mode=error_mode)
    else
        Hist1D(fit(Histogram, A, edges); error_mode=error_mode)
    end
end

"""
    Hist1D(array, wgts::AbstractWeights, edges::AbstractRange, ; kwgs...)
    Hist1D(array, wgts::AbstractWeights, edges::AbstractVector; error_mode=:sqrt)

Create a `Hist1D` with given bin `edges` and vlaues from
array. `wgts` should have the same `size` as `array`.
"""
function Hist1D(A, wgts::AbstractWeights, r::AbstractRange, ; kwgs...)
    @boundscheck @assert size(A) == size(wgts)
    s = step(r)
    start = first(r)
    start2 = start + s / 2
    stop = last(r)
    wgt_zero = zero(eltype(wgts))
    L = length(r) - 1
    counts = zeros(L)
    @inbounds for i in eachindex(A)
        # skip overflow
        c = ifelse(A[i] < start || A[i] > stop, wgt_zero, wgts[i])
        id = round(Int, (A[i] - start2) / s) + 1
        counts[clamp(id, 1, L)] += c
    end
    return Hist1D(Histogram(r, counts); kwgs...)
end
function Hist1D(A, wgts::AbstractWeights, edges::AbstractVector; error_mode=:sqrt)
    @inbounds if _is_uniform_bins(edges)
        s = edges[2] - first(edges)
        r = first(edges):s:last(edges)
        Hist1D(A, wgts, r; error_mode=error_mode)
    else
        Hist1D(fit(Histogram, A, wgts, edges); error_mode=error_mode)
    end
end

"""
    Hist1D(A::AbstractVector{T}; nbins::Integer=StatsBase.sturges(length(A)), error_mode=:sqrt) where T
    Hist1D(A::AbstractVector{T}, wgts::AbstractWeights; nbins::Integer=StatsBase.sturges(length(A)), error_mode=:sqrt) where T

Automatically determine number of bins based on `Sturges` algo.
"""
function Hist1D(
    A::AbstractVector{T}; nbins::Integer=StatsBase.sturges(length(A)), error_mode=:sqrt
) where {T}
    F = float(T)
    lo, hi = extrema(A)
    r = StatsBase.histrange(F(lo), F(hi), nbins)
    return Hist1D(A, r; error_mode=error_mode)
end

function Hist1D(
    A::AbstractVector{T},
    wgts::AbstractWeights;
    nbins::Integer=StatsBase.sturges(length(A)),
    error_mode=:sqrt,
) where {T}
    F = float(T)
    lo, hi = extrema(A)
    r = StatsBase.histrange(F(lo), F(hi), nbins)
    return Hist1D(A, wgts, r; error_mode=error_mode)
end

function Base.show(io::IO, h::Hist1D)
    show(io, UnicodePlots.histogram(h.hist; width=30))
    println()
    println(io, "edges: ", h.hist.edges[1])
    println(io, "bin counts: ", h.hist.weights)
    println(io, "errors: ")
    println(io, "  up  : ", round.(h.errors_up; sigdigits=3))
    println(io, "  down: ", round.(h.errors_down; sigdigits=3))
    return print(io, "error_mode: ", h.error_mode)
end
