mutable struct Hist1D{T<:Real, E} <: AbstractHistogram{T,1,E}
    hist::Histogram{T,1,E}
    errors_up::Vector{Float64}
    errors_down::Vector{Float64}
    error_mode::Symbol #default is :sqrt
    # most concrete inner constructor
    function Hist1D(h::Histogram{T,1,E}, errors_up, errors_down, error_mode=:sqrt) where {T,E}
        h.closed == :left || error("Hist1D must be only closed on :left")
        size(errors_up) == size(errors_down) == size(h.weights) || error("Bin shapes don't match")
        new{T,E}(h, errors_up, errors_down, error_mode)
    end
end


# add error to StatsBase histogram
function Hist1D(h::Histogram; error_mode=:sqrt)
    e_up, e_down = _make_error(h.weights, error_mode)
    Hist1D(h, e_up, e_down, error_mode)
end

# TODO this is too confusing, don't include
# function Hist1D(edges, weights; error_mode=:sqrt)
#     h1 = Histogram(edges, weights)
#     Hist1D(h1; error_mode = error_mode)
# end

# fast uniform
function Hist1D(A, r::AbstractRange; kwgs...)
    s = step(r)
    start = first(r)
    start2 = start+0.5s
    stop = last(r)
    L = length(r) - 1
    counts = zeros(Int, L)
    @inbounds for i in A
        # skip overflow
        c = ifelse(i < start || i > stop, 0, 1)
        id = round(Int, (i-start2)/s) + 1
        counts[clamp(id,1,L)] += c
    end
    Hist1D(Histogram(r, counts); kwgs...)
end

function Hist1D(A, wgts::AbstractWeights, r::AbstractRange, ; kwgs...)
    @assert size(A) == size(wgts)
    s = step(r)
    start = first(r)
    start2 = start+0.5s
    stop = last(r)
    L = length(r) - 1
    counts = zeros(L)
    @inbounds for i in eachindex(A)
        # skip overflow
        c = ifelse(A[i] < start || A[i] > stop, 0.0, wgts[i])
        id = round(Int, (A[i]-start2)/s) + 1
        counts[clamp(id,1,L)] += c
    end
    Hist1D(Histogram(r, counts); kwgs...)
end

# fit then wrap
function Hist1D(A, edges::AbstractVector; error_mode=:sqrt)
    if _is_uniform_bins(edges)
        s = edges[begin+1] - first(edges)
        r = first(edges):s:last(edges)
        Hist1D(A, r; error_mode = error_mode)
    else
        Hist1D(fit(Histogram, A, edges); error_mode=error_mode)
    end
end

function Hist1D(A, wgts::AbstractWeights, edges::AbstractVector; error_mode=:sqrt)
    if _is_uniform_bins(edges)
        s = edges[begin+1] - first(edges)
        r = first(edges):s:last(edges)
        Hist1D(A, wgts, r; error_mode = error_mode)
    else
        Hist1D(fit(Histogram, A, wgts, edges); error_mode=error_mode)
    end
end

function Hist1D(A::AbstractVector{T}; nbins::Integer=StatsBase.sturges(length(A)), error_mode=:sqrt) where T
    F = float(T)
    lo,hi = extrema(A)
    r = StatsBase.histrange(F(lo), F(hi), nbins)
    Hist1D(A, r; error_mode=error_mode)
end

function Hist1D(A::AbstractVector{T}, wgts::AbstractWeights; nbins::Integer=StatsBase.sturges(length(A)), error_mode=:sqrt) where T
    F = float(T)
    lo,hi = extrema(A)
    r = StatsBase.histrange(F(lo), F(hi), nbins)
    Hist1D(A, wgts, r; error_mode=error_mode)
end

function Base.show(io::IO, h::Hist1D)                                        
    # println(io, typeof(h))
    show(io, h.hist)
    println()
    println(io,"errors: ")                                                   
    println(io,"  up  : ",h.errors_up)                                         
    println(io,"  down: ",h.errors_down)                                       
    print(io,"error_mode: ",h.error_mode)                                      
end                                                                          
