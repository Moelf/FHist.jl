
"""
    hists_to_bars(hist1ds)

Given a vector of Hist1D, return `edges` (xs), `heights` (ys),
and `grps` (for grouping) that is useful for plotting stacked
histogram.
"""
function hists_to_bars(hist1ds)
    L = length(hist1ds)
    oneedge = binedges(hist1ds[1])[1:end-1]
    edges = repeat(oneedge, L)
    heights = mapreduce(bincounts, vcat, hist1ds)
    grps = repeat(1:L; inner=length(oneedge))
    
    edges, heights, grps
end

@inline function pearson_err(n::Real)
    s = sqrt(n+0.25)
    s+0.5, s-0.5
end
@inline function sqrt_err(n::Real)
    s = sqrt(n)
    s,s
end

_sturges(x) = StatsBase.sturges(length(x))

@inline function _is_uniform_bins(A::AbstractVector{T}) where T<:Real
    diffs = diff(A)
    diff1 = first(diffs)
    all(isapprox.(diff1, diffs; atol = 1e-9)) #hack
end
function _is_uniform_bins(A::AbstractRange{T}) where T<:Real
    true
end

"""
    valid_rebin_values(h::Union{Hist1D, Hist2D, Hist3D})

Calculates the legal values for rebinning, essentially the prime factors of
the number of bins. For a 1D histogram, a `Set` of numbers is return, for higher
dimensional histograms a `Vector{Set}` for each dimension.
"""
valid_rebin_values(h::Hist1D) = _factor(nbins(h))
valid_rebin_values(h::Union{Hist2D, Hist3D}) = [_factor(x) for x in nbins(h)]

"""
    function _factor(n::Integer)

Helper function to calculate the prime factors of a given integer.
"""
function _factor(n::Integer)
    factors = Set{Int}()
    limit = n
    factor = 2
    while n > 1 & factor < limit
        while n % factor == 0
            n /= factor
            push!(factors, factor)
            limit = sqrt(n)
        end
        factor += 1
    end
    factors
end
