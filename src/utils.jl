
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
    s = sqrt(n + 0.25)
    s + 0.5, s - 0.5
end
@inline function sqrt_err(n::Real)
    s = sqrt(n)
    s, s
end

_sturges(x) = StatsBase.sturges(length(x))

@inline function _is_uniform_bins(A::AbstractVector{T}) where {T<:Real}
    diffs = diff(A)
    diff1 = first(diffs)
    all(isapprox.(diff1, diffs; atol=1e-9)) #hack
end
function _is_uniform_bins(A::AbstractRange{T}) where {T<:Real}
    true
end


Base.convert(::Type{StatsBase.Histogram}, h::Union{Hist1D,Hist2D,Hist3D}) = StatsBase.Histogram(binedges(h), bincounts(h))

"""
    valid_rebin_values(h::Union{Hist1D, Hist2D, Hist3D})

Calculates the legal values for rebinning, essentially the prime factors of
the number of bins. For a 1D histogram, a `Set` of numbers is return, for higher
dimensional histograms a `Vector{Set}` for each dimension.
"""
valid_rebin_values(h::Hist1D) = _factor(nbins(h))
valid_rebin_values(h::Union{Hist2D,Hist3D}) = [_factor(x) for x in nbins(h)]

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

"""
    chi2, ndf = chi2ndf_ww(h1::Hist1D, h2::Hist1D)

Calculate the chi2 and ndf between two weighted histograms. 
This is useful for comparing two histograms to see if they are statistically compatible. 
    
The chi2 is calculated as ROOT does. For more details see the [ROOT implementation](https://root.cern.ch/doc/master/classTH1.html#ab7d63c7c177ccbf879b5dc31f2311b27).
"""
function chi2ndf_ww(h1::Hist1D, h2::Hist1D)
    if nbins(h1) != nbins(h2)
        error("Histograms have different number of bins")
    end

    chi2 = 0.0
    ndf = nbins(h1) - 1

    bins1 = bincounts(h1)
    bins2 = bincounts(h2)

    # if the sum of squares of weights has been defined (via Sumw2), 
    # this function returns the sqrt(sum of w2). otherwise it returns the sqrt(contents) for this bin. 
    errs1 = sumw2(h1)
    errs2 = sumw2(h2)

    sum1 = sum(bins1)
    sum2 = sum(bins2)

    for (cnt1, cnt2, e1sq, e2sq) in zip(bins1, bins2, errs1, errs2)
        if cnt1 * cnt1 == 0.0 && cnt2 * cnt2 == 0.0
            ndf -= 1
            continue
        end

        if e1sq == 0.0 && e2sq == 0.0
            error("Both errors are zero")
        end

        sigma = sum1 * sum1 * e2sq + sum2 * sum2 * e1sq
        delta = sum2 * cnt1 - sum1 * cnt2
        chi2 += delta * delta / sigma
    end

    return chi2, ndf
end