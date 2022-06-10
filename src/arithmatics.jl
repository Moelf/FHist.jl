import Base: ==, +, -, *, /, merge, merge!

#for StatsBase histogram
for op in (:+, :-, :*, :/)
    @eval ($op)(h1::Histogram, h2::Histogram) = 
    (==)(h1.edges,h2.edges) ? Histogram(h1.edges, broadcast($op, h1.weights, h2.weights)) : throw(DimensionMismatch("The bins of the two histograms don't match"))
    @eval ($op)(h1::Histogram, n::Real) = Histogram(h1.edges, broadcast($op, h1.weights, n))
end

for T in (:Hist1D,:Hist2D,:Hist3D)
    for op in (:+, :-)
        @eval function ($op)(h1::($T), h2::($T))
            h1.hist.edges != h2.hist.edges && throw(DimensionMismatch("Edges don't match"))
            h1.overflow != h2.overflow && throw("Can't $op histograms with different overflow settings.")
            _f(counts) = any(x -> x<0, counts)
            _hist = ($op)(h1.hist,  h2.hist)
            ($T)(_hist, h1.sumw2 + h2.sumw2, nentries(h1) + nentries(h2); overflow = h1.overflow)
        end
    end

    @eval function *(h1::($T), num::Real)
        _f(counts) = any(x -> x<0, counts)
        _f(h1.hist.weights) && error("Can't do * when bin count is negative")
        _hist = *(h1.hist, num)
        ($T)(_hist, h1.sumw2 * num^2)
    end

    # https://github.com/aminnj/yahist/blob/4a5767f181ec7fdcc4af18cf15ceedd1c2f89019/yahist/hist1d.py#L427-L430
    @eval function /(h1::($T), h2::($T))
        _f(counts) = any(x -> x<0, counts)
        (_f(h1.hist.weights) || _f(h2.hist.weights)) && error("Can't do / when bin count is negative")
        _hist = /(h1.hist, h2.hist)
        _sumw2 =  @. h1.sumw2 / (h2.hist.weights ^ 2) +
                (sqrt(h2.sumw2) * h1.hist.weights / (h2.hist.weights ^ 2)) ^ 2
                       
        ($T)(_hist, _sumw2)
    end

    @eval function ==(h1::$T, h2::$T)
        h1.hist == h2.hist &&
        h1.sumw2 == h2.sumw2
    end

    @eval function merge!(h1::$T, h2::$T)
        h1.hist.edges != h2.hist.edges && throw(DimensionMismatch("The dimension doesn't match"))
        lock(h1)
        h1.hist.weights .+= h2.hist.weights
        h1.sumw2 .+= h2.sumw2
        unlock(h1)
        ($T)(h1.hist)
    end

    @eval merge(h1::$T, h2::$T) = merge!(deepcopy(h1), h2)
end

function merge(hists...)
    h = deepcopy(first(hists))
    length(hists) == 1 && return h
    foreach(x-> merge!(h, x), hists[2:end])
    h
end
