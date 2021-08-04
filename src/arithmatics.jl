import Base: ==, +, -, *, /, merge, merge!

#for StatsBase histogram
for op in (:+, :-, :*, :/)
    @eval ($op)(h1::Histogram, h2::Histogram) = 
    (==)(h1.edges,h2.edges) ? Histogram(h1.edges, broadcast($op, h1.weights, h2.weights)) : throw(DimensionMismatch("The bins of the two histograms don't match"))
    @eval ($op)(h1::Histogram, n::Real) = Histogram(h1.edges, broadcast($op, h1.weights, n))
end

# + -
for op in (:+, :-)
    for T in (:Hist1D,)
        @eval function ($op)(h1::($T), h2::($T))
            h1.hist.edges != h2.hist.edges && throw(DimensionMismatch("Edges don't match"))
            _f(counts) = any(x -> x<0, counts)
            _hist = ($op)(h1.hist,  h2.hist)
            ($T)(_hist, h1.sumw2 + h2.sumw2)
        end
    end
end

# *
for T in (:Hist1D,)
    @eval function *(h1::($T), num::Real)
        _f(counts) = any(x -> x<0, counts)
        _f(h1.hist.weights) && error("Can't do * when bin count is negative")
        _hist = *(h1.hist, num)
        ($T)(_hist, h1.sumw2 * num^2)
    end
end

# /
# https://github.com/aminnj/yahist/blob/4a5767f181ec7fdcc4af18cf15ceedd1c2f89019/yahist/hist1d.py#L427-L430
for T in (:Hist1D,)
    @eval function /(h1::($T), h2::($T))
        _f(counts) = any(x -> x<0, counts)
        (_f(h1.hist.weights) || _f(h2.hist.weights)) && error("Can't do / when bin count is negative")
        _hist = /(h1.hist, h2.hist)
        _sumw2 =  @. h1.sumw2 / (h2.hist.weights ^ 2) +
                (sqrt(h2.sumw2) * h1.hist.weights / (h2.hist.weights ^ 2)) ^ 2
                       
        ($T)(_hist, _sumw2)
    end
end

function merge!(h1::Hist1D, h2::Hist1D)
    h1.hist.edges != h2.hist.edges && throw(DimensionMismatch("The dimension doesn't match"))
    lock(h1)
    h1.hist.weights .+= h2.hist.weights
    h1.sumw2 .+= h2.sumw2
    unlock(h1)
    Hist1D(h1.hist)
end

merge(h1::Hist1D, h2::Hist1D) = merge!(deepcopy(h1), h2)
function merge(hist1ds...)
    h = deepcopy(first(hist1ds))
    length(hist1ds) == 1 && return h
    foreach(x-> merge!(h, x), hist1ds[2:end])
    h
end

function ==(h1::Hist1D, h2::Hist1D)
    h1.hist == h2.hist &&
    h1.sumw2 == h2.sumw2
end
