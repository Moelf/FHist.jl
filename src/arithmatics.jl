import Base: ==, +, -, *, /, merge, merge!

#for StatsBase histogram
for op in (:+, :-, :*, :/)
    @eval ($op)(h1::Histogram,h2::Histogram) = 
    (==)(h1.edges,h2.edges) ? Histogram(h1.edges, broadcast($op, h1.weights, h2.weights)) : throw(DimensionMismatch("The bins of the two histograms don't match"))
    @eval ($op)(h1::Histogram,n::Real) = Histogram(h1.edges, broadcast($op, h1.weights, n))
end
for op in (:+, :-, :*, :/)
    for T in (:Hist1D,)
        @eval function ($op)(h1::($T), h2::($T))
            h1.hist.edges != h2.hist.edges && throw(DimensionMismatch("Edges don't match"))
            h1.error_mode != h2.error_mode && error("Error mode doesn't exist")
            h3 = ($op)(h1.hist,  h2.hist)
            ($T)(h3; error_mode=h1.error_mode)
        end
    end
end

function merge!(h1::Hist1D, h2::Hist1D)
    h1.hist.edges != h2.hist.edges && throw(DimensionMismatch("The dimension doesn't match"))
    h1.error_mode != h2.error_mode && error("Error mode doesn't exit")

    h1.hist.weights .+= h2.hist.weights
    Hist1D(h1.hist; error_mode = h1.error_mode)
end

merge(h1::Hist1D, h2::Hist1D) = merge!(deepcopy(h1), h2)

function ==(h1::Hist1D, h2::Hist1D)
    h1.hist == h2.hist &&
    h1.errors_up == h2.errors_up &&
    h1.errors_down == h2.errors_down
end
