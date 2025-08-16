import Base: ==, +, -, *, /, merge, merge!

for T in (:Hist1D,:Hist2D,:Hist3D)
    for op in (:+, :-)
        @eval function ($op)(h1::($T), h2::($T))
            edge1 = h1.binedges
            edge1 != h2.binedges && throw(DimensionMismatch("Binedges don't match"))
            h1.overflow != h2.overflow && throw("Can't $op histograms with different overflow settings.")
            newcounts = broadcast($op, bincounts(h1),  bincounts(h2))

            ($T)(; binedges = copy.(edge1), bincounts = newcounts, sumw2 = sumw2(h1) + sumw2(h2), nentries = nentries(h1) + nentries(h2), overflow = h1.overflow)
        end
    end

    @eval function *(h1::($T), num::Real)
        any(<(0), bincounts(h1)) && error("Can't scale (*) a histogram when some bin count is negative")
        newcounts = bincounts(h1) * num

        ($T)(; bincounts = newcounts, binedges = copy.(binedges(h1)), sumw2 = sumw2(h1) * num^2, nentries = nentries(h1), overflow = h1.overflow)
    end
    @eval *(num::Real, h1::($T)) = h1 * num

    # https://github.com/aminnj/yahist/blob/4a5767f181ec7fdcc4af18cf15ceedd1c2f89019/yahist/hist1d.py#L427-L430
    @eval function /(h1::($T), h2::($T))
        _f(counts) = any(x -> x<0, counts)
        counts1 = bincounts(h1)
        counts2 = bincounts(h2)
        edge1 = binedges(h1)
        edge1 != binedges(h2) && throw(DimensionMismatch("Binedges don't match in h1/h2"))
        (_f(counts1) || _f(counts2)) && error("Can't divide (/) when some bin counts are negative")
        h1.overflow != h2.overflow && throw("Can't divide two histograms with different overflow settings")

        newcounts = counts1 ./ counts2
        _sumw2 =  sumw2(h1) ./ (counts2 .^ 2) .+
            (sqrt.(sumw2(h2)) .* counts1 ./ (counts2 .^ 2)) .^ 2
                       
        ($T)(bincounts = newcounts, binedges = edge1, sumw2 = _sumw2, nentries = nentries(h1); overflow=h1.overflow)
    end

    @eval function merge!(h1::$T, h2::$T)
        edge1 = h1.binedges
        edge1 != h2.binedges && throw(DimensionMismatch("The dimension doesn't match"))
        lock(h1)
        bincounts(h1) .+= bincounts(h2)
        sumw2(h1) .+= sumw2(h2)
        h1.nentries[] += nentries(h2)
        unlock(h1)
        h1
    end

    @eval merge(h1::$T, h2::$T) = merge!(deepcopy(h1), h2)

    @eval function merge(hists::$T...)
        h = deepcopy(first(hists))
        length(hists) == 1 && return h
        foreach(x-> merge!(h, x), hists[2:end])
        h
    end
end

"""
    significance(signal, bkg) -> `(significance, error_on_significance)`

Calculate the significance of signal vs. bkg histograms, this function uses a more
accurate algorithm than the naive `S / sqrt(B)`

Ref: https://cds.cern.ch/record/2736148/files/ATL-PHYS-PUB-2020-025.pdf

## Example:
```julia
h1 = Hist1D(rand(1000);  binedges = [0, 0.5])
h2 = Hist1D(rand(10000); binedges = [0, 0.5]);

julia> s1 = significance(h1,h2)
(6.738690967342175, 0.3042424717261312)
```
"""
function significance(signal, bkg)
    S = integral(signal)
    B = integral(bkg)
    Sig = sqrt(2*((S + B) * log(1 + S/B) - S))
    dS = sqrt(sum(sumw2(signal)))
    dB = sqrt(sum(sumw2(bkg)))
    dSigdS = log(1 + S/B) / Sig
    dSigdB = (log(1 + S/B) - S/B) / Sig
    err = sqrt((dSigdS * dS)^2 + (dSigdB * dB)^2)
    return Sig, err
end
