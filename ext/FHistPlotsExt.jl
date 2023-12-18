module FHistPlotsExt
using FHist, RecipesBase, Statistics

isdefined(Base, :get_extension) ? (using Plots) : (using ..Plots)

#
@recipe function f(h::Hist1D)
    seriestype --> :barbins
    N = nentries(h)
    M = round(mean(h); sigdigits= 2)
    S = round(std(h); sigdigits= 2)
    label --> "Entries = $N\nMean = $M\nStd Dev = $S\nOverflow = $(h.overflow)"
    x:= binedges(h)
    y:= bincounts(h)
    ()
end

@recipe function f(h::Hist2D)
    seriestype --> :bins2d
    x := binedges(h)[1]
    y := binedges(h)[2]
    z := (surf = bincounts(h), )
    ()
end
end

