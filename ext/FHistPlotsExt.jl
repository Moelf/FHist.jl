module FHistPlotsExt
using RecipesBase, Statistics

isdefined(Base, :get_extension) ? (using Plots) : (using ..Plots)

#
@recipe function f(h::Hist1D)
    seriestype --> :barbins
    N = nentries(h)
    M = round(mean(h); sigdigits= 2)
    S = round(std(h); sigdigits= 2)
    label --> "Entries = $N\nMean = $M\nStd Dev = $S\nOverflow = $(h.overflow)"
    x:=h.hist.edges[1]
    y:=h.hist.weights
    ()
end

@recipe function f(h::Hist2D)
    seriestype --> :bins2d
    x := h.hist.edges[1]
    y := h.hist.edges[2]
    z := (surf = h.hist.weights, )
    ()
end
end

