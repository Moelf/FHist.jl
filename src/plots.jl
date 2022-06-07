using RecipesBase, Statistics

@recipe function f(h::Hist1D)
    seriestype --> :barbins
    I = round(integral(h); sigdigits= 2)
    M = round(mean(h); sigdigits= 2)
    S = round(std(h); sigdigits= 2)
    label --> "Integral = $I\nMean = $M\nStd Dev = $S\nOverflow = $(h.overflow)"
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

