using RecipesBase

@recipe function f(h::Hist1D)
    seriestype --> :barbins
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

