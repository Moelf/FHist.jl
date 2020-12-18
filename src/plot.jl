@recipe function f(h::Hist1D)
    seriestype --> :barbins
    x:=h.hist.edges[1]
    y:=h.hist.weights
    yerror:=(h.errors_up, h.errors_down)
    ()
end
