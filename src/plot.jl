@recipe function f(h::Hist1D)
    seriestype --> :barbins
    x:=h.hist.edges[1]
    y:=h.hist.weights
    ()
end

function MakieCore.convert_arguments(P::Type{<:AbstractPlot}, h::Hist1D)
    @show P
    MakieCore.convert_arguments(P, rand(1:5, 2,3))
end
