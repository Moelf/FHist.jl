using .Makie

"""
    stackedhist(hs:AbstractVector{<:Hist1D}; errors=true)

Plot a vector of 1D histograms stacked, use `errors` to show or hide error bar in the plot.
"""
@recipe(StackedHist) do scene
    Attributes(
        errors = true
    )
end

function Makie.plot!(input::StackedHist{<:Tuple{AbstractVector{<:Hist1D}}})
    hs = input[1][]
    Nhist = length(hs)
    _e = binedges(first(hs))
    all(==(_e), binedges.(hs)) || throw("binedges must match in stacked histogram")

    centers = bincenters(first(hs))
    Nbin = length(centers)
    xs = repeat(centers; outer=Nhist)
    ys = mapreduce(bincounts, vcat, hs)
    grp = repeat(eachindex(hs); inner=Nbin)
    mes = mapreduce(h -> bincounts(h) .± binerrors(h), (.+), hs)
    totals = Measurements.value.(mes)
    errs = Measurements.uncertainty.(mes)
    
    Makie.barplot!(input, xs, ys;
        stack = grp,
        color = grp,
        gap = 0
    )
    
    if input.errors[]
        errorbars!(input, centers, totals, errs/2, whiskerwidth = 15)
    end
    input
end

Makie.MakieCore.plottype(::Hist1D) = Hist
Makie.convert_arguments(P::Type{<:Stairs}, h::Hist1D) = convert_arguments(P, binedges(h), vcat(0.0, bincounts(h)))
Makie.convert_arguments(P::Type{<:BarPlot}, h::Hist1D) = convert_arguments(P, bincenters(h), bincounts(h))
function Makie.plot!(input::Hist{<:Tuple{<:Hist1D}})
    h = input[1][]
    Makie.barplot!(input, h; gap = 0)
    input
end

# Makie.MakieCore.plottype(::Hist2D) = Heatmap
# MakieCore.convert_arguments(P::Type{<:Heatmap}, h::Hist1D) = convert_arguments(p, ...)
