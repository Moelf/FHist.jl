using .Makie

"""
    stackedhist(hs:AbstractVector{<:Hist1D}; errors=true, color=Makie.wong_colors())

Plot a vector of 1D histograms stacked, use `errors` to show or hide error bar in the plot.

`color` should be a vector of colors that is at least `length(hs)` long. See below example
regarding how to make legends semi-manually.

# Examples
h1 = Hist1D(randn(1000), -3:0.3:3)
h2 = Hist1D(randn(10000), -3:0.3:3)

fig, a, p = stackedhist([h1, h2, h2])
labels = ["ZZ", "Z+jets", "ttbarZ"]
elements = [PolyElement(polycolor = p.attributes.color[][i]) for i in 1:length(labels)]
title = "Processes"

Legend(fig[1,2], elements, labels, title)
fig
"""
@recipe(StackedHist) do scene
    Attributes(
        errors = true,
        color = Makie.wong_colors(),
        labels = nothing,
        whiskerwidth = 15
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
    
    c = input[:color][]
    length(c) < Nhist && throw("provided $(length(c)) colors, not enough for $Nhist histograms")
    Makie.barplot!(input, xs, ys;
        stack = grp,
        color = c[grp],
        gap = 0,
    )
    
    if input[:errors][]
        errorbars!(input, centers, totals, errs/2, whiskerwidth = input[:whiskerwidth][])
    end
    input
end

Makie.MakieCore.plottype(::Hist1D) = Hist
Makie.convert_arguments(P::Type{<:Stairs}, h::Hist1D) = convert_arguments(P, binedges(h), vcat(0.0, bincounts(h)))
Makie.convert_arguments(P::Type{<:BarPlot}, h::Hist1D) = convert_arguments(P, bincenters(h), bincounts(h))
function Makie.plot!(input::Hist{<:Tuple{<:Hist1D}})
    h = input[1][]
    L = haskey(input, :label) ? input[:label][] : nothing
    Makie.barplot!(input, h; gap = 0, label = L)
    input
end

# Makie.MakieCore.plottype(::Hist2D) = Heatmap
# MakieCore.convert_arguments(P::Type{<:Heatmap}, h::Hist1D) = convert_arguments(p, ...)
