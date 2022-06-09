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
    label = haskey(input, :label) ? input[:label][] : nothing
    C = input[:color][]
    # this is a hack, it seems always has non-empty color
    color = C == Makie.RGBA{Float32}(0.0f0,0.0f0,0.0f0,0.6f0) ? Makie.RGB(0/255, 114/255, 178/255) : C
    Makie.barplot!(input, h; gap = 0, label, color=color)
    input
end

function statbox!(fig::Makie.FigureAxisPlot, h)
    f, _, _ = fig
    statbox!(f, h)
    fig
end

function statbox!(fig::Makie.Figure, h)
    N = nentries(h)
    M = round(mean(h); sigdigits= 2)
    S = round(std(h); sigdigits= 2)
    labels = ["Entries = $N", "Mean = $M", "Std Dev = $S", "Overflow = $(h.overflow)" ]
    elements = fill(PolyElement(polycolor = :transparent), 4)
    Legend(fig[1,2], elements, labels)
    fig
end

Makie.MakieCore.plottype(::Hist2D) = Heatmap
Makie.convert_arguments(P::Type{<:Heatmap}, h2d::Hist2D) = convert_arguments(P, bincenters(h2d)..., bincounts(h2d))
