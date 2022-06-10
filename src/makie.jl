using .Makie

"""

## Example
```
with_theme(ATLASTHEME) do
    h1 = Hist1D(randn(10^4))
    hist(h1; label="atlas style histogram")
end
```
"""
const ATLASTHEME = 
Theme(
      Axis = (
              xtickalign=1, ytickalign=1, 
              xminortickalign=1, yminortickalign=1,
              xticksize=10, yticksize=10,
              xminorticksize=6, yminorticksize=6,
              xgridvisible = false, ygridvisible = false,
              xminorticksvisible = true, yminorticksvisible = true,
              xticks = WilkinsonTicks(6), yticks = WilkinsonTicks(6),
              xminorticks = IntervalsBetween(5), yminorticks = IntervalsBetween(5),
             )
     )

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
Makie.convert_arguments(P::Type{<:Errorbars}, h::Hist1D) = convert_arguments(P, bincenters(h), bincounts(h), binerrors(h)/2)
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

"""
    statbox!(fig::Uiont{Figure, AxisFigurePlot}, h::Union{Hist1D, Hist2D}; position = (1,2))

Add a CERN ROOT style statbox to an existing figure.

##Example
```julia
h1 = Hist1D(randn(10^4))
afp = hist(h1; label="a")
statbox!(afp, h1)
```
"""
function statbox!(fig::Makie.FigureAxisPlot, h; position = (1,2))
    f, _, _ = fig
    statbox!(f, h; position)
    fig
end
function statbox!(fig::Makie.Figure, h::Hist1D; position = (1,2))
    N = nentries(h)
    M = round(mean(h); sigdigits= 2)
    S = round(std(h); sigdigits= 2)
    labels = ["Entries = $N", "Mean = $M", "Std Dev = $S", "Overflow = $(h.overflow)" ]
    elements = fill(PolyElement(polycolor = :transparent), 4)
    Legend(getindex(fig, position...), elements, labels)
    fig
end
function statbox!(fig::Makie.Figure, h::Hist2D; position = (1,2))
    N = nentries(h)
    xM, yM = round.(mean(h); sigdigits= 2)
    xS, yS = round.(std(h); sigdigits= 2)
    labels = ["Entries = $N", "Mean x = $xM", "Mean y = $yM", "Std Dev x = $xS", "Std Dev y = $yS", "Overflow = $(h.overflow)" ]
    elements = fill(PolyElement(polycolor = :transparent), 6)
    Legend(getindex(fig, position...), elements, labels)
    fig
end

Makie.MakieCore.plottype(::Hist2D) = Heatmap
function Makie.convert_arguments(P::Type{<:Heatmap}, h2d::Hist2D)
    counts = bincounts(h2d)
    z = zero(eltype(counts))
    convert_arguments(P, bincenters(h2d)..., replace(counts, z => NaN))
end
