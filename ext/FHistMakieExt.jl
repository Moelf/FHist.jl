module FHistMakieExt
using FHist, FHist.Measurements
using Statistics
isdefined(Base, :get_extension) ? (using Makie) : (using ..Makie)

import FHist: stackedhist, stackedhist!

function __init__()
    FHist.ATLASTHEME = Makie.Theme(
        Axis=(
            xtickalign=true, ytickalign=true,
            xticksmirrored=true, yticksmirrored=true,
            xminortickalign=1, yminortickalign=1,
            xticksize=10, yticksize=10,
            xminorticksize=6, yminorticksize=6,
            xgridvisible=false, ygridvisible=false,
            xminorticksvisible=true, yminorticksvisible=true,
        ),
        Colorbar=(
            colormap=:haline,
            highclip=:red,
            lowclip=:black
        )
    )
end

function _clamp_counts!(c_vec)
    min_positive = eps()
    @. c_vec = max(c_vec, min_positive)
    return nothing
end

function _clamp_counts_errors!(c_vec, el_vec, eh_vec)
    # Set the clipping, and make copy of starting counts
    min_positive = eps()
    c_vec_def = copy(c_vec)

    _clamp_counts!(c_vec)

    # clip lower errors
    mask = c_vec - el_vec .< min_positive
    @views el_vec[mask] = c_vec[mask] .- min_positive

    # clip higher errors
    @. eh_vec = max(eh_vec - (c_vec - c_vec_def), min_positive)

    return nothing
end


"""
    stackedhist(hs:AbstractVector{<:Hist1D}; errors=true|:bar|:shade, color=Makie.wong_colors(), gap=-0.01)

Plot a vector of 1D histograms stacked, use `errors` to show or hide error bar in the plot.
`errors = true` and `errors = :bar` are equivalent.

`color` should be a vector of colors that is at least `length(hs)` long. See below example
regarding how to make legends semi-manually.

# Examples
```julia
h1 = Hist1D(randn(1000); binedges = -3:0.3:3)
h2 = Hist1D(randn(10000); binedges = -3:0.3:3)

fig, a, p = stackedhist([h1, h2, h2])
labels = ["ZZ", "Z+jets", "ttbarZ"]
elements = [PolyElement(polycolor = p.attributes.color[][i]) for i in 1:length(labels)]
title = "Processes"

Legend(fig[1,2], elements, labels, title)
fig
```

!!! note
    The `gap` attribute is used to control the gap between the bars, it is set to `-0.01` by default to supress
    artifacts in Cairo backend.
"""
@recipe(StackedHist) do scene
    Attributes(
        error_color=(:black, 0.5),
        color=Makie.wong_colors(),
        labels=nothing,
        whiskerwidth=10,
        gap=-0.01
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
        stack=grp,
        color=c[grp],
        gap=input[:gap],
    )

    error_color = input[:error_color]
    if error_color ∈ (true, :bar)
        errorbars!(input, centers, totals, errs / 2, whiskerwidth=input[:whiskerwidth])
    else
        crossbar!(input, centers, totals, totals .+ errs / 2, totals .- errs / 2;
            gap=input[:gap],
            width=diff(_e),
            show_midline=false,
            color=error_color
        )
    end
    input
end

import FHist: ratiohist, ratiohist!

"""
    ratiohist(h::Hist1D; color=:black, errors=true)

Plot a histogram that represents a ratio (i.e. `h = h1/h3`), you can pass `color` to fix the error bar colors, and use `error` to turn on or off the error bars in the ratio plot.
"""
@recipe(RatioHist) do scene
    Attributes(
        errors=true,
        whiskerwidth=10,
        color=:black
    )
end

function Makie.plot!(input::RatioHist{<:Tuple{<:Hist1D}})
    hratio = input[1][]
    xs = bincenters(hratio)
    ys = bincounts(hratio)

    color = input[:color]

    scatter!(input, xs, ys; color=color)
    if input[:errors][]
        errorbars!(input, xs, ys, binerrors(hratio); color=color, whiskerwidth=input[:whiskerwidth][])
    end
    hlines!(input, 1; color=RGBf(0.2, 0.2, 0.2), linestyle=:dashdot)
    input
end
function Makie.plot!(input::RatioHist{<:Tuple{<:Hist1D,<:Hist1D}})
    hratio = input[1][] / input[2][]
    ratiohist!(input, hratio)
end

Makie.used_attributes(::Type{<:Makie.Plot}, h::Hist1D) = (:clamp_bincounts,)
function Makie.convert_arguments(P::Type{<:Scatter}, h::Hist1D; clamp_bincounts=false)
    ys = copy(bincounts(h))
    if clamp_bincounts
        _clamp_counts!(ys)
    end
    convert_arguments(P, bincenters(h), ys)
end
function Makie.convert_arguments(P::Type{<:BarPlot}, h::Hist1D; clamp_bincounts=false)
    ys = copy(bincounts(h))
    if clamp_bincounts
        _clamp_counts!(ys)
    end
    convert_arguments(P, bincenters(h), ys)
end

Makie.plottype(::Hist1D) = Hist
function Makie.convert_arguments(P::Type{<:Stairs}, h::Hist1D; clamp_bincounts=false)
    edges = binedges(h)
    phantomedge = edges[end] # to bring step back to baseline
    bot = eps()
    bc = copy(bincounts(h))
    if clamp_bincounts
        _clamp_counts!(bc)
    end
    z = zero(eltype(bc))
    nonzero_bincounts = replace(bc, z => bot)
    convert_arguments(P, vcat(edges, phantomedge), vcat(bot, nonzero_bincounts, bot))
end

Makie.used_attributes(::Type{<:Errorbars}, h::Hist1D) = (:clamp_bincounts, :clamp_errors, :error_function)
function Makie.convert_arguments(P::Type{<:Makie.Errorbars}, h::FHist.Hist1D; clamp_bincounts=false, clamp_errors=true, error_function=nothing)
    xs = FHist.bincenters(h)
    ys = copy(FHist.bincounts(h))
    errs = if isnothing(error_function)
        FHist.binerrors(FHist.sqrt, h)
    else
        FHist.binerrors(error_function, h)
    end
    hi_errs, lo_errs = first.(errs), last.(errs)

    if clamp_bincounts && clamp_errors
        _clamp_counts_errors!(ys, lo_errs, hi_errs)
    elseif !clamp_bincounts && clamp_errors
        for i in eachindex(ys, lo_errs)
            if ys[i] - lo_errs[i] <= 0
                lo_errs[i] = ys[i] - eps()
            end
        end
    elseif clamp_bincounts && !clamp_errors
        error("Clamping bincounts without also clamping errors will produce incorrect visualization.")
    end
    convert_arguments(P, xs, ys, lo_errs, hi_errs)
end

function Makie.convert_arguments(P::Type{<:CrossBar}, h::Hist1D)
    cs = bincounts(h)
    es = binerrors(h)
    convert_arguments(P, bincenters(h), cs, cs .- es / 2, cs .+ es / 2)
end

function Makie.plot!(plot::Hist{<:Tuple{<:Hist1D}})
    attributes = Makie.Attributes(plot)
    myhist = plot[1]
    barplot!(plot, attributes, myhist; fillto=eps(), width=diff(binedges(myhist[])))
    plot
end

function Makie.plot!(plot::StepHist{<:Tuple{<:Hist1D}})
    valid_attributes = Makie.Attributes(plot)
    stairs!(plot, valid_attributes, plot[1])
    plot
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
function FHist.statbox!(fig::Makie.FigureAxisPlot, h; position=(1, 2))
    f, _, _ = fig
    statbox!(f, h; position)
    fig
end
function FHist.statbox!(fig::Makie.Figure, h::Hist1D; position=(1, 2))
    N = nentries(h)
    M = round(mean(h); sigdigits=2)
    S = round(std(h); sigdigits=2)
    labels = ["Entries = $N", "Mean = $M", "Std Dev = $S", "Overflow = $(h.overflow)"]
    elements = fill(PolyElement(polycolor=:transparent), 4)
    Legend(getindex(fig, position...), elements, labels)
    fig
end
function FHist.statbox!(fig::Makie.Figure, h::Hist2D; position=(1, 2))
    N = nentries(h)
    xM, yM = round.(mean(h); sigdigits=2)
    xS, yS = round.(std(h); sigdigits=2)
    labels = ["Entries = $N", "Mean x = $xM", "Mean y = $yM", "Std Dev x = $xS", "Std Dev y = $yS", "Overflow = $(h.overflow)"]
    elements = fill(PolyElement(polycolor=:transparent), 6)
    Legend(getindex(fig, position...), elements, labels)
    fig
end

Makie.plottype(::Hist2D) = Heatmap
function Makie.convert_arguments(P::Type{<:Heatmap}, h2d::Hist2D)
    counts = bincounts(h2d)
    z = zero(eltype(counts))
    convert_arguments(P, binedges(h2d)..., replace(counts, z => NaN))
end

function Makie.convert_arguments(::VertexGrid, h2d::Hist2D)
    counts = bincounts(h2d)
    z = zero(eltype(counts))
    (bincenters(h2d)..., replace(counts, z => NaN))
end

_to_endpoints(binedge) = (first(binedge), last(binedge))

Makie.plottype(::Hist3D) = Volume
function Makie.convert_arguments(P::Type{<:Volume}, h::Hist3D)
    counts = bincounts(h)
    z = zero(eltype(counts))
    convert_arguments(P, _to_endpoints.(binedges(h))..., replace(counts, z => NaN))
end

"""
    collabtext!(axis, colabname = "ATLAS", stage = "Preliminary"; position::Union{Symbol, Point2f} = :lt)

Inject collaboration text such as `ATLAS/CMS Preliminary` into the plot. The position `Point2f` is in relative x and y.

## Example
```julia
h1 = Hist1D(randn(10^4))
with_theme(ATLASTHEME) do
    fig, ax, p = tairs(h1)
    errorbars!(h1)
    collabtext!(ax)
    fig
end
```
"""
function FHist.collabtext!(axis, colabname="ATLAS", stage="Preliminary"; position=:lt)
    relative_projection = Makie.camrelative(axis.scene)
    pos = if position isa Symbol
        length(String(position)) != 2 && throw("`position` must be length == 2, support `lt` or `rt`")
        position == :lt ? Point2f(0.04, 0.94) : Point2f(0.70, 0.94)
    else
        position
    end
    text!(relative_projection, "$colabname $stage", position=pos,
        font=[fill("TeX Gyre Heros Bold Italic Makie", length(colabname)); fill("TeX Gyre Heros Makie", length(stage) + 1)]
    )
end

end
