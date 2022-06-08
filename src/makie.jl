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

function Makie.stairs(h::Hist1D; baseline = 0.0, kwargs...)
    Makie.stairs(binedges(h), vcat(baseline, bincounts(h)); kwargs...)
end
function Makie.stairs!(h::Hist1D; baseline = 0.0, kwargs...)
    Makie.stairs!(binedges(h), vcat(baseline, bincounts(h)); kwargs...)
end

# TODO find the correct way of doing this
# MakieCore.plottype(::Hist1D) = Makie.BarPlot
# function MakieCore.convert_arguments(P::Type{<:MakieCore.AbstractPlot}, h::Hist1D)
#     (Makie.Point2f.(bincenters(h), bincounts(h)), )
# end

# MakieCore.plottype(::Hist2D) = Heatmap
