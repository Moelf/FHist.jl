using .Makie

@recipe(StackedHist) do scene
    Theme(
        line_color = :red
    )
end

function Makie.plot!(input::StackedHist{<:Tuple{AbstractVector{<:Hist1D}}}; errors = true, kwargs...)
    hs = input[1][]
    Nhist = length(hs)
    _e = binedges(first(hs))
    all(==(_e), binedges.(hs)) || throw("binedges must match in stacked histogram")

    centers = bincenters(first(hs))
    Nbin = length(centers)
    xs = repeat(centers; outer=Nhist)
    ys = mapreduce(bincounts, vcat, hs)
    grp = repeat(eachindex(hs); inner=Nbin)
    mes = mapreduce(h->(bincounts(h), binerrors(h)), ((c1, e1), (c2, e2)) -> (c1 .± e1) + (c2 .± e2), hs)
    totals = Measurements.value.(mes)
    errs = Measurements.uncertainty.(mes)
    
    
    Makie.barplot!(input, xs, ys,
        stack = grp,
        color = grp; kwargs...
    )
    
    if errors
        errorbars!(input, centers, totals,errs/2,
            color = range(0, 1, length = length(xs)),
            whiskerwidth = 10)
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
