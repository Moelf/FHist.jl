import MakieCore

MakieCore.@recipe(StackedHist, hs) do scene
    Theme(
        line_color = :red
    )
end

function Makie.plot!(myplot::StackedHist{<:Tuple{AbstractVector{<:Hist1D}}})
    Nhist = length(hs)
    _e = binedges(first(hs))
    @assert all(==(_e), binedges.(hs))
    
    centers = bincenters(first(hs))
    xs = repeat(centers; outer=Nhist)
    ys = mapreduce(bincounts, vcat, hs)
    grp = repeat(eachindex(hs); outer=length(centers))
    
    Makie.barplot!(myplot, xs, ys,
        stack = grp,
        color = grp; kwargs...
    )
    myplot
end
# function Makie.stairs(h::Hist1D; baseline = 0.0, kwargs...)
#     Makie.stairs(binedges(h), vcat(baseline, bincounts(h)); kwargs...)
# end
function Makie.stairs!(h::Hist1D; baseline = 0.0, kwargs...)
    Makie.stairs!(binedges(h), vcat(baseline, bincounts(h)); kwargs...)
end

# TODO find the correct way of doing this
# MakieCore.plottype(::Hist1D) = Makie.BarPlot
# function MakieCore.convert_arguments(P::Type{<:MakieCore.AbstractPlot}, h::Hist1D)
#     (Makie.Point2f.(bincenters(h), bincounts(h)), )
# end

# MakieCore.plottype(::Hist2D) = Heatmap
