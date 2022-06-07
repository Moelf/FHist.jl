using .Makie

@recipe(StackedHist) do scene
    Theme(
        line_color = :red
    )
end

function Makie.plot!(input::StackedHist{<:Tuple{AbstractVector{<:Hist1D}}}; kwargs...)
    hs = input[1][]
    Nhist = length(hs)
    _e = binedges(first(hs))
    all(==(_e), binedges.(hs)) || throw("binedges must match in stacked histogram")
    
    centers = bincenters(first(hs))
    Nbin = length(centers)
    xs = repeat(centers; outer=Nhist)
    ys = mapreduce(bincounts, vcat, hs)
    grp = repeat(eachindex(hs); inner=Nbin)
    
    Makie.barplot!(input, xs, ys,
        stack = grp,
        color = grp; kwargs...
    )
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
