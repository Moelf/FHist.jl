"""
    BinEdges <: AbstractVector{Float64}

This type implements a vector-like data structure to be used for histogram bin edges, it can handle both uniform and non-uniform binning in a single type thus reduce the amount of parametric typing in packages like FHist.jl. It would switch to fast implementation of functions like `searchsortedfirst` and `searchsortedlast` based if the binning is uniform or not.
"""
struct BinEdges <: AbstractVector{Float64}
    isuniform::Bool
    nonuniform_edges::Vector{Float64}
    uniform_edges::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}
    inv_step::Float64
    rfirst::Float64
    function BinEdges(edges)
        if isempty(edges)
            throw(ArgumentError("BinEdges cannot be empty"))
        end
        if edges isa AbstractRange
            new(true, Float64[], edges, inv(step(edges)), first(edges))
        else
            if !issorted(edges)
                throw(ArgumentError("BinEdges must be sorted"))
            end
            if !allunique(edges)
                throw(ArgumentError("BinEdges must be unique"))
            end
            new(false, edges, 0.0:-1.0, Inf, -Inf)
        end
    end
end

Base.@constprop :aggressive isuniform(b::BinEdges) = b.isuniform
Base.size(b::BinEdges) = if isuniform(b) size(b.uniform_edges) else size(b.nonuniform_edges) end
Base.getindex(b::BinEdges, i) = if isuniform(b) b.uniform_edges[i] else b.nonuniform_edges[i] end

Base.convert(::Type{BinEdges}, edges::AbstractRange) = BinEdges(edges)
Base.convert(::Type{BinEdges}, edges::AbstractVector) = BinEdges(edges)

function Base.searchsortedlast(r::BinEdges, x::Real)
    if isuniform(r)
        return floor(Int, (x - first(r)) * r.inv_step) + 1
    else
        return searchsortedlast(r.nonuniform_edges, x)
    end
end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, b::BinEdges) 
    if isuniform(b)
        println(io, "Uniform binning: $(first(b)) : $(last(b)) with step $(step(b))")
    else
        println(io, "Non-uniform BinEdges:")
        show(io, mime, b.nonuniform_edges)
    end
end
