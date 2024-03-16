"""
    BinEdges <: AbstractVector{Float64}

    This type implements a vector-like data structure to be used for histogram bin edges, it can handle both uniform and non-uniform binnings in a single type to reduce the amount of parametric types. It would switch to O(1) `searchsortedlast` if the binning is uniform.

!!! note
    Due to the usage of Float64, bin edges shouldn't contain element with absolute value larger than 9007199254740992, which is the `maxintfloat(Float64)`.
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
        elseif edges isa AbstractVector
            if !issorted(edges)
                throw(ArgumentError("BinEdges must be sorted"))
            end
            if !allunique(edges)
                throw(ArgumentError("BinEdges must be unique"))
            end
            if any(x -> abs(x) > maxintfloat(Float64), edges)
                throw(ArgumentError("BinEdges cannot contain element with absolute value larger than $(maxintfloat(Float64))"))
            end
            new(false, edges, 0.0:-1.0, Inf, -Inf)
        end
    end
end

Base.@constprop :aggressive isuniform(b::BinEdges) = b.isuniform
Base.size(b::BinEdges) = if isuniform(b) size(b.uniform_edges) else size(b.nonuniform_edges) end
Base.getindex(b::BinEdges, i) = if isuniform(b) b.uniform_edges[i] else b.nonuniform_edges[i] end
# Base.collect(b::BinEdges) = if isuniform(b) collect(b.uniform_edges) else copy(b.nonuniform_edges) end

Base.convert(::Type{BinEdges}, edges::AbstractRange) = BinEdges(edges)
Base.convert(::Type{BinEdges}, edges::AbstractVector) = BinEdges(edges)

function Base.searchsortedlast(r::BinEdges, x::Real)
    if isuniform(r)
        return floor(Int, (x - r.rfirst) * r.inv_step) + 1
    else
        return searchsortedlast(r.nonuniform_edges, x)
    end
end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, b::BinEdges) 
    if isuniform(b)
        print(io, "Uniform FHist.BinEdges: ")
        show(io, mime, b.uniform_edges)
    else
        print(io, "Non-uniform FHist.BinEdges: ")
        show(io, b.nonuniform_edges)
    end
end
