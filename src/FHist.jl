module FHist

export Hist1D, binedges, bincounts, bincenters, nbins, integral
export sample, lookup, cumulative, normalize, rebin
export unsafe_push!

export Hist2D

using StatsBase, RecipesBase, UnicodePlots, Statistics
import LinearAlgebra: normalize, normalize!
using Base.Threads: SpinLock

@inline function pearson_err(n::Real)
    s = sqrt(n+0.25)
    s+0.5, s-0.5
end
@inline function sqrt_err(n::Real)
    s = sqrt(n)
    s,s
end

@inline function _is_uniform_bins(A::AbstractVector{T}) where T<:Real
    diffs = diff(A)
    diff1 = first(diffs)
    all(isapprox.(diff1, diffs; atol = 1e-9)) #hack
end
function _is_uniform_bins(A::AbstractRange{T}) where T<:Real
    true
end

struct Hist1D{T<:Real,E} <: AbstractHistogram{T,1,E}
    hist::Histogram{T,1,E}
    sumw2::Array{Float64, 1}
    hlock::SpinLock
    function Hist1D(h::Histogram{T,1,E}, sw2 = copy(h.weights)) where {T,E}
        return new{T,E}(h, sw2, SpinLock())
    end
end

struct Hist2D{T<:Real,E} <: AbstractHistogram{T,2,E}
    hist::Histogram{T,2,E}
    sumw2::Array{Float64, 2}
    hlock::SpinLock
    function Hist2D(h::Histogram{T,2,E}, sw2 = copy(h.weights)) where {T,E}
        return new{T,E}(h, sw2, SpinLock())
    end
end

const Hist = Union{Hist1D, Hist2D}

include("./hist1d.jl")
include("./hist2d.jl")
function binerrors(f::Function, h::Hist1D)
    f.(h.hist.weights)
end

include("./arithmatics.jl")
include("./plot.jl")
include("./utils.jl")
end
