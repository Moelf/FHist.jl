module FHist

export Hist1D, binedges, bincounts, bincenters, binerrors, nbins, integral
export sample, lookup, cumulative, normalize, restrict, rebin
export atomic_push!

export Hist2D, project, profile, transpose

export Hist3D

using StatsBase, RecipesBase, UnicodePlots, Statistics
import LinearAlgebra: normalize, normalize!
using Base.Threads: SpinLock

const _default_overflow = false

struct Hist1D{T<:Real,E} <: AbstractHistogram{T,1,E}
    hist::Histogram{T,1,E}
    sumw2::Array{Float64, 1}
    hlock::SpinLock
    overflow::Bool
    function Hist1D(h::Histogram{T,1,E}, sw2 = copy(h.weights); overflow=_default_overflow) where {T,E}
        return new{T,E}(h, sw2, SpinLock(), overflow)
    end
end

struct Hist2D{T<:Real,E} <: AbstractHistogram{T,2,E}
    hist::Histogram{T,2,E}
    sumw2::Array{Float64, 2}
    hlock::SpinLock
    overflow::Bool
    function Hist2D(h::Histogram{T,2,E}, sw2 = copy(h.weights); overflow=_default_overflow) where {T,E}
        return new{T,E}(h, sw2, SpinLock(), overflow)
    end
end

struct Hist3D{T<:Real,E} <: AbstractHistogram{T,3,E}
    hist::Histogram{T,3,E}
    sumw2::Array{Float64, 3}
    hlock::SpinLock
    overflow::Bool
    function Hist3D(h::Histogram{T,3,E}, sw2 = copy(h.weights); overflow=_default_overflow) where {T,E}
        return new{T,E}(h, sw2, SpinLock(), overflow)
    end
end

include("./utils.jl")
include("./hist1d.jl")
include("./hist2d.jl")
include("./hist3d.jl")
include("./arithmatics.jl")
include("./plot.jl")
end
