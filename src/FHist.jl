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

for (H, N) in ((:Hist1D, 1), (:Hist2D, 2), (:Hist3D, 3))
    @eval struct $H{T<:Real,E} <: AbstractHistogram{T,$N,E}
        hist::Histogram{T,$N,E}
        sumw2::Array{Float64, $N}
        hlock::SpinLock
        overflow::Bool
        function $H(h::Histogram{T,$N,E}, sw2 = copy(h.weights); overflow=_default_overflow) where {T,E}
            return new{T,E}(h, sw2, SpinLock(), overflow)
        end
    end
end

include("./utils.jl")
include("./hist1d.jl")
include("./hist2d.jl")
include("./hist3d.jl")
include("./displays.jl")
include("./arithmatics.jl")
include("./plot.jl")
end
