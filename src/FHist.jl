module FHist

export Hist1D, binedges, bincounts, bincenters, binerrors, nbins, integral, nentries, significance
export sample, lookup, cumulative, normalize, restrict, rebin
export atomic_push!

export Hist2D, project, profile, transpose

export Hist3D, collabtext!, statbox!, ATLASTHEME

using StatsBase, UnicodePlots, Statistics, Measurements
import LinearAlgebra: normalize, normalize!
using Base.Threads: SpinLock

using Requires

const _default_overflow = false

for (H, N) in ((:Hist1D, 1), (:Hist2D, 2), (:Hist3D, 3))
    @eval struct $H{T<:Real,E} <: AbstractHistogram{T,$N,E}
        hist::Histogram{T,$N,E}
        sumw2::Array{Float64, $N}
        hlock::SpinLock
        overflow::Bool
        nentries::Ref{Int}
        function $H(h::Histogram{T,$N,E}, sw2 = copy(h.weights), nentries=0; overflow=_default_overflow) where {T,E}
            return new{T,E}(h, sw2, SpinLock(), overflow, Ref(nentries))
        end
    end
end

include("./utils.jl")
include("./hist1d.jl")
include("./hist2d.jl")
include("./hist3d.jl")
include("./displays.jl")
include("./arithmatics.jl")
nentries(h::Union{Hist1D, Hist2D, Hist3D}) = h.nentries[]

function __init__()
    @require Plots="91a5bcdd-55d7-5caf-9e0b-520d859cae80" include("./plots.jl")
    @require Makie="ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a" include("./makie.jl")
    @require CairoMakie="13f3f980-e62b-5c42-98c6-ff1f3baf88f0" include("./cairomakie.jl")
end
end
