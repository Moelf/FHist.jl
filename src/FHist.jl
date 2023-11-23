module FHist

export Hist1D, binedges, bincounts, bincenters, binerrors, nbins, integral, nentries, significance
export sample, lookup, cumulative, normalize, restrict, rebin, bayes_rebin_edges
export atomic_push!

export Hist2D, project, profile, transpose

export Hist3D, collabtext!, statbox!

using StatsBase, Statistics, Measurements
export Weights
import LinearAlgebra: normalize, normalize!
using Base.Threads: SpinLock

if !isdefined(Base, :get_extension)
    using Requires
end
#using Requires

using BayesHistogram
export BayesHistogram

const _default_overflow = false

for (H, N) in ((:Hist1D, 1), (:Hist2D, 2), (:Hist3D, 3))
    @eval struct $H{T<:Real,E} <: AbstractHistogram{T,$N,E}
        hist::Histogram{T,$N,E}
        sumw2::Array{Float64, $N}
        hlock::SpinLock
        overflow::Bool
        nentries::Base.RefValue{Int}
        function $H(h::Histogram{T,$N,E}, sw2::AbstractArray{<:Real, $N} = copy(h.weights), nentries=sum(h.weights); overflow=_default_overflow) where {T,E}
            isinteger(nentries) || @warn "Weights was probably used but StatsBase.Histogram doesn't record # of entries"
            return new{T,E}(h, sw2, SpinLock(), overflow, Ref(round(Int, nentries)))
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

using MakieCore
include("./MakieThemes.jl")
export ATLASTHEME, stackedhist, stackedhist!, ratiohist, ratiohist!

function stackedhist end
function stackedhist! end

function ratiohist end
function ratiohist! end

function statbox! end
function collabtext! end

function __init__()

    @static if !isdefined(Base, :get_extension)
        @require Plots="91a5bcdd-55d7-5caf-9e0b-520d859cae80" include("../ext/FHistPlotsExt.jl")
        @require Makie="ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a" include("../ext/FHistMakieExt.jl")
        @require CairoMakie="13f3f980-e62b-5c42-98c6-ff1f3baf88f0" include("../ext/FHistCairoMakieExt.jl")
    end
    
end
end
