module FHist

export Hist1D, binedges, bincounts, bincenters, binerrors, nbins, integral, nentries, significance
export sample, lookup, cumulative, normalize, restrict, rebin, bayes_rebin_edges, sumw2
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

using BayesHistogram
export BayesHistogram

const _default_overflow = false
_to_tuple(x::Tuple) = x
_to_tuple(x) = tuple(x)

_from_tuple(x::Tuple{Any}) = only(x)
_from_tuple(x) = x


# constructor-like API
for (H, N) in ((:Hist1D, 1), (:Hist2D, 2), (:Hist3D, 3))
    @eval struct $H{T<:Real,E<:NTuple{$N, AbstractVector}} <: AbstractHistogram{T,$N,E}
        binedges::E
        bincounts::Array{T,$N}
        sumw2::Array{Float64,$N}
        nentries::Base.RefValue{Int}
        overflow::Bool
        hlock::SpinLock
        function $H(;
            counttype::Type{T}=Float64,
            binedges,
            bincounts=zeros(counttype, length.(_to_tuple(binedges)) .- 1),
            sumw2=zero(bincounts),
            nentries=0,
            overflow=false) where {T}

            es = _to_tuple(binedges)
            all(length.(es) .- 1 .== size(bincounts) .== size(sumw2)) ||
                throw(DimensionMismatch("Binedges must be tuple of each axes, and each dimension has one more than the corresponding
                dimension of `bincounts`"))

                return new{T,typeof(es)}(es, bincounts, sumw2, Ref(round(Int, nentries)), overflow, SpinLock())
        end
    end

    @eval begin
        bincounts(h::$H) = h.bincounts
        binedges(h::$H) = _from_tuple(h.binedges)
        bincenters(h::$H) = _from_tuple(StatsBase.midpoints.(h.binedges))
        nentries(h::$H) = h.nentries[]
        sumw2(h::$H) = h.sumw2

        function Base.:(==)(h1::$H, h2::$H)
            bincounts(h1) == bincounts(h2) &&
            binedges(h1) == binedges(h2) &&
            nentries(h1) == nentries(h2) &&
            sumw2(h1) == sumw2(h2) &&
            h1.overflow == h2.overflow
        end

        function Base.empty!(h1::$H)
            bincounts(h1) .= false
            sumw2(h1) .= false
        end
    end
end

function auto_bins(ary, ::Val{1}; nbins=nothing)
    xs = only(ary)
    E = eltype(xs)
    F = E <: Number ? float(E) : Float64
    nbins = isnothing(nbins) ? _sturges(xs) : nbins
    lo, hi = minimum(xs), maximum(xs)
    (StatsBase.histrange(F(lo), F(hi), nbins),)
end

function auto_bins(ary, ::Val{2}; nbins=nothing)
    xs, ys = ary
    E = eltype(xs)
    xnbins, ynbins = isnothing(nbins) ? _sturges.(xs) : nbins
    F = E <: Number ? float(E) : Float64
    lo, hi = minimum(xs), maximum(xs)
    loy, hiy = minimum(ys), maximum(ys)
    (StatsBase.histrange(F(lo), F(hi), xnbins),
        StatsBase.histrange(F(loy), F(hiy), ynbins),)
end

function auto_bins(ary, ::Val{3}; nbins=nothing)
    xs, ys, zs = ary
    E = eltype(xs)
    xnbins, ynbins, znbins = isnothing(nbins) ? _sturges.(xs) : nbins
    F = E <: Number ? float(E) : Float64
    lo, hi = minimum(xs), maximum(xs)
    loy, hiy = minimum(ys), maximum(ys)
    loz, hiz = minimum(zs), maximum(zs)
    (StatsBase.histrange(F(lo), F(hi), xnbins),
        StatsBase.histrange(F(loy), F(hiy), ynbins),
        StatsBase.histrange(F(loz), F(hiz), znbins),)
end

# fit-like API
for (H, N) in ((:Hist1D, 1), (:Hist2D, 2), (:Hist3D, 3))
    @eval function $H(_ary::E;
        weights=nothing,
        nbins = nothing,
        binedges = nothing,
        counttype::Type{T}=Float64,
        overflow=false) where {T, E<:NTuple{$N, Any}}

        ary = _to_tuple(_ary)
        binedges = if !isnothing(binedges)
            binedges
        else
            auto_bins(ary, Val($N); nbins)
        end

        h = $H(; binedges=_to_tuple(binedges), overflow=overflow)
        if isnothing(weights)
            for t in zip(ary...)
                push!(h, t...)
            end
        else
            for (t, w) in zip(zip(ary...), weights)
                push!(h, t..., w)
            end
        end
        return h
    end
end

include("./utils.jl")
include("./hist1d.jl")
include("./hist2d.jl")
include("./hist3d.jl")
include("./displays.jl")
include("./arithmatics.jl")

using MakieCore
include("./MakieThemes.jl")
export ATLASTHEME, stackedhist, stackedhist!, ratiohist, ratiohist!

function stackedhist end
function stackedhist! end

function ratiohist end
function ratiohist! end

function statbox! end
function collabtext! end

export h5writehist, h5readhist
function h5writehist end
function h5readhist end

function __init__()

    @static if !isdefined(Base, :get_extension)
        @require Plots="91a5bcdd-55d7-5caf-9e0b-520d859cae80" include("../ext/FHistPlotsExt.jl")
        @require Makie="ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a" include("../ext/FHistMakieExt.jl")
        @require HDF5="f67ccb44-e63f-5c2f-98bd-6dc0ccc4ba2f" include("../ext/FHistHDF5Ext.jl")
    end

end
end
