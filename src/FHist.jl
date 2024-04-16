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

_to_tuple(x::Tuple) = x
_to_tuple(x) = tuple(x)

_from_tuple(x::Tuple{Any}) = only(x)
_from_tuple(x) = x

include("./polybinedges.jl")

for (H, N) in ((:Hist1D, 1), (:Hist2D, 2), (:Hist3D, 3))

    @eval begin
        struct $H{T<:Real} <: AbstractHistogram{T,$N,NTuple{$N, BinEdges}}
            binedges::NTuple{$N, BinEdges}
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

                return new{T}(es, bincounts, sumw2, Ref(round(Int, nentries)), overflow, SpinLock())
            end
        end

        @doc """
            # To make an __empty__ histogram
            
            use the all-keyword-arguments constructor:
            ```julia
            $($H)(;
                counttype=Float64,
                binedges::E
                bincounts = zeros(counttype, length.(_to_tuple(binedges)) .- 1),
                sumw2 = zero(bincounts),
                nentries = 0
                overflow::Bool = false
            ) where {E<:NTuple{$($N),Any}}
            ```

            !!! note
                Everything other than `binedges` are optional (infered from `binedges`).


            # To make an histogram given data (and weights etc.)
            
            use the a positional argument for data and keyword-arguments for the rest:
            ```julia
            $($H)(array::E;
                counttype=Float64,
                binedges=nothing,
                weights=nothing,
                nbins=nothing,
                overflow=false)
            ) where {E<:NTuple{$($N),Any}}
            ```

            !!! note
                Everything other than data (`array`) is optional (infered from data).
            """
        function $H(ary::E;
            counttype::Type{T}=Float64,
            binedges=nothing,
            weights=nothing,
            nbins=nothing,
            overflow=false) where {T,E<:NTuple{$N,Any}}

            length(ary) == $N || throw(DimensionMismatch("Data must be a tuple of $N vectors"))
            isnothing(weights) || length(ary[1]) == length(weights) || throw(DimensionMismatch("Data and weights must have the same length"))

            binedges = if !isnothing(binedges)
                binedges
            else
                auto_bins(ary, Val($N); nbins)
            end

            bs = _to_tuple(binedges)
            h = $H(; counttype, binedges=bs, overflow=overflow)
            _fast_bincounts!(h, ary, binedges, weights)
            return h
        end

        Base.lock(h::$H) = lock(h.hlock)
        Base.unlock(h::$H) = unlock(h.hlock)
        @doc """
                bincounts(h::$($H))
            Get the bin counts (weights) of the histogram.
        """
        bincounts(h::$H) = h.bincounts
        @doc """
            binedges(h)

        Get the bin edges of the histogram

        !!! note
            For 1D histogram, it returns just a vector. For others, it returns a tuple of vectors. If you need a tuple of vectors, use `h.binedges` at your own risk.
        """
        binedges(h::$H) = _from_tuple(h.binedges)

        @doc """
            bincenters(h::$($H))
        Get the bin centers of the histogram

        !!! note
            For 1D histogram, it returns just a vector. For others, it returns a tuple of vectors.
        """
        bincenters(h::$H) = _from_tuple(StatsBase.midpoints.(h.binedges))
        @doc """
                nentries(h::$($H))
            Get the number of times a histogram is filled (`push!`ed)
        """
        nentries(h::$H) = h.nentries[]
        @doc """
                sumw2(h)
            Get the sum of weights squared of the histogram, it has the same shape as `bincounts(h)`.
        """
        sumw2(h::$H) = h.sumw2

        @doc """
                binerrors(f=sqrt, h)
            Get the error (uncertainty) of each bin. By default, calls `sqrt` on `sumw2(h)` bin by bin as an approximation.
        """
        binerrors(f::T, h::$H) where T<:Function = f.(sumw2(h))
        binerrors(h::$H) = binerrors(sqrt, h)

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

# fall back one-shot implementation
function _fast_bincounts!(h::Hist1D, A, binedges, weights)
    xs = A[1]
    if isnothing(weights)
        for x in xs
            push!(h, x)
        end
    else
        for (x, w) in zip(xs, weights)
            push!(h, x, w)
        end
    end
end
function _fast_bincounts!(h::Hist2D, A, binedges, weights)
    xs, ys = A
    if isnothing(weights)
        for (x, y) in zip(xs, ys)
            push!(h, x,y)
        end
    else
        for (x, y, w) in zip(xs, ys, weights)
            push!(h, x, y, w)
        end
    end
end
function _fast_bincounts!(h::Hist3D, A, binedges, weights)
    xs, ys, zs = A
    if isnothing(weights)
        for (x, y, z) in zip(xs, ys, zs)
            push!(h, x, y, z)
        end
    else
        for (x, y, z, w) in zip(xs, ys, zs, weights)
            push!(h, x, y, z, w)
        end
    end
end

include("./utils.jl")
export chi2ndf_ww
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
