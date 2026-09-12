module FHist

export Hist1D, binedges, bincounts, bincenters, binerrors, nbins, integral, nentries, significance
export sample, lookup, cumulative, normalize, restrict, rebin, bayes_rebin_edges, sumw2
export atomic_push!

export Hist2D, project, profile, transpose

export Hist3D, collabtext!, statbox!

using StatsBase, Statistics
export Weights
import LinearAlgebra: normalize, normalize!
using Base.Threads: SpinLock

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
                length(es) == $N || throw(DimensionMismatch("Binedges must be a tuple of $($N) vectors, got $(length(es))"))
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

            `nbins` can be a single integer (used for every axis) or a tuple with one integer per axis.

            Values that fall outside of `binedges` are discarded (and not counted in `nentries`)
            unless `overflow=true`, in which case they are clamped into the first/last bin along
            each axis. `NaN` is treated like `+Inf`.
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
            _fast_bincounts!(h, ary, weights)
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
        bincenters(h::$H) = _from_tuple(map(b -> StatsBase.midpoints(b.edges), h.binedges))
        @doc """
            nentries(h::$($H))
        Get the number of entries that were filled (`push!`ed) into the histogram. Values that
        were discarded because they fell outside of the bin edges (with `overflow=false`) are not
        counted.
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

        @doc raw"""
            effective_entries(h) -> scalar

        Get the number of effective entries for the entire histogram:

        ```math
        n_\text{eff} = \frac{(\sum \text{Weights} )^2}{(\sum \text{Weight}^2 )}
        ```

        This is also equivalent to `integral(hist)^2 / sum(sumw2(hist))`, this is the same as `TH1::GetEffectiveEntries()`
        """
        effective_entries(h::$H) = abs2(integral(h)) / sum(sumw2(h))

        function Base.:(==)(h1::$H, h2::$H)
            bincounts(h1) == bincounts(h2) &&
                h1.binedges == h2.binedges &&
                nentries(h1) == nentries(h2) &&
                sumw2(h1) == sumw2(h2) &&
                h1.overflow == h2.overflow
        end

        Base.hash(h::$H, x::UInt) = hash(h.overflow, hash(nentries(h), hash(sumw2(h), hash(h.binedges, hash(bincounts(h), x)))))

        @doc """
            empty!(h)

        Reset the histogram in place: bin counts, `sumw2` and `nentries` are all set to zero. The
        bin edges and `overflow` setting are kept. Returns `h`.
        """
        function Base.empty!(h::$H)
            bincounts(h) .= false
            sumw2(h) .= false
            h.nentries[] = 0
            return h
        end

        Base.broadcastable(h::$H) = Ref(h)
    end
end

# The bin index along one axis (1-based, `0` when the value is to be discarded), given the
# `BinEdges` of that axis, the number of bins `L` and the `overflow` policy.
@inline function _binindex(b::BinEdges, L::Int, overflow::Bool, x::Real)
    if b.isuniform & !b.twosided
        # fused fast path: one accept test and two loadless guesses, see `_find_bias`
        xf = Float64(x)
        if (xf >= b.rfirst) & (xf < b.rlast)  # false for NaN
            return _biased_lookup(b, xf)
        else
            return overflow ? (xf < b.rfirst ? 1 : L) : 0  # NaN -> L
        end
    end
    i = searchsortedlast(b, x)
    if overflow
        return clamp(i, 1, L)
    else
        return unsigned(i - 1) < unsigned(L) ? i : 0
    end
end

# `_fast_bincounts!` fills a *freshly constructed, empty* histogram with data; this is what
# the "fit like" constructors call. Unweighted fills only touch `bincounts` in the loop and set
# `sumw2 = bincounts` afterwards (valid only because `h` starts empty).
function _fast_bincounts!(h::Hist1D, A, weights)
    xs = A[1]
    b = h.binedges[1]
    L = nbins(h)
    overflow = h.overflow
    counts = bincounts(h)
    n = 0
    if isnothing(weights)
        for x in xs
            i = _binindex(b, L, overflow, x)
            i == 0 && continue
            n += 1
            @inbounds counts[i] += one(eltype(counts))
        end
        sumw2(h) .= counts
    else
        s2 = sumw2(h)
        for (x, w) in zip(xs, weights)
            i = _binindex(b, L, overflow, x)
            i == 0 && continue
            n += 1
            @inbounds counts[i] += w
            @inbounds s2[i] += w^2
        end
    end
    h.nentries[] += n
    return h
end

function _fast_bincounts!(h::Hist2D, A, weights)
    xs, ys = A
    bx, by = h.binedges
    Lx, Ly = nbins(h)
    overflow = h.overflow
    counts = bincounts(h)
    n = 0
    if isnothing(weights)
        for (x, y) in zip(xs, ys)
            ix = _binindex(bx, Lx, overflow, x)
            iy = _binindex(by, Ly, overflow, y)
            (ix == 0 || iy == 0) && continue
            n += 1
            @inbounds counts[ix, iy] += one(eltype(counts))
        end
        sumw2(h) .= counts
    else
        s2 = sumw2(h)
        for (x, y, w) in zip(xs, ys, weights)
            ix = _binindex(bx, Lx, overflow, x)
            iy = _binindex(by, Ly, overflow, y)
            (ix == 0 || iy == 0) && continue
            n += 1
            @inbounds counts[ix, iy] += w
            @inbounds s2[ix, iy] += w^2
        end
    end
    h.nentries[] += n
    return h
end

function _fast_bincounts!(h::Hist3D, A, weights)
    xs, ys, zs = A
    bx, by, bz = h.binedges
    Lx, Ly, Lz = nbins(h)
    overflow = h.overflow
    counts = bincounts(h)
    n = 0
    if isnothing(weights)
        for (x, y, z) in zip(xs, ys, zs)
            ix = _binindex(bx, Lx, overflow, x)
            iy = _binindex(by, Ly, overflow, y)
            iz = _binindex(bz, Lz, overflow, z)
            (ix == 0 || iy == 0 || iz == 0) && continue
            n += 1
            @inbounds counts[ix, iy, iz] += one(eltype(counts))
        end
        sumw2(h) .= counts
    else
        s2 = sumw2(h)
        for (x, y, z, w) in zip(xs, ys, zs, weights)
            ix = _binindex(bx, Lx, overflow, x)
            iy = _binindex(by, Ly, overflow, y)
            iz = _binindex(bz, Lz, overflow, z)
            (ix == 0 || iy == 0 || iz == 0) && continue
            n += 1
            @inbounds counts[ix, iy, iz] += w
            @inbounds s2[ix, iy, iz] += w^2
        end
    end
    h.nentries[] += n
    return h
end

include("./utils.jl")
include("./hist1d.jl")
include("./hist2d.jl")
include("./hist3d.jl")
include("./displays.jl")
include("./arithmatics.jl")

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


end
