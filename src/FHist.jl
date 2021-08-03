module FHist

export Hist1D, sample

using StatsBase, RecipesBase
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

include("./hist1d.jl")
function binerrors(f::Function, h::Hist1D)
    f.(h.hist.weights)
end

include("./arithmatics.jl")
include("./plot.jl")
include("./utils.jl")
end
