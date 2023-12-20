module FHistHDF5Ext

using FHist
using HDF5, StatsBase
using Base.Threads: SpinLock

import FHist: h5readhist, h5writehist

const CURRENT_H5HIST_VERSION = v"1.0"
const SUPPORTED_H5HIST_VERSIONS = [v"1.0"]


"""
    function h5writehist(filename::AbstractString, path::AbstractString, h::Hist1D)

Writes a [`Hist1D`](@ref) instance to an HDF5 file. The histogram can be retrieved
using [`h5readhist`](@ref).

# Examples
h = Hist1D(rand(1000), -3:0.3:3)
h5writehist("foo.h5", "some/path/to/myhist", h)
"""
function h5writehist(filename::AbstractString, path::AbstractString, h::Union{Hist1D, Hist2D, Hist3D})
    h5open(filename, "cw") do f
        g = create_group(f, path)
        write(g, "weights", h.hist.weights)
        write(g, "sumw2", h.sumw2)
        for (dim, edges) in enumerate(h.hist.edges)
            write(g, "edges_$dim", collect(edges))
        end
        attributes(g)["overflow"] = h.overflow
        attributes(g)["nentries"] = h.nentries[]
        attributes(g)["_producer"] = "FHist.jl"
        attributes(g)["_h5hist_version"] = string(CURRENT_H5HIST_VERSION)
    end
    nothing
end


"""
    function h5readhist(filename::AbstractString, path::AbstractString [, H::Type{<: Union{Hist1D, Hist2D, Hist3D}}])
    function h5readhist(f::HDF5.File, path::AbstractString [, H::Type{<: Union{Hist1D, Hist2D, Hist3D}}])

Reads a histogram from an HDF5 file which has been written using
[`h5writehist`](@ref). The type parameter is optional.

# Examples
h = h5readhist("foo.h5", "some/path/to/myhist")
"""
function h5readhist(filename::AbstractString, path::AbstractString, H::Type{<: Union{Hist1D, Hist2D, Hist3D}})
    h5open(filename, "r") do f
        h5readhist(f, path, H)
    end
end
function h5readhist(f::HDF5.File, path::AbstractString)
    dims = length(filter(g->startswith(g, "edges_"), keys(f[path])))
    dims == 1 && return h5readhist(f, path, Hist1D)
    dims == 2 && return h5readhist(f, path, Hist2D)
    dims == 3 && return h5readhist(f, path, Hist3D)
    error("Histograms with $dims dimensions are not supported yet.")
end
function h5readhist(filename::AbstractString, path::AbstractString)
    h5open(filename, "r") do f
        h5readhist(f, path)
    end
end
function h5readhist(f::HDF5.File, path::AbstractString, H::Type{<: Union{Hist1D, Hist2D, Hist3D}})
    version = VersionNumber(read_attribute(f[path], "_h5hist_version"))
    version >= v"2" && error("h5hist $(version) is not supported")
    version > sort(SUPPORTED_H5HIST_VERSIONS)[end] && @warn """
        h5hist $(version) is higher than the currently supperted one, some features might be missing.
        Supported versions: $(join(SUPPORTED_H5HIST_VERSIONS, ", "))
    """

    weights = _read_dset(f["$path/weights"], H)
    dims = parse(Int, match(r"\d+", string(H)).match)
    edges = tuple([f["$path/edges_$dim"][:] for dim in 1:dims]...)

    hist = Histogram(edges, weights)

    sumw2 = _read_dset(f["$path/sumw2"], H)
    overflow = read_attribute(f[path], "overflow")
    nentries = Base.RefValue{Int}(read_attribute(f[path], "nentries"))
    H(hist, sumw2, SpinLock(), overflow, nentries)
end

_read_dset(d::HDF5.Dataset, ::Type{Hist1D}) = d[:]
_read_dset(d::HDF5.Dataset, ::Type{Hist2D}) = d[:, :]
_read_dset(d::HDF5.Dataset, ::Type{Hist3D}) = d[:, :, :]


end
