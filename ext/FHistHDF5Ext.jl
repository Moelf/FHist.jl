module FHistHDF5Ext

using FHist
using HDF5, StatsBase
using Base.Threads: SpinLock

import FHist: h5readhist, h5dumphist


"""
    function h5dumphist(filename::AbstractString, path::AbstractString, h::Hist1D)

Writes a [`Hist1D`](@ref) instance to an HDF5 file. The histogram can be retrieved
using [`h5readhist`](@ref).

# Examples
h = Hist1D(rand(1000), -3:0.3:3)
h5dumphist("foo.h5", "some/path/to/myhist", h)
"""
function h5dumphist(filename::AbstractString, path::AbstractString, h::Union{Hist1D, Hist2D, Hist3D})
    h5open(filename, "cw") do f
        g = create_group(f, path)
        write(g, "weights", h.hist.weights)
        write(g, "sumw2", h.sumw2)
        for (dim, edges) in enumerate(h.hist.edges)
            write(g, "edges_$dim", collect(edges))
        end
        attributes(g)["isdensity"] = h.hist.isdensity
        attributes(g)["closed"] = string(h.hist.closed)
        attributes(g)["overflow"] = string(h.overflow)
        attributes(g)["nentries"] = h.nentries.x
        attributes(g)["_producer"] = "FHist.jl"
    end
    nothing
end


"""
    function h5readhist(filename::AbstractString, path::AbstractString, H::Type{Hist1D})

Reads a histogram from an HDF5 file which has been dumped using [`h5dumphist`](@ref).

# Examples
h = h5readhist("foo.h5", "some/path/to/myhist", Hist1D)
"""
function h5readhist(filename::AbstractString, path::AbstractString, H::Type{<: Union{Hist1D, Hist2D, Hist3D}})
    h5open(filename, "r") do f
        weights = _read_dset(f["$path/weights"], H)
        dims = parse(Int, match(r"\d+", string(H)).match)
        edges = tuple([f["$path/edges_$dim"][:] for dim in 1:dims]...)
        isdensity = read_attribute(f[path], "isdensity")
        closed = Symbol(read_attribute(f[path], "closed"))

        hist = Histogram(edges, weights, closed, isdensity)

        sumw2 = _read_dset(f["$path/sumw2"], H)
        overflow = parse(Bool, read_attribute(f[path], "overflow"))
        nentries = Base.RefValue{Int}(read_attribute(f[path], "nentries"))
        H(hist, sumw2, SpinLock(), overflow, nentries)
    end
end

_read_dset(d::HDF5.Dataset, ::Type{Hist1D}) = d[:]
_read_dset(d::HDF5.Dataset, ::Type{Hist2D}) = d[:, :]
_read_dset(d::HDF5.Dataset, ::Type{Hist3D}) = d[:, :, :]


end