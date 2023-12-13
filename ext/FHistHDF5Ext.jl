module FHistHDF5Ext

using FHist
using HDF5, StatsBase, Base.Threads: SpinLock

import FHist: h5readhist, h5dumphist


"""
    function h5dumphist(filename::AbstractString, path::AbstractString, h::Hist1D)

Writes a [`Hist1D`](@ref) instance to an HDF5 file. The histogram can be retrieved
using [`h5readhist`](@ref).

# Examples
h = Hist1D(rand(1000), -3:0.3:3)
h5dumphist("foo.h5", "some/path/to/myhist", h)
"""
function h5dumphist(filename::AbstractString, path::AbstractString, h::Hist1D)
    h5open(filename, "cw") do f

        g = create_group(f, path)
        write(g, "weights", h.hist.weights)
        write(g, "sumw2", h.sumw2)
        write(g, "edges", collect(h.hist.edges[1]))
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
function h5readhist(filename::AbstractString, path::AbstractString, H::Type{Hist1D})
    h5open(filename, "r") do f
        weights = f["$path/weights"][:]
        edges = f["$path/edges"][:]
        isdensity = read_attribute(f[path], "isdensity")
        closed = Symbol(read_attribute(f[path], "closed"))

        hist = Histogram(edges, weights, closed, isdensity)

        sumw2 = f["$path/sumw2"][:]
        overflow = parse(Bool, read_attribute(f[path], "overflow"))
        nentries = Base.RefValue{Int}(read_attribute(f[path], "nentries"))
        H(hist, sumw2, SpinLock(), overflow, nentries)
    end
end


end
