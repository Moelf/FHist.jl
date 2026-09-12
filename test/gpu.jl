# The GPU extension is exercised on the KernelAbstractions `CPU()` backend, which runs the very
# same kernels as the GPU backends.
using Test, FHist, KernelAbstractions, GPUArraysCore, Random

@testset "GPU extension (CPU backend)" begin
    ext = Base.get_extension(FHist, :FHistGPUExt)
    @test ext !== nothing

    Random.seed!(11)
    xs = randn(10^5)
    ws = rand(10^5)
    # values exactly on edges, out of range, non-finite
    append!(xs, [-3.0, -2.7, -2.2, 0.0, 0.3, 0.7, 2.9, 3.0, 10.0, -10.0, NaN, Inf, -Inf])
    append!(ws, ones(13))
    xs32 = Float32.(xs)
    ws32 = Float32.(ws)

    edges_list = (
        -3:0.1:3,                       # range
        collect(-3:0.1:3),              # uniform vector
        [-3.0, -2.5, -2.2, -1.0, 0.3, 0.35, 0.7, 1.0, 3.0],  # non-uniform
        range(-3, 3, length=20001),     # too many bins for shared memory -> global atomics
    )
    for edges in edges_list, overflow in (false, true), (data, w) in ((xs, ws), (xs32, ws32))
        big = length(edges) > 10_000
        for naive in (big ? (true,) : (false, true))
            href = Hist1D(data; binedges=edges, overflow)
            c = gpu_bincounts(data; binedges=edges, overflow, naive, counttype=Float64)
            @test c isa Vector{Float64}
            @test c == bincounts(href)

            hrefw = Hist1D(data; binedges=edges, overflow, weights=w)
            cw = gpu_bincounts(data; binedges=edges, overflow, naive, counttype=Float64, weights=w)
            @test cw ≈ bincounts(hrefw) rtol = 1e-6

            # full Hist1D (bincounts, sumw2 and nentries), unweighted and weighted
            h = Hist1D(; binedges=edges, overflow)
            ext._gpu_hist1d!(h, data, nothing; naive)
            @test h == href
            hw = Hist1D(; binedges=edges, overflow)
            ext._gpu_hist1d!(hw, data, w; naive)
            @test nentries(hw) == nentries(hrefw)
            @test bincounts(hw) ≈ bincounts(hrefw) rtol = 1e-6
            @test sumw2(hw) ≈ sumw2(hrefw) rtol = 1e-6
        end
    end

    # in-place accumulation
    c = zeros(Float32, 60)
    gpu_bincounts!(c, xs; binedges=-3:0.1:3)
    gpu_bincounts!(c, xs; binedges=-3:0.1:3)
    @test c == 2 .* bincounts(Hist1D(xs; binedges=-3:0.1:3))
    @test_throws DimensionMismatch gpu_bincounts!(zeros(Float32, 10), xs; binedges=-3:0.1:3)

    # count types
    @test gpu_bincounts(xs; binedges=-3:0.1:3) isa Vector{Float32}
    @test gpu_bincounts(xs; binedges=-3:0.1:3, counttype=Int32) == bincounts(Hist1D(xs; binedges=-3:0.1:3))
    # empty data
    @test gpu_bincounts(Float64[]; binedges=0:1) == [0.0]
    # bad inputs
    @test_throws DimensionMismatch gpu_bincounts(xs; binedges=-3:0.1:3, weights=ws[1:10])
    @test_throws ArgumentError gpu_bincounts(xs; binedges=range(-3, 3, length=20001), naive=false)
    @test_throws ArgumentError gpu_bincounts(xs; binedges=[0.0, 1.0 + 1e-12, 1.0 + 2e-12], edgetype=Float32)  # both round up to nextfloat(1f0)
end
