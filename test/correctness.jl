using Test, FHist, StatsBase, Statistics, Random

@testset "BinEdges" begin
    b = FHist.BinEdges(0:0.1:1)
    @test FHist.isuniform(b)
    @test length(b) == 11
    @test b == collect(0:0.1:1)
    @test first(b) == 0.0 && last(b) == 1.0
    @test diff(b) == diff(collect(0:0.1:1))
    @test repr(b) == repr(0.0:0.1:1.0)
    @test occursin("Uniform", repr("text/plain", b))
    @test FHist.BinEdges(b) === b
    @test copy(b) == b && copy(b) !== b
    @test hash(b) == hash(FHist.BinEdges(collect(0:0.1:1)))

    # uniform edges given as a vector still get O(1) lookup, exact values are kept
    v = [0.1 * i for i in 0:10]  # 0.30000000000000004 etc.
    bv = FHist.BinEdges(v)
    @test FHist.isuniform(bv)
    @test bv == v
    @test repr(bv) == repr(v)
    @test occursin("Uniform", repr("text/plain", bv))

    # mildly non-uniform edges may still qualify for the (exact) O(1) lookup, the display
    # label reflects the actual spacing though
    bm = FHist.BinEdges([0, 0.5, 0.6, 1])
    @test occursin("Non-uniform", repr("text/plain", bm))
    @test all(x -> searchsortedlast(bm, x) == searchsortedlast([0, 0.5, 0.6, 1], x), -0.1:0.01:1.1)
    bn = FHist.BinEdges([0, 0.01, 0.02, 0.03, 10])
    @test !FHist.isuniform(bn)
    @test occursin("Non-uniform", repr("text/plain", bn))

    # Int vectors and ranges
    @test FHist.BinEdges([0, 1, 2]) == [0.0, 1.0, 2.0]
    @test FHist.BinEdges(0:5) == 0.0:1.0:5.0
    @test FHist.isuniform(FHist.BinEdges(0:5))

    # validation
    @test_throws ArgumentError FHist.BinEdges(10:-1:0)
    @test_throws ArgumentError FHist.BinEdges(range(0, 0, length=5))
    @test_throws ArgumentError FHist.BinEdges(1:1)
    @test_throws ArgumentError FHist.BinEdges([1.0])
    @test_throws ArgumentError FHist.BinEdges(Float64[])
    @test_throws ArgumentError FHist.BinEdges([0, 1, 1])
    @test_throws ArgumentError FHist.BinEdges([0, 2, 1])
    @test_throws ArgumentError FHist.BinEdges([0, 1, NaN])
    @test_throws ArgumentError FHist.BinEdges([0, 1, Inf])
    @test_throws ArgumentError FHist.BinEdges([0, 1e300])
    @test_throws ArgumentError Hist1D(; binedges=10:-1:0)

    # the edges are copied, mutating the input doesn't affect the histogram
    v = [0.0, 1.0, 2.0]
    h = Hist1D(; binedges=v)
    v[2] = 1.5
    @test binedges(h) == [0.0, 1.0, 2.0]

    # exactness of the lookup for values exactly on the edges, uniform vs reference
    for edges in (-3:0.1:3, 0:0.1:1, collect(-3:0.1:3), [0.1 * i for i in -30:30], -6:2:6, 0:5, range(-1, 1, length=1001))
        b = FHist.BinEdges(edges)
        ref = collect(b)
        for x in ref
            @test searchsortedlast(b, x) == searchsortedlast(ref, x)
            @test searchsortedlast(b, prevfloat(x)) == searchsortedlast(ref, prevfloat(x))
            @test searchsortedlast(b, nextfloat(x)) == searchsortedlast(ref, nextfloat(x))
        end
        for x in (NaN, Inf, -Inf, 1e300, -1e300, first(ref) - 1, last(ref) + 1)
            @test searchsortedlast(b, x) == searchsortedlast(ref, x)
        end
    end
end

@testset "Exact edge values land in the right bin" begin
    # values sitting exactly on an edge belong to the bin on their right, like StatsBase
    for edges in (-3:0.1:3, collect(-3:0.1:3))
        data = collect(-3:0.1:3)[1:end-1]  # every lower edge exactly once
        h = Hist1D(data; binedges=edges)
        @test bincounts(h) == ones(60)
        @test convert(Histogram, h) == fit(Histogram, data, -3:0.1:3)
        h2 = Hist1D(; binedges=edges)
        push!.(h2, data)
        @test h2 == h
    end
    h = Hist1D([1, 2, 3, 4, 5]; binedges=1:6)
    @test bincounts(h) == ones(5)
    @test bincounts(Hist1D([0.3, 0.7]; binedges=0:0.1:1)) == bincounts(Hist1D([0.3, 0.7]; binedges=collect(0:0.1:1)))
    @test bincounts(Hist1D([0.3, 0.7]; binedges=0:0.1:1)) == fit(Histogram, [0.3, 0.7], 0:0.1:1).weights
end

@testset "NaN / Inf handling is uniform across binnings" begin
    for edges in (0:0.5:1, [0, 0.5, 1.0])
        for (ov, expected) in ((false, [0.0, 0.0]), (true, [1.0, 3.0]))  # -Inf -> first bin, NaN/Inf/1e300 -> last bin
            h = Hist1D(; binedges=edges, overflow=ov)
            for v in (NaN, Inf, -Inf, 1e300)
                push!(h, v)
            end
            @test bincounts(h) == expected
            @test nentries(h) == (ov ? 4 : 0)
        end
        h = Hist1D([NaN, 0.25, Inf]; binedges=edges)
        @test bincounts(h) == [1.0, 0.0]
        @test nentries(h) == 1
        h2 = Hist2D(([NaN, 0.25, 0.75], [0.25, 0.25, Inf]); binedges=(edges, edges), overflow=true)
        @test bincounts(h2) == [1 0; 1 1]  # NaN x -> last x bin, Inf y -> last y bin
        h3 = Hist3D(([0.25], [0.25], [NaN]); binedges=(edges, edges, edges))
        @test integral(h3) == 0 && nentries(h3) == 0
    end
end

@testset "nentries is consistent across dimensions" begin
    # out-of-range entries are not counted (overflow=false), in every dimension
    h1 = Hist1D([0.5, 5.0]; binedges=0:1)
    h2 = Hist2D(([0.5, 5.0], [0.5, 0.5]); binedges=(0:1, 0:1))
    h3 = Hist3D(([0.5, 5.0], [0.5, 0.5], [0.5, 0.5]); binedges=(0:1, 0:1, 0:1))
    @test nentries(h1) == nentries(h2) == nentries(h3) == 1
    @test integral(h1) == integral(h2) == integral(h3) == 1
    # ...and pushed individually
    h2b = Hist2D(; binedges=(0:1, 0:1))
    push!(h2b, 0.5, 0.5)
    push!(h2b, 5.0, 0.5)
    @test h2b == h2
    h3b = Hist3D(; binedges=(0:1, 0:1, 0:1))
    push!(h3b, 0.5, 0.5, 0.5)
    push!(h3b, 5.0, 0.5, 0.5)
    @test h3b == h3
    # with overflow everything is counted
    @test nentries(Hist2D(([0.5, 5.0], [0.5, 0.5]); binedges=(0:1, 0:1), overflow=true)) == 2
    @test nentries(Hist3D(([0.5, 5.0], [0.5, 0.5], [0.5, 0.5]); binedges=(0:1, 0:1, 0:1), overflow=true)) == 2
    # nentries and integral agree for unweighted data in all cases
    Random.seed!(5)
    x, y, z = randn(1000), randn(1000), randn(1000)
    for ov in (false, true)
        for h in (Hist1D(x; binedges=-1:0.5:1, overflow=ov),
                  Hist2D((x, y); binedges=(-1:0.5:1, [-1, 0, 1.0]), overflow=ov),
                  Hist3D((x, y, z); binedges=(-1:0.5:1, -1:0.5:1, [-1, 0, 1.0]), overflow=ov))
            @test nentries(h) == integral(h)
            @test sumw2(h) == bincounts(h)
        end
    end
end

@testset "Constructor matches StatsBase for every binning type" begin
    Random.seed!(6)
    a = randn(10^4)
    b = randn(10^4)
    c = randn(10^4)
    w = rand(10^4)
    for edges in (-3:0.5:3, collect(-3:0.5:3), [-3, -1, -0.5, 0, 2, 3.0])
        h = Hist1D(a; binedges=edges)
        @test bincounts(h) == fit(Histogram, a, collect(edges)).weights
        hw = Hist1D(a; binedges=edges, weights=w)
        @test bincounts(hw) ≈ fit(Histogram, a, weights(w), collect(edges)).weights
        # sumw2 for weighted data is the sum of squared weights
        hw2 = Hist1D(; binedges=edges)
        for (x, ww) in zip(a, w)
            push!(hw2, x, ww)
        end
        @test sumw2(hw) ≈ sumw2(hw2)
        @test bincounts(hw) ≈ bincounts(hw2)
        @test nentries(hw) == nentries(hw2)

        h2 = Hist2D((a, b); binedges=(edges, edges))
        @test bincounts(h2) == fit(Histogram, (a, b), (collect(edges), collect(edges))).weights
        h2w = Hist2D((a, b); binedges=(edges, edges), weights=w)
        @test bincounts(h2w) ≈ fit(Histogram, (a, b), weights(w), (collect(edges), collect(edges))).weights
        h3 = Hist3D((a, b, c); binedges=(edges, edges, edges))
        @test bincounts(h3) == fit(Histogram, (a, b, c), (collect(edges), collect(edges), collect(edges))).weights
        h3w = Hist3D((a, b, c); binedges=(edges, edges, edges), weights=w)
        @test bincounts(h3w) ≈ fit(Histogram, (a, b, c), weights(w), (collect(edges), collect(edges), collect(edges))).weights
    end
    # Int count type and Float32 data
    @test bincounts(Hist1D(Float32.(a); binedges=-3:0.5:3, counttype=Int)) == fit(Histogram, Float32.(a), -3:0.5:3).weights
end

@testset "Auto binning" begin
    Random.seed!(7)
    x, y, z = randn(1000), randn(1000), randn(1000)
    # 3D auto binning used to produce a single bin per axis
    h3 = Hist3D((x, y, z))
    @test all(nbins(h3) .> 5)
    @test nbins(h3) == (nbins(Hist1D(x)), nbins(Hist1D(y)), nbins(Hist1D(z)))
    @test nbins(Hist3D((x[1:2], y[1:2], z[1:2]))) isa NTuple{3,Int}
    # nbins as a single integer or a tuple
    @test nbins(Hist2D((x, y); nbins=5)) == nbins(Hist2D((x, y); nbins=(5, 5)))
    @test nbins(Hist3D((x, y, z); nbins=4)) == nbins(Hist3D((x, y, z); nbins=(4, 4, 4)))
    @test nbins(Hist1D(x; nbins=7)) == nbins(Hist1D(x; nbins=(7,)))
    @test_throws ArgumentError Hist2D((x, y); nbins=(5, 5, 5))
    # each axis uses its own range
    h2 = Hist2D((x, 10 .* y .+ 100))
    @test first(binedges(h2)[2]) >= 60 && last(binedges(h2)[2]) <= 140
end

cumulative_or_self(h::Hist1D) = cumulative(h)
cumulative_or_self(h) = h
@testset "Arithmetic keeps uniform bin edges" begin
    h = Hist1D(randn(100); binedges=-3:0.5:3)
    h2 = Hist2D((randn(100), randn(100)); binedges=(-3:1:3, -3:1:3))
    h3 = Hist3D((randn(100), randn(100), randn(100)); binedges=(-3:1:3, -3:1:3, -3:1:3))
    for hh in (h, h2, h3)
        for r in (hh + hh, hh - hh, hh * 2, 2 * hh, hh / hh, normalize(hh), merge(hh, hh), cumulative_or_self(hh))
            @test all(FHist.isuniform, r.binedges)
            @test all(b -> b.isrange, r.binedges)
            @test binedges(r) == binedges(hh)
        end
    end
    @test first(split(repr(h + h), '\n')) == first(split(repr(h), '\n')) == "edges: -3.0:0.5:3.0"
    # scaling a histogram with negative bins is allowed (it used to throw)
    hneg = h - 2 * h
    @test any(<(0), bincounts(hneg))
    @test bincounts(hneg * 2) == 2 .* bincounts(hneg)
    @test sumw2(hneg * 2) == 4 .* sumw2(hneg)
    @test integral(normalize(hneg; width=false)) ≈ 1
    # division error propagation
    ha = Hist1D([0.5, 1.5, 1.5, 2.5]; binedges=0:3)
    hb = Hist1D([0.5, 1.5, 2.5, 2.5]; binedges=0:3)
    @test sumw2(ha / hb) ≈ [2.0, 6.0, 0.375]
end
@testset "empty!, hash and ==" begin
    h = Hist1D(randn(100); binedges=-3:3)
    @test empty!(h) === h
    @test nentries(h) == 0
    @test all(iszero, bincounts(h)) && all(iszero, sumw2(h))
    @test h == Hist1D(; binedges=-3:3)
    h2 = Hist2D((randn(10), randn(10)); binedges=(-3:3, -3:3))
    @test empty!(h2) === h2 && nentries(h2) == 0
    h3 = Hist3D((randn(10), randn(10), randn(10)); binedges=(-3:3, -3:3, -3:3))
    @test empty!(h3) === h3 && nentries(h3) == 0

    ha = Hist1D([0.5]; binedges=0:1)
    hb = Hist1D([0.5]; binedges=[0, 1.0])
    @test ha == hb
    @test hash(ha) == hash(hb)
    @test length(Set([ha, hb])) == 1
    @test hash(ha) != hash(Hist1D([0.5]; binedges=0:1, overflow=true))
    @test ha != Hist1D([0.5]; binedges=0:1, overflow=true)
end

@testset "append! for all dimensions" begin
    Random.seed!(8)
    x, y, z, w = randn(100), randn(100), randn(100), rand(100)
    h2 = Hist2D(; binedges=(-3:3, -3:3))
    @test append!(h2, x, y) === h2
    @test h2 == Hist2D((x, y); binedges=(-3:3, -3:3))
    h2w = Hist2D(; binedges=(-3:3, -3:3))
    append!(h2w, x, y, w)
    @test h2w == Hist2D((x, y); binedges=(-3:3, -3:3), weights=w)
    h3 = Hist3D(; binedges=(-3:3, -3:3, -3:3))
    @test append!(h3, x, y, z) === h3
    @test h3 == Hist3D((x, y, z); binedges=(-3:3, -3:3, -3:3))
    h3w = Hist3D(; binedges=(-3:3, -3:3, -3:3))
    append!(h3w, x, y, z, w)
    @test h3w == Hist3D((x, y, z); binedges=(-3:3, -3:3, -3:3), weights=w)
    @test_throws DimensionMismatch append!(h2, x, y[1:10])
    @test_throws DimensionMismatch append!(h3, x, y, z, w[1:10])
    h1 = Hist1D(; binedges=-3:3)
    @test append!(h1, x) === h1
    # a failing push! (InexactError for Int counts) must not leave the lock held
    hi = Hist1D(; counttype=Int, binedges=-3:3)
    @test_throws InexactError append!(hi, [0.5], [0.5])
    @test append!(hi, [0.5]) === hi
end

@testset "integral / normalize with width in 2D/3D" begin
    h2 = Hist2D(; binedges=([0, 1, 3.0], [0, 2.0]))
    push!(h2, 0.5, 1.0, 2.0)  # bin area 2
    push!(h2, 2.0, 1.0, 3.0)  # bin area 4
    @test integral(h2) == 5
    @test integral(h2; width=true) == 2 * 2 + 3 * 4
    n2 = normalize(h2)
    @test integral(n2) ≈ 1
    @test all(FHist.isuniform.(n2.binedges) .== FHist.isuniform.(h2.binedges))
    n2w = normalize(h2; width=true)
    @test integral(n2w; width=true) ≈ 1
    @test bincounts(n2w) ≈ reshape([2 / 5 / 2, 3 / 5 / 4], 2, 1)
    @test sumw2(n2w) ≈ sumw2(n2) ./ reshape([2, 4], 2, 1) .^ 2

    h3 = Hist3D(; binedges=([0, 1, 3.0], [0, 2.0], [0, 0.5]))
    push!(h3, 0.5, 1.0, 0.25, 2.0)  # volume 1
    push!(h3, 2.0, 1.0, 0.25, 3.0)  # volume 2
    @test integral(h3) == 5
    @test integral(h3; width=true) == 2 * 1 + 3 * 2
    @test integral(normalize(h3)) ≈ 1
    @test integral(normalize(h3; width=true); width=true) ≈ 1
end

@testset "rebin / restrict for all dimensions" begin
    Random.seed!(9)
    x, y, z = rand(1000), rand(1000), rand(1000)
    r = 0:0.1:1
    v = [0.0, 0.1, 0.3, 0.6, 1.0]

    # curried integer rebin works for every dimension
    h1 = Hist1D(x; binedges=r)
    h2 = Hist2D((x, y); binedges=(r, r))
    h3 = Hist3D((x, y, z); binedges=(r, r, r))
    @test (h2 |> rebin(5)) == rebin(h2, 5, 5)
    @test (h3 |> rebin(5)) == rebin(h3, 5, 5, 5)
    @test (h3 |> rebin(5, 2, 10)) == rebin(h3, 5, 2, 10)

    # 3D integer rebin
    r3 = rebin(h3, 2, 5, 10)
    @test nbins(r3) == (5, 2, 1)
    @test integral(r3) == integral(h3)
    @test sum(sumw2(r3)) == sum(sumw2(h3))
    @test nentries(r3) == nentries(h3)
    @test binedges(r3) == (0:0.2:1, 0:0.5:1, [0, 1.0])
    @test all(b -> b.isrange, r3.binedges)
    @test r3 == Hist3D((x, y, z); binedges=(0:0.2:1, 0:0.5:1, 0:1))
    @test_throws ErrorException rebin(h3, 3)
    @test rebin(h3, 1) == h3

    # edge based rebin, 1D/2D/3D
    @test rebin(h1, [0.0, 0.2, 0.5, 1.0]) == Hist1D(x; binedges=[0.0, 0.2, 0.5, 1.0])
    hsub = rebin(h1, [0.2, 0.5])
    @test bincounts(hsub) == bincounts(Hist1D(x[0.2 .<= x .< 0.5]; binedges=[0.2, 0.5]))
    @test binedges(hsub) == [0.2, 0.5]
    @test nentries(hsub) == nentries(h1)  # nentries is carried over as is
    @test rebin(h2, [0.0, 0.5, 1.0], [0.0, 0.2, 1.0]) == Hist2D((x, y); binedges=([0.0, 0.5, 1.0], [0.0, 0.2, 1.0]))
    @test rebin(h3, [0.0, 0.5, 1.0], [0.0, 0.2, 1.0], [0.0, 1.0]) == Hist3D((x, y, z); binedges=([0.0, 0.5, 1.0], [0.0, 0.2, 1.0], [0.0, 1.0]))
    @test_throws ArgumentError rebin(h2, [0.0, 0.55, 1.0], [0.0, 1.0])
    # the overflow flag is dropped when the new edges do not span the full range
    ho = Hist1D(x; binedges=r, overflow=true)
    @test rebin(ho, [0.0, 0.5, 1.0]).overflow
    @test !rebin(ho, [0.2, 0.5]).overflow
    @test rebin(Hist2D((x, y); binedges=(r, r), overflow=true), 2).overflow

    # rebinning keeps the exact original edge values (no re-derivation from first/last)
    hv = Hist1D(x; binedges=v)
    @test binedges(rebin(hv, 2)) == [0.0, 0.3, 1.0]
    hv2 = Hist2D((x, y); binedges=(v, v))
    @test binedges(rebin(hv2, 2)) == ([0.0, 0.3, 1.0], [0.0, 0.3, 1.0])

    # 3D restrict
    rs = restrict(h3, 0.2, 0.55, -Inf, Inf, 0.5, Inf)
    @test nbins(rs) == (4, 10, 5)
    @test binedges(rs)[1] == 0.2:0.1:0.6
    @test binedges(rs)[3] == 0.5:0.1:1.0
    @test integral(rs) == count((0.2 .<= x .< 0.6) .& (0.5 .<= z))
    @test rs == (h3 |> restrict(0.2, 0.55, -Inf, Inf, 0.5, Inf))
    @test restrict(h3) == h3
    @test_throws AssertionError restrict(h3, 5, 6)
    # 1D/2D restrict keep ranges as ranges
    @test binedges(restrict(h1, 0.2, 0.55)) == 0.2:0.1:0.6
    @test restrict(h1, 0.2, 0.55).binedges[1].isrange
    @test all(b -> b.isrange, restrict(h2, 0.2, 0.55).binedges)
    @test binedges(restrict(hv, 0.2, 0.7)) == [0.1, 0.3, 0.6]  # bin centers 0.2 and 0.45
end

@testset "Statistics for Hist3D" begin
    Random.seed!(10)
    x, y, z = randn(10^5), randn(10^5) .+ 1, 2 .* randn(10^5)
    h3 = Hist3D((x, y, z); binedges=(-5:0.1:5, -5:0.1:5, -10:0.1:10))
    mx, my, mz = mean(h3)
    @test mx ≈ 0 atol = 0.05
    @test my ≈ 1 atol = 0.05
    @test mz ≈ 0 atol = 0.1
    sx, sy, sz = std(h3)
    @test sx ≈ 1 atol = 0.05
    @test sz ≈ 2 atol = 0.1
    @test median(h3) isa NTuple{3,Float64}
    @test mean(h3)[1] == mean(project(project(h3, :z), :x))
    @test mean(h3)[3] == mean(project(project(h3, :x), :y))
end

@testset "Display" begin
    h = Hist1D(randn(100); binedges=-3:1.0:3)
    @test occursin("edges: -3.0:1.0:3.0", repr(h))
    hv = Hist1D(randn(100); binedges=[-3, 0, 3.0])
    @test occursin("edges: [-3.0, 0.0, 3.0]", repr(hv))
    h2 = Hist2D((randn(10), randn(10)); binedges=(-3:3, -3:3))
    @test occursin("edges: (-3.0:1.0:3.0, -3.0:1.0:3.0)", repr(h2))
    @test occursin("-3.0:1.0:3.0", repr("text/html", h2))
end

@testset "Misc regressions" begin
    @test FHist._factor(12) == Set([2, 3])
    @test FHist._factor(97) == Set([97])
    @test FHist._factor(1) == Set{Int}()
    @test FHist._factor(360) == Set([2, 3, 5])
    @test FHist._factor(49) == Set([7])
    # Histogram conversion keeps ranges
    h = Hist1D(randn(100); binedges=-3:3)
    @test convert(Histogram, h).edges[1] isa AbstractRange
    @test convert(Histogram, Hist1D(randn(100); binedges=[-3, 0, 3.0])).edges[1] == [-3, 0, 3.0]
    # 0-allocation push! for all dimensions
    function _pushloop2(h, a, b)
        for i in eachindex(a)
            push!(h, a[i], b[i])
        end
    end
    function _pushloop3(h, a, b, c)
        for i in eachindex(a)
            push!(h, a[i], b[i], c[i])
        end
    end
    a, b, c = randn(100), randn(100), randn(100)
    h2 = Hist2D(; binedges=(-3:3, [-3, 0, 3.0]))
    _pushloop2(h2, a, b)
    @test (@allocated _pushloop2(h2, a, b)) == 0
    h3 = Hist3D(; binedges=(-3:3, -3:3, [-3, 0, 3.0]))
    _pushloop3(h3, a, b, c)
    @test (@allocated _pushloop3(h3, a, b, c)) == 0
    # non-uniform edges also allocate nothing
    hn = Hist1D(; binedges=[-3, -1, 0, 0.5, 3.0])
    function _pushloop(h, a)
        for v in a
            push!(h, v)
        end
    end
    _pushloop(hn, a)
    @test (@allocated _pushloop(hn, a)) == 0
    # lookup on the edge / consistent with push!
    hl = Hist1D([0.3]; binedges=0:0.1:1)
    @test lookup(hl, 0.3) == 1
    @test lookup(hl, prevfloat(0.3)) == 0
    # sample keyword
    @test length(sample(hl; n=3)) == 3
end
