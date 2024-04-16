using FHist, StatsBase, Statistics, HDF5
using Test

@testset "Fast route" begin
    a = rand(10^5)
    r1 = 0:0.1:1
    r2 = collect(r1)
    for r in (r1, r2)
        h1 = Hist1D(a; binedges=r)
        @test nentries(h1) == 10^5
        sth1 = fit(Histogram, a, r)
        @test all(bincounts(h1) .≈ sth1.weights)
    end

    @test_throws DimensionMismatch Hist1D(rand(10); weights=rand(9))
    @test_throws DimensionMismatch Hist2D((rand(10), rand(10)); weights=rand(9))
end

@testset "Slow route" begin
    a = rand(10^5)
    r1 = [0, 0.3, 0.4, 0.8, 1]
    h1 = Hist1D(a; binedges=r1)
    sth1 = fit(Histogram, a, r1)
    @test nentries(h1) == 10^5
    @test all(bincounts(h1) .≈ sth1.weights)
end

@testset "push! route" begin
    a = rand(10^5)
    r1 = [0, 0.3, 0.4, 0.8, 1]
    h1 = Hist1D(a; binedges=r1)
    h2 = Hist1D(; counttype=Int, binedges=r1)
    @test nentries(h1) == 10^5
    sth1 = fit(Histogram, a, r1)
    for ele in a
        push!(h2, ele)
    end

    @test all(bincounts(h1) .≈ sth1.weights)
    @test all(bincounts(h1) .≈ bincounts(h2))
end

@testset "1D append!" begin
    a = rand(10^5)
    r1 = [0, 0.3, 0.4, 0.8, 1]
    h1 = Hist1D(; binedges=r1)
    h2 = Hist1D(; binedges=r1)
    h3 = Hist1D(; binedges=r1)
    for ele in a
        push!(h1, ele)
    end

    append!(h2, a)
    @test all(bincounts(h1) .≈ bincounts(h2))

    append!(h3, a, zero(a) .+ 0.1)
    @test all(bincounts(h1) .* 0.1 .≈ bincounts(h3))
end

@testset "Weighted" begin
    a = rand(10^5)
    wgts1 = ones(10^5) |> weights
    wgts = rand(10^5) |> weights
    r1 = 0:0.1:1
    r2 = collect(r1)

    h1 = Hist1D(a)
    h2 = Hist1D(a; weights=wgts1)
    @test nentries(h1) == 10^5
    @test nentries(h2) == 10^5
    @test h1 == h2
    for r in (r1, r2)
        h1 = Hist1D(a; weights=wgts, binedges=r)
        sth1 = fit(Histogram, a, weights(wgts), r)
        @test all(bincounts(h1) .≈ sth1.weights)
    end
end

@testset "Hist2D" begin
    x = rand(10)
    y = rand(10)
    h = Hist2D((x, y))
    @test integral(h) == 10
    @test nentries(h) == 10

    rx, ry = 0:0.1:1, 0:0.2:1
    wgts = weights(2 * ones(length(x)))
    h = Hist2D((x, y); weights=wgts, binedges=(rx, ry))
    @test nentries(h) == 10
    @test integral(h) == sum(wgts)
    @test nbins(h) == (length(rx) - 1, length(ry) - 1)

    @test bincenters(Hist2D((x, y); binedges=(0:1, 0:1))) == ([0.5], [0.5])
    @test bincenters(Hist2D((x, y); weights=wgts, binedges=(0:1, 0:1))) == ([0.5], [0.5])
    @test nbins(Hist2D((x, y); binedges=([0, 0.5, 1], [0, 0.5, 1]))) == (2, 2)
    @test nbins(Hist2D((x, y); binedges=([0, 0.3, 1], [0, 0.3, 1]))) == (2, 2)
    @test nbins(Hist2D((x, y); weights=wgts, binedges=([0, 0.5, 1], [0, 0.5, 1]))) == (2, 2)
    @test nbins(Hist2D((x, y); weights=wgts, binedges=([0, 0.3, 1], [0, 0.3, 1]))) == (2, 2)

    @test integral(Hist2D((x, y); weights=wgts, nbins=(5, 5))) == sum(wgts)
    @test integral(Hist2D((x, y); nbins=(5, 5))) == length(x)

    h1 = Hist2D((x, y); weights=wgts, binedges=(rx, ry))
    h2 = Hist2D(; counttype=Float64, binedges=(rx, ry))
    h3 = Hist2D(; counttype=Float64, binedges=(rx, ry))
    push!.(h2, x, y, wgts)
    atomic_push!.(h3, x, y, wgts)
    @test h1 == h2
    @test h1 == h3
end

@testset "Hist3D" begin
    x = rand(10)
    y = rand(10)
    z = rand(10)
    h = Hist3D((x, y, z))
    @test integral(h) == 10
    @test nentries(h) == 10

    rx, ry, rz = 0:0.1:1, 0:0.2:1, 0:0.5:1
    wgts = weights(2 * ones(length(x)))
    h = Hist3D((x, y, z); weights=wgts, binedges=(rx, ry, rz))
    @test nentries(h) == 10
    @test integral(h) == sum(wgts)
    @test nbins(h) == (length(rx) - 1, length(ry) - 1, length(rz) - 1)

    @test bincenters(Hist3D((x, y, z); binedges=(0:1, 0:1, 0:1))) == ([0.5], [0.5], [0.5])
    @test bincenters(Hist3D((x, y, z); weights=wgts, binedges=(0:1, 0:1, 0:1))) == ([0.5], [0.5], [0.5])
    @test nbins(Hist3D((x, y, z), binedges=([0, 0.5, 1], [0, 0.5, 1], [0, 0.5, 1]))) == (2, 2, 2)
    @test nbins(Hist3D((x, y, z), binedges=([0, 0.3, 1], [0, 0.3, 1], [0, 0.3, 1]))) == (2, 2, 2)
    @test nbins(Hist3D((x, y, z), weights=wgts, binedges=([0, 0.5, 1], [0, 0.5, 1], [0, 0.5, 1]))) == (2, 2, 2)
    @test nbins(Hist3D((x, y, z), weights=wgts, binedges=([0, 0.3, 1], [0, 0.3, 1], [0, 0.3, 1]))) == (2, 2, 2)

    @test integral(Hist3D((x, y, z); weights=wgts, nbins=(5, 5, 5))) == sum(wgts)
    @test integral(Hist3D((x, y, z); nbins=(5, 5, 5))) == length(x)

    h1 = Hist3D((x, y, z); weights=wgts, binedges=(rx, ry, rz))
    h2 = Hist3D(; counttype=Float64, binedges=(rx, ry, rz))
    h3 = Hist3D(; counttype=Float64, binedges=(rx, ry, rz))
    push!.(h2, x, y, z, wgts)
    atomic_push!.(h3, x, y, z, wgts)
    @test h1 == h2
    @test h1 == h3
end

@testset "Iterable fall back" begin
    # 1D
    a = Iterators.flatten([[1, 2], [], [3]])
    b = collect(a)
    r = 0:0.1:1
    @test Hist1D(a; binedges=r) == Hist1D(b; binedges=r)

    # 2D
    x1 = Iterators.flatten(rand(10))
    y1 = Iterators.flatten(rand(10))
    h1 = Hist2D((x1, y1); binedges=(r, r))
    x2 = collect(x1)
    y2 = collect(y1)
    h2 = Hist2D((x2, y2); binedges=(r, r))
    @test h1 == h2
end

@testset "Special bins" begin
    # integer values and integer binning
    a = floor.(Int, abs.(randn(10^6)))
    h1 = Hist1D(a; binedges=0:5)
    h2 = fit(Histogram, a, 0:5)
    @test convert(Histogram, h1) == h2

    # integer values and integer binning
    a = floor.(Int, abs.(randn(10^6)))
    h1 = Hist1D(a; binedges=-6:2:6)
    h2 = fit(Histogram, a, -6:2:6)
    @test convert(Histogram, h1) == h2
end

@testset "Normalize" begin
    a = rand(10^5)
    wgts1 = 2 .* ones(10^5) |> weights
    h1 = Hist1D(a; weights=wgts1)
    @test integral(h1) ≈ sum(wgts1) atol = 1e-8
    @test integral(normalize(h1; width=false)) ≈ 1 atol = 1e-8
    @test nentries(normalize(h1; width=false)) == length(a)

    h1 = Hist2D((a, a); weights=wgts1, binedges=(0:0.1:1, 0:0.1:1))
    @test integral(h1) ≈ sum(wgts1) atol = 1e-8
    @test integral(normalize(h1)) ≈ 1 atol = 1e-8

    h1 = Hist3D((a, a, a); weights=wgts1, binedges=(0:0.1:1, 0:0.1:1, 0:0.1:1))
    @test integral(h1) ≈ sum(wgts1) atol = 1e-8
    @test integral(normalize(h1)) ≈ 1 atol = 1e-8
end

@testset "Normalize with width" begin
    t = Hist1D(; binedges=[0, 1, 2, 4])
    push!(t, 0.5, 0.5)
    push!(t, 0.5, 0.5)
    push!(t, 2.5, 0.5)
    # 0.5 + 0.5 + 0.5 = 1.5

    @test integral(t) == 1.5
    @test integral(t; width=true) == 2

    nt = normalize(t; width=false)
    @test bincounts(nt) == bincounts(t) ./ 1.5
    @test integral(nt) == 1 # self-consistent requirement

    ntw = normalize(t; width=true)
    @test bincounts(ntw) == bincounts(t) ./ 2
    @test integral(ntw; width=true) == 1 #self-consistent
end


@testset "Cumulative" begin
    a = rand(10^5)
    wgts1 = 2 .* ones(10^5) |> weights
    h1 = Hist1D(a; weights=wgts1)
    @test argmax(bincounts(cumulative(h1; forward=true))) == nbins(h1)
    @test argmax(bincounts(cumulative(h1; forward=false))) == 1
    h1 = Hist1D([1, 2, 2, 3, 3, 3]; binedges=1:4)
    @test bincounts(cumulative(h1, forward=true)) == [1, 3, 6]
    @test bincounts(cumulative(h1, forward=true)) == sumw2(cumulative(h1, forward=true))
    @test bincounts(cumulative(h1, forward=false)) == [6, 5, 3]
    @test bincounts(cumulative(h1, forward=false)) == sumw2(cumulative(h1, forward=false))
end

@testset "Bin errors" begin
    h1 = Hist1D(randn(100); binedges=-3:3)
    @test h1.sumw2 == binerrors(identity, h1)
    @test sqrt.(h1.sumw2) == binerrors(h1)
    h1 = Hist2D((randn(100), randn(100)); binedges=(-3:3, -3:3))
    @test h1.sumw2 == binerrors(identity, h1)
    @test sqrt.(h1.sumw2) == binerrors(h1)
    h1 = Hist3D((randn(100), randn(100), randn(100)); binedges=(-3:3, -3:3, -3:3))
    @test h1.sumw2 == binerrors(identity, h1)
    @test sqrt.(h1.sumw2) == binerrors(h1)

    h2 = Hist1D(; counttype=Float64, binedges=0:1)
    push!(h2, 0.5, 0.5)
    push!(h2, 0.5, 0.5)

    @test h2.sumw2 == [0.5]
    @test binerrors(h2) == [sqrt(0.5)]
end

@testset "Statistics" begin
    r = -3:0.1:3
    N = 10^5
    h1 = Hist1D(randn(N); binedges=r)
    @test mean(h1) ≈ 0 atol = 0.05
    @test std(h1) ≈ 1 atol = 0.1
    @test nbins(h1) == length(r) - 1
    @test median(h1) == quantile(h1, 0.5)
    @test quantile(h1, 0.0) == first(bincenters(h1))
    @test quantile(h1, 1.0) == last(bincenters(h1))

end

@testset "Lookup" begin
    h1 = Hist1D(randn(100))
    @test lookup.(Ref(h1), bincenters(h1)) == bincounts(h1)
    @test ismissing(lookup(h1, last(binedges(h1)) + 0.1))
    @test ismissing(lookup(h1, first(binedges(h1)) - 0.1))

    h1 = Hist2D((randn(100), randn(100)); binedges=(-3:3, -3:3))
    cx, cy = bincenters(h1)
    points = [(x, y) for x in cx, y in cy]
    @test map(p -> lookup(h1, p...), points) == bincounts(h1)
    @test ismissing(lookup(h1, 10, 10))

    h1 = Hist3D((randn(100), randn(100), randn(100)); binedges=(-3:3, -3:3, -3:3))
    cx, cy, cz = bincenters(h1)
    points = [(x, y, z) for x in cx, y in cy, z in cz]
    @test map(p -> lookup(h1, p...), points) == bincounts(h1)
    @test ismissing(lookup(h1, 10, 10, 10))
end

@testset "Sample" begin
    @test mean(FHist.sample(Hist1D(rand(10^5); binedges=0:0.1:1), n=10^5)) ≈ 0.5 atol = 0.1
end

@testset "Empty" begin
    a = rand(10^5)
    r = 0:0.1:1
    h1 = Hist1D(a; binedges=r)
    empty!(h1)
    @test maximum(bincounts(h1)) == 0
    @test maximum(sumw2(h1)) == 0
    @test binedges(h1) == r

    h1 = Hist2D((randn(10), randn(10)); binedges=(-3:3, -3:3))
    empty!(h1)
    @test maximum(bincounts(h1)) == 0
    @test maximum(sumw2(h1)) == 0

    h1 = Hist3D((randn(10), randn(10), randn(10)); binedges=(-3:3, -3:3, -3:3))
    empty!(h1)
    @test maximum(bincounts(h1)) == 0
    @test maximum(sumw2(h1)) == 0
end

@testset "Atomic push" begin
    if get(ENV, "CI", "false") == "true"
        @test Threads.nthreads() > 1
    end
    a = randn(10^5)
    h1 = Hist1D(a; binedges=-3:0.2:3)
    h2 = Hist1D(; counttype=Int, binedges=-3:0.2:3)
    h3 = Hist1D(; counttype=Int, binedges=-3:0.2:3)
    h4 = Hist1D(; counttype=Int, binedges=-3:0.2:3)
    for i in a
        push!(h2, i)
        atomic_push!(h3, i)
    end
    Threads.@threads for i in a
        atomic_push!(h4, i)
    end
    @test h1 == h2
    @test h1 == h3
    @test h1 == h4
end

@testset "Broadcast" begin
    h1 = Hist1D(; counttype=Int, binedges=0:3)

    push!.(h1, [0.5, 1.5])
    @test lookup.(h1, [0.5, 1.5, 2.5]) == [1, 1, 0]

    empty!(h1)
    push!.(h1, [0.5, 1.5], 2)
    @test lookup.(h1, [0.5, 1.5, 2.5]) == [2, 2, 0]
end

@testset "Arithmetic" begin
    @testset "Unweighted regular binning" begin
        h1 = Hist1D([0.5, 1.5, 1.5, 2.5]; binedges=0:3)
        h2 = Hist1D([0.5, 1.5, 2.5, 2.5]; binedges=0:3)

        h = h1 + h2
        @test bincounts(h) == [2.0, 3.0, 3.0]
        @test sumw2(h) == [2.0, 3.0, 3.0]
        @test nentries(h) == nentries(h1) + nentries(h2)

        h = h1 - h2
        @test bincounts(h) == [0.0, 1.0, -1.0]
        @test sumw2(h) == [2.0, 3.0, 3.0]

        h = h1 / h2
        @test bincounts(h) == [1.0, 2.0, 0.5]
        @test sumw2(h) == [2.0, 6.0, 0.375]

        h = h1 * 2
        @test bincounts(h) == [2.0, 4.0, 2.0]
        @test sumw2(h) == [4.0, 8.0, 4.0]

        h = h1 / (h1 + h2 * 2)
        @test bincounts(h) ≈ [0.333333, 0.5, 0.2] atol = 1e-6
        @test sumw2(h) ≈ [0.17284, 0.21875, 0.0544] atol = 1e-6
    end

    @testset "Weighted regular binning" begin
        h1 = Hist1D([0.5, 1.5, 1.5, 2.5]; weights=[3, 3, 2, 1], binedges=0:3)
        h2 = Hist1D([0.5, 1.5, 2.5, 2.5]; weights=[3, 3, 2, 1], binedges=0:3)

        h = h1 / (h1 + h2 * 2)
        @test bincounts(h) ≈ [0.333333, 0.454545, 0.142857] atol = 1e-6
        @test sumw2(h) ≈ [0.17284, 0.191107, 0.029155] atol = 1e-6
    end

    @testset "Weighted irregular binning" begin
        h1 = Hist1D([0.5, 1.5, 1.5, 2.5]; weights=[3, 3, 2, 1], binedges=[0, 1, 2, 4])
        h2 = Hist1D([0.5, 1.5, 2.5, 2.5]; weights=[3, 3, 2, 1], binedges=[0, 1, 2, 4])

        h = h1 / (h1 + h2 * 2)
        @test bincounts(h) ≈ [0.333333, 0.454545, 0.142857] atol = 1e-6
        @test sumw2(h) ≈ [0.17284, 0.191107, 0.029155] atol = 1e-6
    end

    @testset "Hist2D" begin
        h1 = Hist2D(([0.5, 1.5, 1.5, 2.5], [0.5, 0.5, 0.5, 0.5]); binedges=(0:3, 0:1))
        h2 = Hist2D(([0.5, 1.5, 2.5, 2.5], [0.5, 0.5, 0.5, 0.5]); binedges=(0:3, 0:1))
        h = h1 / (h1 + h2 * 2)
        @test vec(bincounts(h)) ≈ [0.333333, 0.5, 0.2] atol = 1e-6
        @test vec(sumw2(h)) ≈ [0.17284, 0.21875, 0.0544] atol = 1e-6
    end

    @testset "Hist3D" begin
        h1 = Hist3D(([0.5, 1.5, 1.5, 2.5], [0.5, 0.5, 0.5, 0.5], [0.5, 0.5, 0.5, 0.5]); binedges=(0:3, 0:1, 0:1))
        h2 = Hist3D(([0.5, 1.5, 2.5, 2.5], [0.5, 0.5, 0.5, 0.5], [0.5, 0.5, 0.5, 0.5]); binedges=(0:3, 0:1, 0:1))
        h = h1 / (h1 + h2 * 2)
        @test vec(bincounts(h)) ≈ [0.333333, 0.5, 0.2] atol = 1e-6
        @test vec(sumw2(h)) ≈ [0.17284, 0.21875, 0.0544] atol = 1e-6
    end

end

@testset "Merging" begin
    h1 = Hist1D(randn(100); binedges=-3:3)
    h2 = Hist1D(randn(100); binedges=-3:3)
    @test merge(h1, h2) == h1 + h2
    @test merge(h1, h2, h1) == h1 + h2 + h1

    h1 = Hist2D((randn(10), randn(10)); binedges=(-3:3, -3:3))
    h2 = Hist2D((randn(10), randn(10)); binedges=(-3:3, -3:3))
    @test merge(h1, h2) == h1 + h2

    h1 = Hist3D((randn(10), randn(10), randn(10)); binedges=(-3:3, -3:3, -3:3))
    h2 = Hist3D((randn(10), randn(10), randn(10)); binedges=(-3:3, -3:3, -3:3))
    @test merge(h1, h2) == h1 + h2
end

@testset "Rebinning" begin
    h1 = Hist1D(rand(10^2); binedges=0:0.1:1)
    @test h1 == rebin(h1, 1)
    @test nentries(h1) == nentries(rebin(h1, 1))
    @test integral(h1) == integral(rebin(h1, 5))
    @test sum(h1.sumw2) == sum(rebin(h1, 5).sumw2)
    @test binedges(rebin(h1, 5)) == [0, 0.5, 1.0]
    @test Set([2, 5]) == FHist.valid_rebin_values(h1)

    h2 = Hist1D(rand(10^2); binedges=[0.0, 0.1, 0.7, 0.9, 1.0])
    @test h2 == rebin(h2, 1)
    @test integral(h2) == integral(rebin(h2, 2))
    @test nentries(h2) == nentries(rebin(h2, 1))
    @test sum(h2.sumw2) == sum(rebin(h2, 2).sumw2)
    @test binedges(rebin(h2, 2)) == [0, 0.7, 1.0]
    @test Set([2]) == FHist.valid_rebin_values(h2)

    @test rebin(h1, 2) == (h1 |> rebin(2))

    h1 = Hist2D((rand(10^2), rand(10^2)); binedges=(0:0.1:1, 0:0.1:1))
    @test h1 == rebin(h1, 1, 1)
    @test integral(h1) == integral(rebin(h1, 5))
    @test sum(h1.sumw2) == sum(rebin(h1, 5).sumw2)
    @test binedges(rebin(h1, 5)) == ([0, 0.5, 1.0], [0, 0.5, 1.0])
    @test Set([2]) == FHist.valid_rebin_values(h2)

    bins = [0.0, 0.1, 0.7, 0.9, 1.0]
    h2 = Hist2D((rand(10^2), rand(10^2)); binedges=(bins, bins))
    @test h2 == rebin(h2, 1) == (h2 |> rebin(1, 1))
    @test integral(h2) == integral(rebin(h2, 2))
    @test sum(h2.sumw2) == sum(rebin(h2, 2).sumw2)
    @test binedges(rebin(h2, 2)) == ([0, 0.7, 1.0], [0, 0.7, 1.0])
    @test [Set([2]), Set([2])] == FHist.valid_rebin_values(h2)

    h2 = Hist2D((rand(10^2), rand(10^2)); binedges=(0:0.1:1, 0:0.5:1))
    @test nbins(rebin(h2, 10, 2)) == (1, 1)
    @test_throws ErrorException rebin(h2, 2, 10)
    @test [Set([5, 2]), Set([2])] == FHist.valid_rebin_values(h2)
end

@testset "Profile" begin
    xy = collect(hcat([[-2.0, 1.5], [-2.0, -3.5], [-2.0, 1.5], [0.0, -2.0], [0.0, -2.0], [0.0, 0.0], [0.0, 2.0], [0.0, 4.0], [2.0, 1.5]]...)')
    h = Hist2D((xy[:, 1], xy[:, 2]); binedges=(-5:2:5, -5:2:5))

    hx = profile(h, :x)
    hy = profile(h, :y)

    @test hx == (h |> profile(:x))
    @test hx.sumw2 == [0.0, 2.6666666666666665, 1.088, 0.0, 0.0]
    @test hy.sumw2 == [0.0, 0.0, 0.0, 0.6875, 0.0]
    @test bincounts(hx) == [0.0, 0.0, 0.4, 2.0, 0.0]
    @test bincounts(hy) == [-2.0, 0.0, 0.0, -0.5, 0.0]
    @test binedges(hx) == -5:2:5
    @test binedges(hy) == -5:2:5
end

@testset "Projection" begin
    xs = rand(10)
    ys = rand(10)
    zs = rand(10)
    r = 0:0.1:1
    h1x = Hist1D(xs; binedges=r)
    h1y = Hist1D(ys; binedges=r)
    h2 = Hist2D((xs, ys); binedges=(r, r))
    @test project(h2, :x) == (h2 |> project(:x))
    @test h1x == project(h2, :x)
    @test h1y == project(h2, :y)
    h3 = Hist3D((xs, ys, zs); binedges=(r, r, r))
    @test (h3 |> project(:z)) == h2
    @test (h3 |> project(:y) |> project(:x)) == (h2 |> project(:x))
end

@testset "Transpose" begin
    h1 = Hist2D((randn(10), randn(10)); binedges=(-3:3, -3:3))
    t = FHist.transpose
    @test (t ∘ t)(h1) == h1
end


@testset "Repr" begin
    for h1 in (Hist1D(randn(100); binedges=-3:3),
        Hist2D((randn(10), randn(10)); binedges=(-3:3, -3:3)),
        Hist3D((randn(10), randn(10), randn(10)); binedges=(-3:3, -3:3, -3:3)),
        profile(Hist2D((randn(10^5), randn(10^5)); binedges=(-3:0.1:3, -5:0.1:5)), :x))
        if h1 isa Hist3D
            @test all(occursin.(["Hist3D", "edges=", "integral="], repr(h1)))
        else
            @test all(occursin.(["edges:", "total count:", "bin counts:"], repr(h1)))
        end
        @test !occursin("<svg", repr(h1))
        @test all(occursin.(["edges:", "total count:", "bin counts:", "<svg"], repr("text/html", h1)))
        @test all(occursin.(["edges=", "integral=", "Hist"], repr(h1, context=:compact => true)))
    end
end

@testset "Restrict" begin
    h = Hist1D(randn(500); binedges=-5:0.2:5)
    hleft = restrict(h, -Inf, 0.0)
    hright = restrict(h, 0.0, Inf)

    @test h == restrict(h)
    @test nentries(h) == nentries(restrict(h))
    @test restrict(h, -1, 1) == (h |> restrict(-1, 1))
    @test integral(hleft) + integral(hright) == integral(h)
    @test nbins(hleft) + nbins(hright) == nbins(h)
    @test sum(hleft.sumw2) + sum(hright.sumw2) == sum(h.sumw2)

    @test all(-1 .<= bincenters(restrict(h, -1, 1)) .<= 1)
    @test_throws AssertionError restrict(h, 10, Inf)
end

# https://github.com/Moelf/FHist.jl/issues/81
@testset "0-allocation" begin
    data = randn(100)
    h = Hist1D(; binedges=-10:1:10)

    # warm up
    for d in data
        push!(h, d)
    end

    aloc = @allocated for d in data
        push!(h, d)
    end

    @test aloc == 0
end

@testset "Overflow" begin
    a = randn(10^5)
    b = randn(10^5)
    c = randn(10^5)
    ϵ = 1e-6
    edges = -3:0.5:3
    aclip = clamp.(a, nextfloat(first(edges)), prevfloat(last(edges)))
    bclip = clamp.(b, nextfloat(first(edges)), prevfloat(last(edges)))
    cclip = clamp.(c, nextfloat(first(edges)), prevfloat(last(edges)))
    sh1 = fit(Histogram, aclip, edges)
    h1 = Hist1D(a; binedges=edges, overflow=true)
    @test convert(Histogram, h1) == sh1
    @test h1.sumw2 == bincounts(h1)
    # normalize has approximation built-in
    @test normalize(h1; width=true) != normalize(h1; width=false)

    h1 = Hist1D(a; overflow=true)
    before = last(bincounts(h1))
    push!(h1, nextfloat(last(binedges(h1))))
    after = last(bincounts(h1))
    @test after == before + 1

    sh2 = fit(Histogram, (aclip, bclip), (edges, edges))
    h2 = Hist2D((a, b); binedges=(edges, edges), overflow=true)
    @test convert(Histogram, h2) == sh2
    @test h2.sumw2 == bincounts(h2)
    @test integral(project(h2, :x)) == length(aclip)
    @test rebin(h2, 2).overflow == true

    sh3 = fit(Histogram, (aclip, bclip, cclip), (edges, edges, edges))
    h3 = Hist3D((a, b, c); binedges=(edges, edges, edges), overflow=true)
    @test convert(Histogram, h3) == sh3
    @test h3.sumw2 == bincounts(h3)
    @test project(h3, :x).overflow == h3.overflow
end

@testset "Utils" begin
    @test all(FHist.pearson_err(10) .≈ (3.7015621187164243, 2.7015621187164243))
    @test FHist.sqrt_err(9) == (3, 3)
    @test FHist._is_uniform_bins([1, 2, 3, 4, 5]) == true
    @test FHist._is_uniform_bins([1, 2, 3, 4, 5 + 1e-8]) == false
    @test FHist._is_uniform_bins(1:10) == true
    h = Hist1D(randn(10^3))
    @test FHist.hists_to_bars([h]) == (binedges(h)[1:end-1], bincounts(h), ones(nbins(h)))
end

@testset "2D Restrict" begin
    h = Hist2D((randn(500), randn(500)); binedges=(-5:0.2:5, -5:0.2:5))
    hleftx = restrict(h, -Inf, 0.0)
    hrightx = restrict(h, 0.0, Inf)

    @test h == restrict(h)
    @test nentries(h) == nentries(restrict(h))
    @test restrict(h, -1, 1, -Inf, Inf) == (h |> restrict(-1, 1, -Inf, Inf))
    @test integral(hleftx) + integral(hrightx) == integral(h)
    @test nbins(hleftx)[1] + nbins(hrightx)[1] == nbins(h)[1]
    @test sum(sumw2(hleftx)) + sum(sumw2(hrightx)) == sum(sumw2(h))

    @test all(-1 .<= bincenters(restrict(h, -1, 1) |> project(:x)) .<= 1)
    @test_throws AssertionError restrict(h, 10, Inf)

    hlefty = restrict(h, -Inf, Inf, -Inf, 0.0)
    hrighty = restrict(h, -Inf, Inf, 0.0, Inf)
    @test integral(hlefty) + integral(hrighty) == integral(h)
    @test nbins(hlefty)[2] + nbins(hrighty)[2] == nbins(h)[2]
    @test sum(sumw2(hlefty)) + sum(sumw2(hrighty)) == sum(sumw2(h))
end

@testset "Simple Significance" begin
    h1 = Hist1D(rand(1000); binedges=[0, 1.0])
    h2 = Hist1D(rand(10000); binedges=[0, 1.0])

    @test all(significance(h1, h2) .≈ (9.839916447569484, 0.30998654607114046))
end

include("hdf5.jl")

@testset "Chi2/NDF Weighted Weighted" begin
    d1 = [0.8388269402538583, 0.8958351745183146, 0.2316508164636788, 0.24083718705423285, 0.6245648655679755, 0.35869169913644683, 0.3382060145949445, 0.8061120399619636, 0.5209547814793946, 0.6259242010751945, 0.13123460680118015, 0.35617814386171154, 0.8157429580311278, 0.7786035467788149, 0.4980043281225496, 0.7808633174242755, 0.5768808914502541, 0.48916076104285633, 0.37954007607599927, 0.8996479721776526, 0.643643125252326, 0.9275003933412703, 0.5190487717423468, 0.6740078417541109, 0.9906402647669642, 0.6872471970344323, 0.6332170150890514, 0.7236950610523046, 0.7731161630782892, 0.524178702877647, 0.3091730738196147, 0.6110641705814107, 0.2667196494175682, 0.17812152206392873, 0.4581133428773738, 0.17804225069082424, 0.17344588002118944, 0.9765209912310101, 0.399763986727947, 0.13612451667867131, 0.9623519210879821, 0.47767066031218175, 0.38920310048078, 0.03786397526044594, 0.25890678350921215, 0.5846256514039073, 0.8564875931354698, 0.42745063558744223, 0.0445270703394508, 0.054082027476026195, 0.8354425787850815, 0.6434848261990777, 0.43156307224589163, 0.803162381096874, 0.7963488255945479, 0.45627491532556963, 0.22973217345431007, 0.7923533783398489, 0.44440063296809473, 0.2948068442630375, 0.778228385066929, 0.5204506886078755, 0.2572661880014232, 0.5000750530805148, 0.2119338142281626, 0.464053314692788, 0.9659922356783114, 0.650827098441123, 0.13787816988249268, 0.9716358599471212, 0.7192875673032975, 0.6108497639467901, 0.07524171065183871, 0.8236576677815411, 0.31022589463485906, 0.27074094784709946, 0.1503445762757155, 0.4674221288758059, 0.9942773324640871, 0.9406755509503956, 0.988599843628264, 0.1790449658223855, 0.8268698960761557, 0.7558740264440693, 0.8093830277002685, 0.4109747044794374, 0.4649605562169661, 0.7890482404663629, 0.5707208785678005, 0.5227638801081067, 0.5876098657907026, 0.900340840906362, 0.21252447657755724, 0.8918261492365075, 0.052363365131436135, 0.2963595430945126, 0.06949244066351312, 0.6342954226390055, 0.007590286027126747, 0.9939196749931867]
    d2 = [0.4329643029449205, 0.5368359152978799, 0.5151223481131831, 0.5508341031089367, 0.7267171835180046, 0.41991243760664254, 0.9603414573202825, 0.1004671260163057, 0.8150981879987487, 0.09711405181103161, 0.5890105023757488, 0.671592042424567, 0.8317019198463826, 0.1555419687498485, 0.6169089753293562, 0.9395569902982173, 0.5261669475127486, 0.48785621928301093, 0.7245487688377991, 0.20741285002910448, 0.03557626943692116, 0.6463586019733072, 0.9579471703235853, 0.26748780159595087, 0.7784864998384394, 0.422413629754108, 0.21113075405463266, 0.6332454202183904, 0.9476424579365361, 0.9994375547556938, 0.2678189993558624, 0.8703287160200308, 0.5743783411066758, 0.6665083498396543, 0.8197076941621113, 0.3011396826781949, 0.4897232228285663, 0.8135747533449076, 0.8345151391140637, 0.023283948842404145, 0.712298679596815, 0.15079180894034394, 0.7052855628779243, 0.7495656467069672, 0.010865120106437476, 0.16325894416710096, 0.7599013402667105, 0.8301112578052123, 0.5023094377028445, 0.22239437527462236, 0.16863518859223925, 0.7438943952050405, 0.9891435768700272, 0.7788703722738602, 0.14619653767255114, 0.020377086227996388, 0.5264255527185473, 0.05798883630134788, 0.615823195046625, 0.8534188109811262, 0.4446247933218007, 0.545268910527669, 0.2835957118441361, 0.4917290509779251, 0.9406648289603791, 0.979655851838159, 0.9403753601997822, 0.8108246000548059, 0.6812763184913636, 0.11738140600176217, 0.3775827508062394, 0.24643283278439598, 0.8342097292998591, 0.7216352033083174, 0.938186406454954, 0.48482288626461, 0.2582773882299224, 0.6910712673217132, 0.3522468915691054, 0.6548466942622202, 0.05448464596306324, 0.7045699557104349, 0.10174689515958124, 0.7699654453913868, 0.6895079285734036, 0.2792998054885498, 0.19938759899543934, 0.7911257281914083, 0.08685911239432831, 0.3811565467136796, 0.36085873228400145, 0.18596959299944849, 0.22919674563496406, 0.2663089061351366, 0.5838862706508355, 0.14721179600535128, 0.272687753386062, 0.7649251190779611, 0.024904366177311066, 0.39902276006186543]
    
    h1 = Hist1D(d1; binedges=0:0.01:1)
    h2 = Hist1D(d2; binedges=0:0.01:1)

    chi2, ndf = chi2ndf_ww(h1, h2)

    @test chi2 ≈ 81.6 atol = 0.1
    @test ndf == 87
end