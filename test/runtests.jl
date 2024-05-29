using FHist, StatsBase, Statistics, HDF5
using Test

@testset "Fast route" begin
    a = rand(10^5)
    r1 = 0:0.1:1
    r2 = collect(r1)
    for r in (r1, r2)
        h1 = Hist1D(a; binedges = r)
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
    h1 = Hist1D(a; binedges = r1)
    sth1 = fit(Histogram, a, r1)
    @test nentries(h1) == 10^5
    @test all(bincounts(h1) .≈ sth1.weights)
end

@testset "push! route" begin
    a = rand(10^5)
    r1 = [0, 0.3, 0.4, 0.8, 1]
    h1 = Hist1D(a; binedges = r1)
    h2 = Hist1D(; counttype = Int, binedges = r1)
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
    h1 = Hist1D(; binedges = r1)
    h2 = Hist1D(; binedges = r1)
    h3 = Hist1D(; binedges = r1)
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
    h2 = Hist1D(a; weights = wgts1)
    @test nentries(h1) == 10^5
    @test nentries(h2) == 10^5
    @test h1 == h2
    for r in (r1,r2)
        h1 = Hist1D(a; weights = wgts, binedges=r)
        sth1 = fit(Histogram, a, weights(wgts), r)
        @test all(bincounts(h1) .≈ sth1.weights)
    end
end

@testset "Hist2D" begin
    x = rand(10)
    y = rand(10)
    h = Hist2D((x,y))
    @test integral(h) == 10
    @test nentries(h) == 10

    rx, ry = 0:0.1:1, 0:0.2:1
    wgts = weights(2*ones(length(x)))
    h = Hist2D((x,y); weights = wgts, binedges = (rx, ry))
    @test nentries(h) == 10
    @test integral(h) == sum(wgts)
    @test nbins(h) == (length(rx)-1, length(ry)-1)

    @test bincenters(Hist2D((x,y); binedges = (0:1,0:1))) == ([0.5], [0.5])
    @test bincenters(Hist2D((x,y); weights = wgts, binedges = (0:1,0:1))) == ([0.5], [0.5])
    @test nbins(Hist2D((x,y); binedges =([0,0.5,1],[0,0.5,1]))) == (2,2)
    @test nbins(Hist2D((x,y); binedges =([0,0.3,1],[0,0.3,1]))) == (2,2)
    @test nbins(Hist2D((x,y); weights = wgts, binedges =([0,0.5,1],[0,0.5,1]))) == (2,2)
    @test nbins(Hist2D((x,y); weights = wgts, binedges =([0,0.3,1],[0,0.3,1]))) == (2,2)

    @test integral(Hist2D((x,y); weights = wgts, nbins=(5,5))) == sum(wgts)
    @test integral(Hist2D((x,y); nbins=(5,5))) == length(x)

    h1 = Hist2D((x,y); weights=wgts, binedges = (rx, ry))
    h2 = Hist2D(; counttype = Float64, binedges=(rx, ry))
    h3 = Hist2D(; counttype = Float64, binedges=(rx, ry))
    push!.(h2, x, y, wgts)
    atomic_push!.(h3, x, y, wgts)
    @test h1 == h2
    @test h1 == h3
end

@testset "Hist3D" begin
    x = rand(10)
    y = rand(10)
    z = rand(10)
    h = Hist3D((x,y,z))
    @test integral(h) == 10
    @test nentries(h) == 10

    rx, ry, rz = 0:0.1:1, 0:0.2:1, 0:0.5:1
    wgts = weights(2*ones(length(x)))
    h = Hist3D((x,y,z); weights = wgts, binedges = (rx, ry, rz))
    @test nentries(h) == 10
    @test integral(h) == sum(wgts)
    @test nbins(h) == (length(rx)-1, length(ry)-1, length(rz)-1)

    @test bincenters(Hist3D((x,y,z); binedges = (0:1,0:1,0:1))) == ([0.5], [0.5], [0.5])
    @test bincenters(Hist3D((x,y,z); weights = wgts, binedges = (0:1,0:1,0:1))) == ([0.5], [0.5], [0.5])
    @test nbins(Hist3D((x,y,z), binedges = ([0,0.5,1],[0,0.5,1],[0,0.5,1]))) == (2,2,2)
    @test nbins(Hist3D((x,y,z), binedges = ([0,0.3,1],[0,0.3,1],[0,0.3,1]))) == (2,2,2)
    @test nbins(Hist3D((x,y,z), weights = wgts, binedges = ([0,0.5,1],[0,0.5,1],[0,0.5,1]))) == (2,2,2)
    @test nbins(Hist3D((x,y,z), weights = wgts, binedges = ([0,0.3,1],[0,0.3,1],[0,0.3,1]))) == (2,2,2)

    @test integral(Hist3D((x,y,z); weights = wgts, nbins=(5,5,5))) == sum(wgts)
    @test integral(Hist3D((x,y,z); nbins=(5,5,5))) == length(x)

    h1 = Hist3D((x,y,z); weights = wgts, binedges = (rx, ry, rz))
    h2 = Hist3D(; counttype = Float64, binedges = (rx, ry, rz))
    h3 = Hist3D(; counttype = Float64, binedges = (rx, ry, rz))
    push!.(h2, x, y, z, wgts)
    atomic_push!.(h3, x, y, z, wgts)
    @test h1 == h2
    @test h1 == h3
end

@testset "Iterable fall back" begin
    # 1D
    a = Iterators.flatten([[1,2], [], [3]])
    b = collect(a)
    r = 0:0.1:1
    @test Hist1D(a; binedges = r) == Hist1D(b; binedges = r)
    
    # 2D
    x1 = Iterators.flatten(rand(10))
    y1 = Iterators.flatten(rand(10))
    h1 = Hist2D((x1,y1); binedges = (r,r))
    x2 = collect(x1)
    y2 = collect(y1)
    h2 = Hist2D((x2,y2); binedges = (r,r))
    @test h1 == h2
end

@testset "Special bins" begin
    # integer values and integer binning
    a = floor.(Int,abs.(randn(10^6)))
    h1 = Hist1D(a; binedges = 0:5)
    h2 = fit(Histogram, a, 0:5)
    @test convert(Histogram, h1) == h2

    # integer values and integer binning
    a = floor.(Int,abs.(randn(10^6)))
    h1 = Hist1D(a; binedges = -6:2:6)
    h2 = fit(Histogram, a, -6:2:6)
    @test convert(Histogram, h1) == h2
end

@testset "Normalize" begin
    a = rand(10^5)
    wgts1 = 2 .* ones(10^5) |> weights
    h1 = Hist1D(a; weights = wgts1)
    @test integral(h1) ≈ sum(wgts1) atol=1e-8
    @test integral(normalize(h1; width=false)) ≈ 1 atol=1e-8
    @test nentries(normalize(h1; width=false)) == length(a)

    h1 = Hist2D((a,a); weights = wgts1, binedges = (0:0.1:1,0:0.1:1))
    @test integral(h1) ≈ sum(wgts1) atol=1e-8
    @test integral(normalize(h1)) ≈ 1 atol=1e-8

    h1 = Hist3D((a,a,a); weights = wgts1, binedges = (0:0.1:1,0:0.1:1,0:0.1:1))
    @test integral(h1) ≈ sum(wgts1) atol=1e-8
    @test integral(normalize(h1)) ≈ 1 atol=1e-8
end

@testset "Normalize with width" begin
    t = Hist1D(; binedges = [0, 1, 2, 4])
    push!(t, 0.5, 0.5);
    push!(t, 0.5, 0.5);
    push!(t, 2.5, 0.5);
    # 0.5 + 0.5 + 0.5 = 1.5

    @test integral(t) == 1.5
    @test integral(t; width=true) == 2

    nt = normalize(t; width=false)
    @test bincounts(nt) == bincounts(t) ./ 1.5
    @test integral(nt) == 1 # self-consistent requirement

    ntw = normalize(t; width = true)
    @test bincounts(ntw) == bincounts(t) ./ 2
    @test integral(ntw; width=true) == 1 #self-consistent
end


@testset "Cumulative" begin
    a = rand(10^5)
    wgts1 = 2 .* ones(10^5) |> weights
    h1 = Hist1D(a; weights = wgts1)
    @test argmax(bincounts(cumulative(h1; forward=true))) == nbins(h1)
    @test argmax(bincounts(cumulative(h1; forward=false))) == 1
    h1 = Hist1D([1,2,2,3,3,3]; binedges = 1:4)
    @test bincounts(cumulative(h1, forward=true)) == [1, 3, 6]
    @test bincounts(cumulative(h1, forward=true)) == sumw2(cumulative(h1, forward=true))
    @test bincounts(cumulative(h1, forward=false)) == [6, 5, 3]
    @test bincounts(cumulative(h1, forward=false)) == sumw2(cumulative(h1, forward=false))
end

@testset "Bin errors" begin
    h1 = Hist1D(randn(100); binedges = -3:3)
    @test h1.sumw2 == binerrors(identity, h1)
    @test sqrt.(h1.sumw2) == binerrors(h1)
    h1 = Hist2D((randn(100), randn(100)); binedges = (-3:3,-3:3))
    @test h1.sumw2 == binerrors(identity, h1)
    @test sqrt.(h1.sumw2) == binerrors(h1)
    h1 = Hist3D((randn(100), randn(100), randn(100)); binedges = (-3:3,-3:3,-3:3))
    @test h1.sumw2 == binerrors(identity, h1)
    @test sqrt.(h1.sumw2) == binerrors(h1)

    h2 = Hist1D(; counttype = Float64, binedges = 0:1)
    push!(h2, 0.5, 0.5)
    push!(h2, 0.5, 0.5)

    @test h2.sumw2 == [0.5]
    @test binerrors(h2) == [sqrt(0.5)]
end

@testset "Statistics" begin
    r = -3:0.1:3
    N = 10^5
    h1 = Hist1D(randn(N); binedges = r)
    @test mean(h1) ≈ 0 atol=0.05
    @test std(h1) ≈ 1 atol=0.1
    @test nbins(h1) == length(r)-1
    @test median(h1) == quantile(h1, 0.5)
    @test quantile(h1, 0.0) == first(bincenters(h1))
    @test quantile(h1, 1.0) == last(bincenters(h1))

end

@testset "Lookup" begin
    h1 = Hist1D(randn(100))
    @test lookup.(Ref(h1), bincenters(h1)) == bincounts(h1)
    @test ismissing(lookup(h1, last(binedges(h1)) + 0.1))
    @test ismissing(lookup(h1, first(binedges(h1)) - 0.1))

    h1 = Hist2D((randn(100), randn(100)); binedges = (-3:3,-3:3))
    cx, cy = bincenters(h1)
    points = [(x,y) for x in cx, y in cy]
    @test map(p->lookup(h1,p...), points) == bincounts(h1)
    @test ismissing(lookup(h1, 10, 10))

    h1 = Hist3D((randn(100), randn(100), randn(100)); binedges = (-3:3,-3:3,-3:3))
    cx, cy, cz = bincenters(h1)
    points = [(x,y,z) for x in cx, y in cy, z in cz]
    @test map(p->lookup(h1,p...), points) == bincounts(h1)
    @test ismissing(lookup(h1, 10, 10, 10))
end

@testset "Sample" begin
    @test mean(FHist.sample(Hist1D(rand(10^5); binedges = 0:0.1:1), n=10^5)) ≈ 0.5 atol=0.1
end

@testset "Empty" begin
    a = rand(10^5)
    r = 0:0.1:1
    h1 = Hist1D(a; binedges = r)
    empty!(h1)
    @test maximum(bincounts(h1)) == 0
    @test maximum(sumw2(h1)) == 0
    @test binedges(h1) == r

    h1 = Hist2D((randn(10),randn(10)); binedges = (-3:3,-3:3))
    empty!(h1)
    @test maximum(bincounts(h1)) == 0
    @test maximum(sumw2(h1)) == 0

    h1 = Hist3D((randn(10),randn(10),randn(10)); binedges = (-3:3,-3:3,-3:3))
    empty!(h1)
    @test maximum(bincounts(h1)) == 0
    @test maximum(sumw2(h1)) == 0
end

@testset "Atomic push" begin
    if get(ENV, "CI", "false") == "true"
        @test Threads.nthreads() > 1
    end
    a = randn(10^5)
    h1 = Hist1D(a; binedges = -3:0.2:3)
    h2 = Hist1D(;counttype = Int, binedges =-3:0.2:3)
    h3 = Hist1D(;counttype = Int, binedges =-3:0.2:3)
    h4 = Hist1D(;counttype = Int, binedges =-3:0.2:3)
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
    h1 = Hist1D(; counttype = Int, binedges = 0:3)

    push!.(h1, [0.5, 1.5])
    @test lookup.(h1, [0.5,1.5,2.5]) == [1, 1, 0]

    empty!(h1)
    push!.(h1, [0.5, 1.5], 2)
    @test lookup.(h1, [0.5,1.5,2.5]) == [2, 2, 0]
end

@testset "Arithmetic" begin
    @testset "Unweighted regular binning" begin
        h1 = Hist1D([0.5,1.5,1.5,2.5]; binedges = 0:3)
        h2 = Hist1D([0.5,1.5,2.5,2.5]; binedges = 0:3)

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

        h = h1/(h1+h2*2)
        @test bincounts(h) ≈ [0.333333, 0.5, 0.2] atol=1e-6
        @test sumw2(h) ≈ [0.17284, 0.21875, 0.0544] atol=1e-6
    end

    @testset "Weighted regular binning" begin
        h1 = Hist1D([0.5,1.5,1.5,2.5]; weights = [3,3,2,1], binedges = 0:3)
        h2 = Hist1D([0.5,1.5,2.5,2.5]; weights = [3,3,2,1], binedges = 0:3)

        h = h1/(h1+h2*2)
        @test bincounts(h) ≈ [0.333333, 0.454545, 0.142857] atol=1e-6
        @test sumw2(h) ≈ [0.17284, 0.191107, 0.029155] atol=1e-6
    end

    @testset "Weighted irregular binning" begin
        h1 = Hist1D([0.5,1.5,1.5,2.5]; weights = [3,3,2,1], binedges = [0,1,2,4])
        h2 = Hist1D([0.5,1.5,2.5,2.5]; weights = [3,3,2,1], binedges = [0,1,2,4])

        h = h1/(h1+h2*2)
        @test bincounts(h) ≈ [0.333333, 0.454545, 0.142857] atol=1e-6
        @test sumw2(h) ≈ [0.17284, 0.191107, 0.029155] atol=1e-6
    end

    @testset "Hist2D" begin
        h1 = Hist2D(([0.5,1.5,1.5,2.5], [0.5,0.5,0.5,0.5]); binedges = (0:3,0:1))
        h2 = Hist2D(([0.5,1.5,2.5,2.5], [0.5,0.5,0.5,0.5]); binedges = (0:3,0:1))
        h = h1/(h1+h2*2)
        @test vec(bincounts(h)) ≈ [0.333333, 0.5, 0.2] atol=1e-6
        @test vec(sumw2(h)) ≈ [0.17284, 0.21875, 0.0544] atol=1e-6
    end

    @testset "Hist3D" begin
        h1 = Hist3D(([0.5,1.5,1.5,2.5], [0.5,0.5,0.5,0.5], [0.5,0.5,0.5,0.5]); binedges = (0:3,0:1,0:1))
        h2 = Hist3D(([0.5,1.5,2.5,2.5], [0.5,0.5,0.5,0.5], [0.5,0.5,0.5,0.5]); binedges = (0:3,0:1,0:1))
        h = h1/(h1+h2*2)
        @test vec(bincounts(h)) ≈ [0.333333, 0.5, 0.2] atol=1e-6
        @test vec(sumw2(h)) ≈ [0.17284, 0.21875, 0.0544] atol=1e-6
    end

end

@testset "Merging" begin
    h1 = Hist1D(randn(100); binedges = -3:3)
    h2 = Hist1D(randn(100); binedges = -3:3)
    @test merge(h1,h2) == h1+h2
    @test merge(h1,h2,h1) == h1+h2+h1

    h1 = Hist2D((randn(10),randn(10)); binedges = (-3:3,-3:3))
    h2 = Hist2D((randn(10),randn(10)); binedges = (-3:3,-3:3))
    @test merge(h1,h2) == h1+h2

    h1 = Hist3D((randn(10),randn(10),randn(10)); binedges = (-3:3,-3:3,-3:3))
    h2 = Hist3D((randn(10),randn(10),randn(10)); binedges = (-3:3,-3:3,-3:3))
    @test merge(h1,h2) == h1+h2
end

@testset "Rebinning" begin
    h1 = Hist1D(rand(10^2); binedges = 0:0.1:1)
    @test h1 == rebin(h1, 1)
    @test nentries(h1) == nentries(rebin(h1, 1))
    @test integral(h1) == integral(rebin(h1, 5))
    @test sum(h1.sumw2) == sum(rebin(h1, 5).sumw2)
    @test binedges(rebin(h1, 5)) == [0, 0.5, 1.0]
    @test Set([2, 5]) == FHist.valid_rebin_values(h1)

    h2 = Hist1D(rand(10^2); binedges = [0.0, 0.1, 0.7, 0.9, 1.0])
    @test h2 == rebin(h2, 1)
    @test integral(h2) == integral(rebin(h2, 2))
    @test nentries(h2) == nentries(rebin(h2, 1))
    @test sum(h2.sumw2) == sum(rebin(h2, 2).sumw2)
    @test binedges(rebin(h2, 2)) == [0, 0.7, 1.0]
    @test Set([2]) == FHist.valid_rebin_values(h2)

    @test rebin(h1, 2) == (h1 |> rebin(2))

    h1 = Hist2D((rand(10^2),rand(10^2)); binedges = (0:0.1:1,0:0.1:1))
    @test h1 == rebin(h1, 1, 1)
    @test integral(h1) == integral(rebin(h1, 5))
    @test sum(h1.sumw2) == sum(rebin(h1, 5).sumw2)
    @test binedges(rebin(h1, 5)) == ([0, 0.5, 1.0], [0, 0.5, 1.0])
    @test Set([2]) == FHist.valid_rebin_values(h2)

    bins = [0.0, 0.1, 0.7, 0.9, 1.0]
    h2 = Hist2D((rand(10^2),rand(10^2)); binedges = (bins,bins))
    @test h2 == rebin(h2, 1) == (h2 |> rebin(1, 1))
    @test integral(h2) == integral(rebin(h2, 2))
    @test sum(h2.sumw2) == sum(rebin(h2, 2).sumw2)
    @test binedges(rebin(h2, 2)) == ([0, 0.7, 1.0], [0, 0.7, 1.0])
    @test [Set([2]), Set([2])] == FHist.valid_rebin_values(h2)

    h2 = Hist2D((rand(10^2),rand(10^2)); binedges = (0:0.1:1,0:0.5:1))
    @test nbins(rebin(h2, 10, 2)) == (1, 1)
    @test_throws ErrorException rebin(h2, 2, 10)
    @test [Set([5, 2]), Set([2])] == FHist.valid_rebin_values(h2)
end

@testset "Profile" begin
    xy = collect(hcat([[-2.0, 1.5], [-2.0, -3.5], [-2.0, 1.5], [0.0, -2.0], [0.0, -2.0], [0.0, 0.0], [0.0, 2.0], [0.0, 4.0], [2.0, 1.5]]...)')
    h = Hist2D((xy[:,1],xy[:,2]); binedges = (-5:2:5,-5:2:5))

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
    h1x = Hist1D(xs; binedges =  r)
    h1y = Hist1D(ys; binedges =  r)
    h2 = Hist2D((xs,ys); binedges = (r,r))
    @test project(h2, :x) == (h2 |> project(:x))
    @test h1x == project(h2, :x)
    @test h1y == project(h2, :y)
    h3 = Hist3D((xs,ys,zs); binedges = (r,r,r))
    @test (h3 |> project(:z)) == h2
    @test (h3 |> project(:y) |> project(:x)) == (h2 |> project(:x))
end

@testset "Transpose" begin
    h1 = Hist2D((randn(10),randn(10)); binedges = (-3:3,-3:3))
    t = FHist.transpose
    @test (t∘t)(h1) == h1
end


@testset "Repr" begin
    for h1 in (Hist1D(randn(100); binedges = -3:3),
               Hist2D((randn(10),randn(10)); binedges = (-3:3,-3:3)),
               Hist3D((randn(10),randn(10),randn(10)); binedges = (-3:3,-3:3,-3:3)),
               profile(Hist2D((randn(10^5),randn(10^5)); binedges = (-3:0.1:3, -5:0.1:5)), :x))
        if h1 isa Hist3D
            @test all(occursin.(["Hist3D", "edges=", "integral="], repr(h1)))
        else
            @test all(occursin.(["edges:", "total count:", "bin counts:"], repr(h1)))
        end
        @test !occursin("<svg", repr(h1))
        @test all(occursin.(["edges:", "total count:", "bin counts:", "<svg"], repr("text/html", h1)))
        @test all(occursin.(["edges=", "integral=", "Hist"], repr(h1, context=:compact=>true)))
    end
end

@testset "Restrict" begin
    h = Hist1D(randn(500); binedges = -5:0.2:5)
    hleft = restrict(h, -Inf, 0.0)
    hright = restrict(h, 0.0, Inf)

    @test h == restrict(h)
    @test nentries(h) == nentries(restrict(h))
    @test restrict(h, -1, 1) == (h |> restrict(-1,1))
    @test integral(hleft) + integral(hright) == integral(h)
    @test nbins(hleft) + nbins(hright) == nbins(h)
    @test sum(hleft.sumw2) + sum(hright.sumw2) == sum(h.sumw2)

    @test all(-1 .<= bincenters(restrict(h,-1,1)) .<= 1)
    @test_throws AssertionError restrict(h, 10, Inf)
end

# https://github.com/Moelf/FHist.jl/issues/81
@testset "0-allocation" begin
    data = randn(100)
    h = Hist1D(; binedges =-10:1:10)

    # warm up
    for d in data
        push!(h,d)
    end

    aloc = @allocated for d in data
        push!(h,d)
    end

    @test aloc == 0
end

@testset "Overflow" begin
    a = randn(10^5)
    b = randn(10^5)
    c = randn(10^5)
    ϵ = 1e-6
    edges = -3:0.5:3
    aclip = clamp.(a,nextfloat(first(edges)),prevfloat(last(edges)))
    bclip = clamp.(b,nextfloat(first(edges)),prevfloat(last(edges)))
    cclip = clamp.(c,nextfloat(first(edges)),prevfloat(last(edges)))
    sh1 = fit(Histogram, aclip, edges)
    h1 = Hist1D(a; binedges = edges, overflow=true)
    @test convert(Histogram, h1) == sh1
    @test h1.sumw2 == bincounts(h1)
    # normalize has approximation built-in
    @test normalize(h1; width=true) != normalize(h1; width=false)

    h1 = Hist1D(a; overflow=true)
    before = last(bincounts(h1))
    push!(h1, nextfloat(last(binedges(h1))))
    after = last(bincounts(h1))
    @test after == before + 1

    sh2 = fit(Histogram, (aclip,bclip), (edges,edges))
    h2 = Hist2D((a,b); binedges = (edges,edges), overflow=true)
    @test convert(Histogram, h2) == sh2
    @test h2.sumw2 == bincounts(h2)
    @test integral(project(h2, :x)) == length(aclip)
    @test rebin(h2, 2).overflow == true

    sh3 = fit(Histogram, (aclip,bclip,cclip), (edges,edges,edges))
    h3 = Hist3D((a,b,c); binedges = (edges,edges,edges), overflow=true)
    @test convert(Histogram, h3) == sh3
    @test h3.sumw2 == bincounts(h3)
    @test project(h3, :x).overflow == h3.overflow
end

@testset "Utils" begin
    @test all(FHist.pearson_err(10) .≈ (3.7015621187164243, 2.7015621187164243))
    @test FHist.sqrt_err(9) == (3,3)
    @test FHist._is_uniform_bins([1,2,3,4,5]) == true
    @test FHist._is_uniform_bins([1,2,3,4,5+1e-8]) == false
    @test FHist._is_uniform_bins(1:10) == true
    h = Hist1D(randn(10^3))
    @test FHist.hists_to_bars([h]) == (binedges(h)[1:end-1], bincounts(h), ones(nbins(h)))
end

@testset "2D Restrict" begin
    h = Hist2D((randn(500), randn(500)); binedges = (-5:0.2:5, -5:0.2:5))
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
    h1 = Hist1D(rand(1000);  binedges = [0, 1.0])
    h2 = Hist1D(rand(10000); binedges = [0, 1.0]);

    @test all(significance(h1,h2) .≈ (9.839916447569484, 0.30998654607114046))
end

@testset "SafeLog" begin

    FHist.SafeLog.set_safe_log(true)
    @test FHist.SafeLog.SAFE_LOG[] == true
    
    h = Hist1D(bincounts=[0.5,-0.5],binedges=[1,2,3],sumw2=[1.01,1.01])
    bc = copy(bincounts(h))
    FHist.SafeLog._clip_counts!(bc)
    @test all(bc .>= eps())

    bc = copy(bincounts(h))
    err_up = copy(binerrors(h))
    err_dn = copy(binerrors(h))
    FHist.SafeLog._clip_counts!(bc, err_dn, err_up)
    @test all(@. bc >= eps())
    @test all(@. bc - err_dn >= 0)
    @test all( abs.(bc .+ err_up  .- (bincounts(h) .+ binerrors(h))) .<= eps())

end

include("hdf5.jl")
