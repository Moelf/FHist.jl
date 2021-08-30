using FHist, StatsBase, Statistics
using Test

@testset "The basics" begin
    a = rand(10^3)
    sth1 = fit(Histogram, a)
    h1sqrt = Hist1D(sth1)
    h1p = Hist1D(sth1)
    @test h1sqrt.hist == sth1
    @test h1p.hist == sth1

    h1 = Hist1D(a)
    sth1 = fit(Histogram, a)
    @test h1.hist.weights ≈ sth1.weights
end

@testset "Fast route" begin
    a = rand(10^5)
    r1 = 0:0.1:1
    r2 = collect(r1)
    for r in (r1,r2)
        h1 = Hist1D(a, r)
        sth1 = fit(Histogram, a, r)
        @test all(h1.hist.weights .≈ sth1.weights)
    end
end

@testset "Slow route" begin
    a = rand(10^5)
    r1 = [0, 0.3, 0.4, 0.8, 1]
    h1 = Hist1D(a, r1)
    sth1 = fit(Histogram, a, r1)
    @test all(h1.hist.weights .≈ sth1.weights)
end

@testset "push! route" begin
    a = rand(10^5)
    r1 = [0, 0.3, 0.4, 0.8, 1]
    h1 = Hist1D(a, r1)
    h2 = Hist1D(Int; bins = r1)
    sth1 = fit(Histogram, a, r1)
    for ele in a
        push!(h2, ele)
    end

    @test all(h1.hist.weights .≈ sth1.weights)
    @test all(h1.hist.weights .≈ h2.hist.weights)
end

@testset "Weighted" begin
    a = rand(10^5)
    wgts1 = ones(10^5) |> weights
    wgts = rand(10^5) |> weights
    r1 = 0:0.1:1
    r2 = collect(r1)

    h1 = Hist1D(a)
    h2 = Hist1D(a, wgts1)
    @test h1==h2
    for r in (r1,r2)
        h1 = Hist1D(a, weights(wgts), r)
        sth1 = fit(Histogram, a, weights(wgts), r)
        @test all(h1.hist.weights .≈ sth1.weights)
    end
end

@testset "Hist2D" begin
    x = rand(10)
    y = rand(10)
    h = Hist2D((x,y))
    @test integral(h) == 10

    rx, ry = 0:0.1:1, 0:0.2:1
    wgts = weights(2*ones(length(x)))
    h = Hist2D((x,y), wgts, (rx, ry))
    @test integral(h) == sum(wgts)
    @test nbins(h) == (length(rx)-1, length(ry)-1)

    @test bincenters(Hist2D((x,y), (0:1,0:1))) == ([0.5], [0.5])
    @test bincenters(Hist2D((x,y), wgts, (0:1,0:1))) == ([0.5], [0.5])
    @test nbins(Hist2D((x,y), ([0,0.5,1],[0,0.5,1]))) == (2,2)
    @test nbins(Hist2D((x,y), ([0,0.3,1],[0,0.3,1]))) == (2,2)
    @test nbins(Hist2D((x,y), wgts, ([0,0.5,1],[0,0.5,1]))) == (2,2)
    @test nbins(Hist2D((x,y), wgts, ([0,0.3,1],[0,0.3,1]))) == (2,2)

    @test integral(Hist2D((x,y), wgts, nbins=(5,5))) == sum(wgts)
    @test integral(Hist2D((x,y), nbins=(5,5))) == length(x)

    h1 = Hist2D((x,y), wgts, (rx, ry))
    h2 = Hist2D(Float64; bins=(rx, ry))
    push!.(h2, x, y, wgts)
    @test h1 == h2
end


@testset "Special bins" begin
    # integer values and integer binning
    a = floor.(Int,abs.(randn(10^6)))
    h1 = Hist1D(a, 0:5)
    h2 = Hist1D(fit(Histogram, a, 0:5))
    @test h1 == h2

    # integer values and integer binning
    a = floor.(Int,abs.(randn(10^6)))
    h1 = Hist1D(a, -6:2:6)
    h2 = Hist1D(fit(Histogram, a, -6:2:6))
    @test h1 == h2

    # test floating point rounding when
    # coercing explicit bins into StepRange
    bins = collect(-3:0.1:3.1)
    h = Hist1D(rand(100), bins)
    @test binedges(h) isa AbstractRange
    @test collect(binedges(h)) == bins
end

@testset "Normalize" begin
    a = rand(10^5)
    wgts1 = 2 .* ones(10^5) |> weights
    h1 = Hist1D(a, wgts1)
    @test integral(h1) ≈ sum(wgts1) atol=1e-8
    @test integral(normalize(h1)) ≈ 1 atol=1e-8

    h1 = Hist2D((a,a), wgts1, (0:0.1:1,0:0.1:1))
    @test integral(h1) ≈ sum(wgts1) atol=1e-8
    @test integral(normalize(h1)) ≈ 1 atol=1e-8
end


@testset "Cumulative" begin
    a = rand(10^5)
    wgts1 = 2 .* ones(10^5) |> weights
    h1 = Hist1D(a, wgts1)
    @test argmax(bincounts(cumulative(h1; forward=true))) == nbins(h1)
    @test argmax(bincounts(cumulative(h1; forward=false))) == 1
    h1 = Hist1D([1,2,2,3,3,3],1:4)
    @test bincounts(cumulative(h1, forward=true)) == [1, 3, 6]
    @test bincounts(cumulative(h1, forward=true)) == cumulative(h1, forward=true).sumw2
    @test bincounts(cumulative(h1, forward=false)) == [6, 5, 3]
    @test bincounts(cumulative(h1, forward=false)) == cumulative(h1, forward=false).sumw2
end


@testset "Statistics" begin
    r = -3:0.1:3
    N = 10^5
    h1 = Hist1D(randn(N), r)
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

    h1 = Hist2D((randn(100), randn(100)), (-3:3,-3:3))
    cx, cy = bincenters(h1)
    # (x,y) tuples of bin centers with same shape as counts
    tbc = ((x,y)->(x,y)).(cx, cy')
    f = x->lookup(h1,x...)
    @test f.(tbc) == bincounts(h1)
    @test ismissing(lookup(h1, 10, 10))
end

@testset "Sample" begin
    @test mean(FHist.sample(Hist1D(rand(10^5), 0:0.1:1), n=10^5)) ≈ 0.5 atol=0.1
end

@testset "Empty" begin
    a = rand(10^5)
    r = 0:0.1:1
    h1 = Hist1D(a, r)
    empty!(h1)
    @test maximum(h1.hist.weights) == 0
    @test maximum(h1.sumw2) == 0
    @test h1.hist.edges[1] == r

    h1 = Hist2D((randn(10),randn(10)), (-3:3,-3:3))
    empty!(h1)
    @test maximum(h1.hist.weights) == 0
    @test maximum(h1.sumw2) == 0
end

@testset "Unsafe push" begin
    a = randn(10^5)
    h1 = Hist1D(a, -3:0.2:3)
    h2 = Hist1D(Int; bins=-3:0.2:3)
    h3 = Hist1D(Int; bins=-3:0.2:3)
    for i in a
        push!(h2, i)
        FHist.unsafe_push!(h3, i)
    end
    @test h1 == h2
    @test h1 == h3
end

@testset "Broadcast" begin
    lookup(h::Hist1D, value) = h.hist.weights[FHist._edge_binindex(h.hist.edges[1], value)]

    h1 = Hist1D(Int; bins=0:3)

    push!.(h1, [0.5, 1.5])
    @test lookup.(h1, [0.5,1.5,2.5]) == [1, 1, 0]

    empty!(h1)
    push!.(h1, [0.5, 1.5], 2)
    @test lookup.(h1, [0.5,1.5,2.5]) == [2, 2, 0]
end

@testset "Arithmetic" begin
    @testset "Unweighted regular binning" begin
        h1 = Hist1D([0.5,1.5,1.5,2.5], 0:3)
        h2 = Hist1D([0.5,1.5,2.5,2.5], 0:3)

        h = h1 + h2
        @test h.hist.weights == [2.0, 3.0, 3.0]
        @test h.sumw2 == [2.0, 3.0, 3.0]

        h = h1 - h2
        @test h.hist.weights == [0.0, 1.0, -1.0]
        @test h.sumw2 == [2.0, 3.0, 3.0]

        h = h1 / h2
        @test h.hist.weights == [1.0, 2.0, 0.5]
        @test h.sumw2 == [2.0, 6.0, 0.375]

        h = h1 * 2
        @test h.hist.weights == [2.0, 4.0, 2.0]
        @test h.sumw2 == [4.0, 8.0, 4.0]

        h = h1/(h1+h2*2)
        @test h.hist.weights ≈ [0.333333, 0.5, 0.2] atol=1e-6
        @test h.sumw2 ≈ [0.17284, 0.21875, 0.0544] atol=1e-6
    end

    @testset "Weighted regular binning" begin
        h1 = Hist1D([0.5,1.5,1.5,2.5], weights([3,3,2,1]), 0:3)
        h2 = Hist1D([0.5,1.5,2.5,2.5], weights([3,3,2,1]), 0:3)

        h = h1/(h1+h2*2)
        @test h.hist.weights ≈ [0.333333, 0.454545, 0.142857] atol=1e-6
        @test h.sumw2 ≈ [0.17284, 0.191107, 0.029155] atol=1e-6
    end

    @testset "Weighted irregular binning" begin
        h1 = Hist1D([0.5,1.5,1.5,2.5], weights([3,3,2,1]), [0,1,2,4])
        h2 = Hist1D([0.5,1.5,2.5,2.5], weights([3,3,2,1]), [0,1,2,4])

        h = h1/(h1+h2*2)
        @test h.hist.weights ≈ [0.333333, 0.454545, 0.142857] atol=1e-6
        @test h.sumw2 ≈ [0.17284, 0.191107, 0.029155] atol=1e-6
    end

    @testset "Hist2D" begin
        h1 = Hist2D(([0.5,1.5,1.5,2.5], [0.5,0.5,0.5,0.5]), (0:3,0:1))
        h2 = Hist2D(([0.5,1.5,2.5,2.5], [0.5,0.5,0.5,0.5]), (0:3,0:1))
        h = h1/(h1+h2*2)
        @test vec(h.hist.weights) ≈ [0.333333, 0.5, 0.2] atol=1e-6
        @test vec(h.sumw2) ≈ [0.17284, 0.21875, 0.0544] atol=1e-6
    end

end

@testset "Merging" begin
    h1 = Hist1D(randn(100), -3:3)
    h2 = Hist1D(randn(100), -3:3)
    @test merge(h1,h2) == h1+h2

    h1 = Hist2D((randn(10),randn(10)), (-3:3,-3:3))
    h2 = Hist2D((randn(10),randn(10)), (-3:3,-3:3))
    @test merge(h1,h2) == h1+h2
end

@testset "Rebinning" begin
    h1 = Hist1D(rand(10^2), 0:0.1:1)
    @test h1 == rebin(h1, 1)
    @test integral(h1) == integral(rebin(h1, 5))
    @test sum(h1.sumw2) == sum(rebin(h1, 5).sumw2)
    @test binedges(rebin(h1, 5)) == [0, 0.5, 1.0]

    h2 = Hist1D(rand(10^2), [0.0, 0.1, 0.7, 0.9, 1.0])
    @test h2 == rebin(h2, 1)
    @test integral(h2) == integral(rebin(h2, 2))
    @test sum(h2.sumw2) == sum(rebin(h2, 2).sumw2)
    @test binedges(rebin(h2, 2)) == [0, 0.7, 1.0]

    @test rebin(h1, 2) == (h1 |> rebin(2))

    h1 = Hist2D((rand(10^2),rand(10^2)), (0:0.1:1,0:0.1:1))
    @test h1 == rebin(h1, 1, 1)
    @test integral(h1) == integral(rebin(h1, 5))
    @test sum(h1.sumw2) == sum(rebin(h1, 5).sumw2)
    @test binedges(rebin(h1, 5)) == ([0, 0.5, 1.0], [0, 0.5, 1.0])

    bins = [0.0, 0.1, 0.7, 0.9, 1.0]
    h2 = Hist2D((rand(10^2),rand(10^2)), (bins,bins))
    @test h2 == rebin(h2, 1) == (h2 |> rebin(1, 1))
    @test integral(h2) == integral(rebin(h2, 2))
    @test sum(h2.sumw2) == sum(rebin(h2, 2).sumw2)
    @test binedges(rebin(h2, 2)) == ([0, 0.7, 1.0], [0, 0.7, 1.0])

    h2 = Hist2D((rand(10^2),rand(10^2)), (0:0.1:1,0:0.5:1))
    @test nbins(rebin(h2, 10, 2)) == (1, 1)
    @test_throws AssertionError rebin(h2, 2, 10)
end

@testset "Profile" begin
    xy = collect(hcat([[-2.0, 1.5], [-2.0, -3.5], [-2.0, 1.5], [0.0, -2.0], [0.0, -2.0], [0.0, 0.0], [0.0, 2.0], [0.0, 4.0], [2.0, 1.5]]...)')
    h = Hist2D((xy[:,1],xy[:,2]), (-5:2:5,-5:2:5))

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
    r = 0:0.1:1
    h1x = Hist1D(xs, r)
    h1y = Hist1D(ys, r)
    h2 = Hist2D((xs,ys), (r,r))
    @test project(h2, :x) == (h2 |> project(:x))
    @test h1x == project(h2, :x)
    @test h1y == project(h2, :y)
end

@testset "Transpose" begin
    h1 = Hist2D((randn(10),randn(10)), (-3:3,-3:3))
    t = FHist.transpose
    @test (t∘t)(h1) == h1
end


@testset "Repr" begin
    for h1 in (Hist1D(randn(100), -3:3),
              Hist2D((randn(10),randn(10)), (-3:3,-3:3)))
        @test all(occursin.(["edges:", "total count:", "bin counts:"], repr(h1)))
        @test !occursin("<svg", repr(h1))
        @test all(occursin.(["edges:", "total count:", "bin counts:", "<svg"], repr("text/html", h1)))
    end
end

@testset "Restrict" begin
    h = Hist1D(randn(500), -5:0.2:5)
    hleft = restrict(h, -Inf, 0.0)
    hright = restrict(h, 0.0, Inf)

    @test h == restrict(h)
    @test restrict(h, -1, 1) == (h |> restrict(-1,1))
    @test integral(hleft) + integral(hright) == integral(h)
    @test nbins(hleft) + nbins(hright) == nbins(h)
    @test sum(hleft.sumw2) + sum(hright.sumw2) == sum(h.sumw2)

    @test all(-1 .<= bincenters(restrict(h,-1,1)) .<= 1)
    @test_throws AssertionError restrict(h, 10, Inf)
end

