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

@testset "Normalize" begin
    a = rand(10^5)
    wgts1 = 2 .* ones(10^5) |> weights
    h1 = Hist1D(a, wgts1)
    @test integral(h1) ≈ sum(wgts1) atol=1e-8
    @test integral(normalize(h1)) ≈ 1 atol=1e-8
end


@testset "Cumulative" begin
    a = rand(10^5)
    wgts1 = 2 .* ones(10^5) |> weights
    h1 = Hist1D(a, wgts1)
    @test argmax(bincounts(cumulative(h1; forward=true))) == nbins(h1)
    @test argmax(bincounts(cumulative(h1; forward=false))) == 1
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

@testset "Broadcasted push" begin
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

end

@testset "Merging" begin
    h1 = Hist1D(randn(100), -3:3)
    h2 = Hist1D(randn(100), -3:3)
    @test merge(h1,h2) == h1+h2
end

@testset "Rebinning" begin
    h1 = Hist1D(rand(10^2), 0:0.1:1)
    @test h1 == rebin(h1, 1)
    @test integral(h1) == integral(rebin(h1, 5))
    @test sum(h1.sumw2) == sum(rebin(h1, 5).sumw2)

    h2 = Hist1D(rand(10^2), [0.0, 0.1, 0.7, 0.9, 1.0])
    @test h2 == rebin(h2, 1)
    @test integral(h2) == integral(rebin(h2, 2))
    @test sum(h2.sumw2) == sum(rebin(h2, 2).sumw2)

    @test rebin(h1, 2) == (h1 |> rebin(2))
end

@testset "Repr" begin
    h1 = Hist1D(randn(100), -3:3)
    @test all(occursin.(["edges:", "total count:", "bin counts:"], repr(h1)))
end

