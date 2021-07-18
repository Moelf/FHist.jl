using FHist, StatsBase
using Test

@testset "The basics" begin
    a = rand(10^3)
    sth1 = fit(Histogram, a)
    h1sqrt = Hist1D(sth1)
    h1p = Hist1D(sth1; error_mode=:pearson)
    @test h1sqrt.hist == sth1
    @test h1p.hist == sth1
    @test h1sqrt.errors_up == sqrt.(sth1.weights)
    @test h1p.errors_up != sqrt.(sth1.weights)
    p_errs = FHist.pearson_err.(sth1.weights)
    pe_up = [x[1] for x in p_errs]
    @test h1p.errors_up == pe_up

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
    update_error!(h2)
    @test h1 == h2
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
