using Test
using FHist
using FHist.Statistics

@testset "Negative content" begin
    h0 = Hist1D([1, 3, 4, 5, 6, 7, 3]; binedges = range(0, 1, 3))
    original_mean = mean(h0)
    h0.bincounts .-= mean(h0.bincounts)  # This makes sum of weights ≈ 0
    @test abs(sum(h0.bincounts)) < 1e-10  # Verify weights sum to approximately zero
    @test isnan(std(h0))   # Std should be NaN with near-zero sum of weights
end
