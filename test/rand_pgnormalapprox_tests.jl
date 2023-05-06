#----------------------------------#
# Sampling method tests
#----------------------------------#

# Test that the rand_pgnormalapprox function return 0 for b = 0
function test_rand_b0_normalapprox()
    @test PolyaGammaHybridSamplers.rand_pgnormalapprox(0, 1.0, Random.GLOBAL_RNG) ≈ 0
end

# Run tests
@testset "Sampling method tests for normal approx" begin
    test_rand_b0_normalapprox()
end

