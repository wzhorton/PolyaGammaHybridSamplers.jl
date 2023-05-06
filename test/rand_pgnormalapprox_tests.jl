#----------------------------------#
# Sampling method tests
#----------------------------------#

# Test that the rand_pgnormalapprox function return 0 for b = 0
function test_rand_b0_normalapprox()
    @test PolyaGammaHybridSamplers.rand_pgnormalapprox(0, 1.0, Random.GLOBAL_RNG) ≈ 0
end

# No tests on Base.rand as it is just a wrapper for the rand_ functions
# No test on rand_pghybrid as it is just a wrapper for the other rand_ functions

# Run tests
@testset "Sampling method tests" begin
    test_rand_b0_normalapprox()
end

