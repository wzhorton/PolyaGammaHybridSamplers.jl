#----------------------------------#
# Sampling method tests
#----------------------------------#

# Test that a simulation of 1,000,000 draws from a PG(250, 2) distribution has the correct mean and variance
function test_rand_pgnormalapprox_mean_var()
    b = 250
    z = 2.0
    n = 1000000
    s = PolyaGammaHybridSamplers.PolyaGammaHybridSampler(b, z, PolyaGammaHybridSamplers.NORMALAPPROX)
    draws = rand(s, n)
    @test mean(draws) ≈ mean(s) rtol=0.005
    @test var(draws) ≈ var(s) rtol=0.005
end

# Run tests
@testset "Sampling method tests for normal approx" begin
    test_rand_pgnormalapprox_mean_var()
end

