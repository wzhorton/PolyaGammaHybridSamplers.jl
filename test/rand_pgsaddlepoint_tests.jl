#----------------------------------#
# Sampling method tests
#----------------------------------#

# Test that a simulation of 1,000,000 draws from a PG(50, 2) distribution has the correct mean and variance
function test_rand_pgsaddlepoint_mean_var()
    b = 150
    z = 2.0
    n = 1000000
    s = PolyaGammaHybridSamplers.PolyaGammaHybridSampler(b, z, PolyaGammaHybridSamplers.SADDLEPOINT)
    draws = rand(s, n)
    @test mean(draws) ≈ mean(s) rtol=0.005
    @test var(draws) ≈ var(s) rtol=0.005
end

# Run tests
@testset "Sampling method tests for Saddlepoint" begin
    test_rand_pgsaddlepoint_mean_var()
end
