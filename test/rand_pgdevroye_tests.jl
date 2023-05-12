#----------------------------------#
# Sampling method tests
#----------------------------------#

# Test that a simulation of 100,000 draws from a PG(2, 1) distribution has the correct mean and variance
function test_rand_pgdevroye_mean_var()
    b = 2
    z = 1.0
    n = 100000
    s = PolyaGammaHybridSamplers.PolyaGammaHybridSampler(b, z, PolyaGammaHybridSamplers.DEVROYE)
    draws = rand(s, n)
    @test mean(draws) ≈ mean(s) atol=0.003
    @test var(draws) ≈ var(s) atol=0.003
end

# Run tests
@testset "Sampling method tests for Devroye" begin
    test_rand_pgdevroye_mean_var()
end

