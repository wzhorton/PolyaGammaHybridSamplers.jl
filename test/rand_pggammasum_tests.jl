#----------------------------------#
# Sampling method tests
#----------------------------------#

# Test that a simulation of 100,000 draws from a PG(0.5, 1) distribution has the correct mean and variance
function test_rand_pggammasum_mean_var()
    b = 0.5
    z = 1.0
    n = 100000
    s = PolyaGammaHybridSamplers.PolyaGammaHybridSampler(b, z, PolyaGammaHybridSamplers.GAMMASUM)
    draws = rand(s, n)
    @test mean(draws) ≈ mean(s) atol=0.003
    @test var(draws) ≈ var(s) atol=0.003
end

# Test that a simulation of 100,000 draws from a PG(1.5, 1) distribution has the correct mean and variance
function test_rand_devroyeplusgamma_mean_var()
    b = 1.5
    z = 1.0
    n = 100000
    s = PolyaGammaHybridSamplers.PolyaGammaHybridSampler(b, z, PolyaGammaHybridSamplers.DEVROYEPLUSGAMMASUM)
    draws = rand(s, n)
    @test mean(draws) ≈ mean(s) atol=0.003
    @test var(draws) ≈ var(s) atol=0.003
end

# Run tests
@testset "Sampling method tests for devroye plus gamma sum" begin
    test_rand_pggammasum_mean_var()
    test_rand_devroyeplusgamma_mean_var()
end
