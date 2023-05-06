#----------------------------------#
# Sampling method tests
#----------------------------------#

# Test that the rand_pgnormalapprox function return 0 for b = 0
function test_rand_b0_devroye()
    @test PolyaGammaHybridSamplers.rand_pgdevroye(0, 1.0, Random.GLOBAL_RNG) ≈ 0
end

# Test that a simulation of 10,000 draws from a PG(2, 1) distribution has the correct mean and variance
function test_rand_pgdevroye_mean_var()
    b = 2
    z = 1.0
    n = 100000
    s = PolyaGammaHybridSamplers.PolyaGammaHybridSampler(b, z, PolyaGammaHybridSamplers.DEVROYE)
    draws = rand(s, n)
    @test mean(draws) ≈ PolyaGammaHybridSamplers.pg_mean(b, z) atol=0.003
    @test var(draws) ≈ PolyaGammaHybridSamplers.pg_var(b, z) atol=0.003
end

# Run tests
@testset "Sampling method tests for Devroye" begin
    test_rand_b0_devroye()
    test_rand_pgdevroye_mean_var()
end

