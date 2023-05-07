#----------------------------------#
# Sampling method tests
#----------------------------------#

# Test that the rand_pgsaddlepoint function return 0 for b = 0
function test_rand_b0_saddlepoint()
    @test PolyaGammaHybridSamplers.rand_pgdevroye(0, 1.0, Random.GLOBAL_RNG) ≈ 0
end

# Test that a simulation of 1,000,000 draws from a PG(50, 2) distribution has the correct mean and variance
function test_rand_pgsaddlepoint_mean_var()
    b = 150
    z = 2.0
    n = 1000000
    s = PolyaGammaHybridSamplers.PolyaGammaHybridSampler(b, z, PolyaGammaHybridSamplers.SADDLEPOINT)
    draws = rand(s, n)
    @test mean(draws) ≈ PolyaGammaHybridSamplers.pg_mean(b, z) rtol=0.005
    @test var(draws) ≈ PolyaGammaHybridSamplers.pg_var(b, z) rtol=0.005
end

# Run tests
@testset "Sampling method tests for Saddlepoint" begin
    test_rand_b0_saddlepoint()
    test_rand_pgsaddlepoint_mean_var()
end
