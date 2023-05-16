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

function test_pgsaddlepoint_underflow()
    b = 17
    z = -6.128797345912032e16
    s = PolyaGammaHybridSamplers.PolyaGammaHybridSampler(b, z, PolyaGammaHybridSamplers.SADDLEPOINT)
    @test rand(s) ≈ 8.158207422758152e-18
end

# Run tests
@testset "Sampling method tests for Saddlepoint" begin
    test_rand_pgsaddlepoint_mean_var()
    test_pgsaddlepoint_underflow()
end
