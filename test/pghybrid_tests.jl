#----------------------------------#
# PolyaGammaHybridSampler struct tests
#----------------------------------#

# Test that the constructor creates the object with the correct fields
function test_constructor()
    s = PolyaGammaHybridSampler(1, 2.0, SADDLEPOINT)
    @test s.b == 1
    @test s.z == 2.0
    @test s.method == SADDLEPOINT
end

# Test that the outer constructor creates the object with the correct fields
function test_outer_constructor()
    s = PolyaGammaHybridSampler(1, 2.0)
    @test s.b == 1
    @test s.z == 2.0
end

# Test that an error is thrown if b is negative but not when it is 0
function test_domain_constructor_error()
    @test_throws DomainError PolyaGammaHybridSampler(-1, 2.0)
    @test PolyaGammaHybridSampler(0, 1.0).b == 0
end

# Run tests
@testset "Struct tests" begin
    test_constructor()
    test_outer_constructor()
    test_domain_constructor_error()
end


#----------------------------------#
# Sampleable interface tests
#----------------------------------#

# Test that the object implements parts of the Sampleable interface
function test_sampleable()
    s = PolyaGammaHybridSampler(1, 2.0)
    # To do: implement mgf and cf
end

# Test that the mean and variance are correct for b=1, z=0 and b=2, z=1
function test_mean_var()
    s = PolyaGammaHybridSampler(1, 0)
    @test mean(s) ≈ 1/4
    @test var(s) ≈ 1/24
    s = PolyaGammaHybridSampler(2, 1)
    @test mean(s) ≈ 0.46211715726000976
    @test var(s) ≈ 0.06889329077704603
    s = PolyaGammaHybridSampler(500, 0)
    @test mean(s) ≈ 125
    @test var(s) ≈ 125/6
end

# Run tests
@testset "Sampleable interface tests" begin
    test_sampleable()
    test_mean_var()
end


#----------------------------------#
# Sampling method tests
#----------------------------------#

# Test that the rand_pghybrid function return 0 for b = 0
function test_rand_b0_hybrid()
    @test rand(PolyaGammaHybridSampler(0, 1.0)) ≈ 0
end

# Test that the correct methods are used for various values of b
function test_rand_methods()
    @test PolyaGammaHybridSampler(0.5, 1.0).method == GAMMASUM
    @test PolyaGammaHybridSampler(2, 1.0).method == DEVROYE
    @test PolyaGammaHybridSampler(2.5, 1.0).method == DEVROYEPLUSGAMMASUM
    @test PolyaGammaHybridSampler(13, 1.0).method == SADDLEPOINT
    @test PolyaGammaHybridSampler(170, 1.0).method == NORMALAPPROX
end

# Run tests
@testset "Sampling method tests for hybrid" begin
    test_rand_b0_hybrid()
    test_rand_methods()
end


