"""
    PolyaGammaHybridSampler(b::Integer, z::Real, method::AbstractString)

A sampler for the Polya-Gamma distribution using a hybrid of the saddlepoint approximation, the Devroye method, 
and the normal approximation, each of which are discussed in Windle et al. (2014) -- see README.md for details.
The method used is determined by the the `method` parameter, which defaults to `"hybrid"` when not specified.

# Arguments
- `b::Integer`: The shape parameter of the Polya-Gamma distribution. Must be positive.
- `z::Real`: The exponential tilting parameter of the Polya-Gamma distribution.
- `method::PGSamplingMethod : An Enum object specifying the method used to sample from the Polya-Gamma distribution. 
    Must be one of `PGSamplingMethod.hybrid`, `PGSamplingMethod.devroye`, 
    `PGSamplingMethod.saddlepoint`, or `PGSamplingMethod.normalapprox`.

# Returns
- A `PolyaGammaHybridSampler` object which can be sampled using `rand` or `rand!`.
"""

#----------------------------------#
# Define sampler objects
#----------------------------------#

# Define the sampler type
@enum PGSamplingMethod begin
    hybrid
    devroye
    saddlepoint
    normalapprox
end

struct PolyaGammaHybridSampler{T <: Real, N <: Integer} <: Sampleable{Univariate, Continuous}
    b::N
    z::T
    method::PGSamplingMethod

    function PolyaGammaHybridSampler(b::N, z::T, method::PGSamplingMethod) where {T <: Real, N <: Integer}
        if b < zero(b) 
            throw(DomainError("b must be positive"))
        end
        new{T,N}(b, z, method)
    end
end

# Define the outer constuctors 
function PolyaGammaHybridSampler(b::Integer, z::Real)
    PolyaGammaHybridSampler(b, z, PGSamplingMethod.HYBRID)
end

function PolyaGammaHybridSampler(b::Integer, z::Real, method::PGSamplingMethod)
    PolyaGammaHybridSampler(b, z, method)
end


# Define mean and variance
function Distributions.mean(s::PolyaGammaHybridSampler)
    return s.b * inv(2.0*s.z) * tanh(s.z/2.0)
end

function Distributions.var(s::PolyaGammaHybridSampler)
    return s.b * inv(4.0*s.z^2) * (sinh(s.z)-s.z) * sech(s.z/2.0)^2
end

#----------------------------------#
# Define rand hybrid methods
#----------------------------------#

function Base.rand(rng::AbstractRNG, sampler::PolyaGammaHybridSampler)
    if sampler.method == PGSamplingMethod.HYBRID
        return rand_pghybrid(sampler.b, sampler.z, rng)
    elseif sampler.method == PGSamplingMethod.DEVROYE
        return rand_pgdevroye(sampler.b, sampler.z, rng)
    elseif sampler.method == PGSamplingMethod.SADDLEPOINT
        return rand_pgsaddlepoint(sampler.b, sampler.z, rng)
    else # sampler.method = PGSamplingMethod.NORMALAPPROX
        return rand_pgnormalapprox(sampler.b, sampler.z, rng)
    end
end

function rand_pghybrid(b::Integer, z::Real, rng::AbstractRNG)
    if b > 170
        return pgnormalapprox(b, z, rng)
    elseif b > 13
        return pgsaddlepoint(b, z, rng)
    elseif b > 0
        return pgdevroye(b, z, rng)
    else # b = 0 -> degenerate distribution at 0
        return zero(b)
    end
end