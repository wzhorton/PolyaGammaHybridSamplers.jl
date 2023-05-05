"""
    PolyaGammaHybridSampler(b::Integer, z::Real, method::AbstractString)

A sampler for the Polya-Gamma distribution using a hybrid of the saddlepoint approximation, the Devroye method, 
and the normal approximation, each of which are discussed in Windle et al. (2014) -- see README.md for details.
The method used is determined by the the `method` parameter, which defaults to `"hybrid"` when not specified.

# Arguments
- `b::Integer`: The shape parameter of the Polya-Gamma distribution. Must be positive.
- `z::Real`: The exponential tilting parameter of the Polya-Gamma distribution.
- `method::PGSamplingMethod : An Enum object specifying the method used to sample from the Polya-Gamma distribution. 
    Must be one of `HYBRID`, `DEVROYE`, `SADDLEPOINT`, or `NORMALAPPROX`. Defaults to `HYBRID` when not specified.

# Returns
- A `PolyaGammaHybridSampler` object which can be sampled using `rand` or `rand!`.
"""

#----------------------------------#
# Define sampler objects
#----------------------------------#

# Define the sampler type
@enum PGSamplingMethod HYBRID DEVROYE SADDLEPOINT NORMALAPPROX

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
    PolyaGammaHybridSampler(b, z, HYBRID)
end

# Define mean and variance
function Distributions.mean(s::PolyaGammaHybridSampler)
    if iszero(s.z)
        return zero(s.z)
    else
        return s.b * tanh(s.z/2.0) * tanh(s.z/2.0)
    end
end

function Distributions.var(s::PolyaGammaHybridSampler)
    if iszero(s.z)
        return inv(24.0)
    else
        return s.b * inv(4.0*s.z^3) * (sinh(s.z)-s.z) * sech(s.z/2.0)^2
    end
end

#----------------------------------#
# Define rand hybrid methods
#----------------------------------#

function Base.rand(rng::AbstractRNG, s::PolyaGammaHybridSampler)
    if s.method == HYBRID
        return rand_pghybrid(s.b, s.z, rng)
    elseif s.method == DEVROYE
        return rand_pgdevroye(s.b, s.z, rng)
    elseif s.method == SADDLEPOINT
        return rand_pgsaddlepoint(s.b, s.z, rng)
    else # sampler.method = NORMALAPPROX
        return rand_pgnormalapprox(s.b, s.z, rng)
    end
end

function rand_pghybrid(b::Integer, z::Real, rng::AbstractRNG)
    if b > 170
        return rand_pgnormalapprox(b, z, rng)
    elseif b > 13
        return rand_pgsaddlepoint(b, z, rng)
    elseif b > 0
        return rand_pgdevroye(b, z, rng)
    else # b = 0 -> degenerate distribution at 0
        return zero(b)
    end
end