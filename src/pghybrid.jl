"""
    PolyaGammaHybridSampler(b::Real, z::Real, [method::PGSamplingMethod])

A sampler for the Polya-Gamma distribution using a hybrid of the saddlepoint approximation, the Devroye method, 
the normal approximation, and the sum of gammas approximation, each of which are discussed in 
Windle et al. (2014) -- see README.md for details.

# Arguments
- `b::Real`: The shape parameter of the Polya-Gamma distribution. Must be positive.
- `z::Real`: The exponential tilting parameter of the Polya-Gamma distribution.
- `method::PGSamplingMethod : An Enum object specifying the method used to sample from the Polya-Gamma distribution. 
    Must be one of `DEVROYE`, `SADDLEPOINT`, `NORMALAPPROX`, `GAMMASUM`, or `DEVROYEPLUSGAMMASUM`.
    If omitted, the method is chosen automatically based on the value of `b`.

# Returns
- A `PolyaGammaHybridSampler` object which can be sampled using `rand` or `rand!`.

# Examples
```julia
julia> using PolyaGammaHybridSamplers
julia> s = PolyaGammaHybridSampler(1, 1.0)
julia> rand(s)
```

# Notes
- Automatic selection criteria: 
    - `b > 170` -> `NORMALAPPROX`
    - `b >= 13` -> `SADDLEPOINT`
    - `b >= 1 && isinteger(b)` -> `DEVROYE`
    - `b > 1 && !isinteger(b)` -> `DEVROYEPLUSGAMMASUM`
    - `b >= 0` -> `GAMMASUM`
    - `b = 0` -> degenerate distribution at 0
"""
struct PolyaGammaHybridSampler{T <: Real, N <: Real} <: Sampleable{Univariate, Continuous}
    b::N
    z::T
    method::PGSamplingMethod

    function PolyaGammaHybridSampler(b::N, z::T, method::PGSamplingMethod) where {T <: Real, N <: Real}
        if b < zero(b) 
            throw(DomainError("b must be positive"))
        end
        new{T,N}(b, z, method)
    end
end


#----------------------------------#
# Define interface components
#----------------------------------#

# Determine which method to use
function PolyaGammaHybridSampler(b::Real, z::Real)
    if b >= 170
        return PolyaGammaHybridSampler(b, z, NORMALAPPROX)
    elseif b >= 13
        return PolyaGammaHybridSampler(b, z, SADDLEPOINT)
    elseif b >= 1 && isinteger(b)
        return PolyaGammaHybridSampler(Int(b), z, DEVROYE)
    elseif b > 1
        return PolyaGammaHybridSampler(b, z, DEVROYEPLUSGAMMASUM)
    else # 1 > b >= 0
        return PolyaGammaHybridSampler(b, z, GAMMASUM) 
    end
end

# Define mean and variance
function Distributions.mean(s::PolyaGammaHybridSampler)
    if iszero(s.z)
        return s.b / 4.0
    else
        return s.b * 0.5 * inv(s.z) * tanh(0.5*s.z)
    end
end

function Distributions.var(s::PolyaGammaHybridSampler)
    if iszero(s.z)
        return s.b / 24.0
    elseif isinf(sinh(s.z))
        return zero(s.z) # extremely large z -> 0
    else 
        return s.b * 0.25 * inv(s.z^3) * (sinh(s.z)-s.z) * sech(0.5*s.z)^2
    end
end

# Sampling
function Base.rand(rng::AbstractRNG, s::PolyaGammaHybridSampler)
    if iszero(s.b) # b = 0 -> degenerate distribution at 0
        return zero(s.z)
    elseif s.method == NORMALAPPROX
        return rand_pgnormalapprox(s.b, s.z, rng)
    elseif s.method == SADDLEPOINT
        return rand_pgsaddlepoint(s.b, s.z, rng)
    elseif s.method == DEVROYE
        return rand_pgdevroye(s.b, s.z, rng)
    elseif s.method == GAMMASUM
        return rand_pggammasum(s.b, s.z, rng)
    else # s.method == DEVROYEPLUSGAMMASUM
        return rand_pgdevroye(Int(trunc(s.b)), s.z, rng) + rand_pggammasum(s.b - trunc(s.b), s.z, rng)
    end
end